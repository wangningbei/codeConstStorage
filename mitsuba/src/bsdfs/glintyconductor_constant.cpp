/*
This file is part of Mitsuba, a physically based rendering system.

Copyright (c) 2007-2014 by Wenzel Jakob and others.

This is an implementation of Paper "	
Beibei Wang, Miloš Hašan, Nicolas Holzschuch, Ling-Qi Yan. Example-Based Microstructure Rendering with Constant Storage. ACM Transactions on Graphics, 2020"

Mitsuba is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 3
as published by the Free Software Foundation.

Mitsuba is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>

#include <fstream>
#include "microfacet.h"
#include "ior.h"

#include "glintconst/blending.h"

MTS_NAMESPACE_BEGIN

class GlintyConductorConstant : public BSDF {
public:
	GlintyConductorConstant(const Properties &props) : BSDF(props) {
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));

		std::string materialName = props.getString("material", "Cu");

		Spectrum intEta, intK;
		if (boost::to_lower_copy(materialName) == "none") {
			intEta = Spectrum(0.0f);
			intK = Spectrum(1.0f);
		}
		else {
			intEta.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
			intK.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".k.spd")));
		}

		Float extEta = lookupIOR(props, "extEta", "air");

		m_eta = props.getSpectrum("eta", intEta) / extEta;
		m_k = props.getSpectrum("k", intK) / extEta;

		MicrofacetDistribution distr(props);
		m_type = distr.getType();
		m_sampleVisible = distr.getSampleVisible();

		m_alphaU = new ConstantFloatTexture(distr.getAlphaU());
		if (distr.getAlphaU() == distr.getAlphaV())
			m_alphaV = m_alphaU;
		else
			m_alphaV = new ConstantFloatTexture(distr.getAlphaV());

		m_intrinsic = props.getFloat("intrinsic", 0.1f);
		m_tiles = props.getFloat("tiles", 1.0f);
		m_texture_scale = props.getFloat("textureScale", 1.0f);

		m_visTexels = props.getInteger("visTex", 30);

		fs::path filename;
		filename = Thread::getThread()->getFileResolver()->resolve(props.getString("nameNormalMap"));

		Log(EInfo, "Loading normal map \"%s\"", filename.filename().string().c_str());

		if (!fs::exists(filename))
			Log(EError, "Normal map \"%s\" could not be found!", filename.string().c_str());
		fs::path fullPath = filename.parent_path() / filename.filename();

		SLog(EDebug, "%s", fullPath.string());
		m_name_nor_map = fullPath.string();

		m_blendType = props.getInteger("blendType", 0);

		m_exampleScale = props.getFloat("exampleScale", 1.0f);

		m_blender = new BlendCPU(m_name_nor_map, m_blendType, m_exampleScale, m_exam_resU, m_exam_resV);

		m_target_resU = m_exam_resU * m_texture_scale;
		m_target_resV = m_exam_resV * m_texture_scale;

		m_target_resU = props.getInteger("targetWith", m_target_resU);
		m_target_resV = props.getInteger("targetHeight", m_target_resV);

		int res_target_width = props.getInteger("targetWith", m_target_resU);
		int res_target_height = props.getInteger("targetHeight", m_target_resV);

		m_blender->setTargetInform(res_target_width, res_target_height);
		m_blender->setIntrinsic(m_intrinsic);

		printf("Construct the flakes.\n");
		m_blender->buildFlakes();

		SLog(EDebug, "Start to precompute MinMax Table.");

		if(m_blendType == 0)
			m_blender->m_example_normal_Gaussian->precomputeMinMaxTable();
		else
			m_blender->m_example_normal->precomputeMinMaxTable();

		//this is the memory cost
		size_t memoryCost = m_blender->getSize();
		float memCost = (int)memoryCost / (1024.0f * 1024.0f);
		SLog(EDebug, "The memory cost is %f", memCost);

	}

	~GlintyConductorConstant()
	{
		delete m_blender;
	}

	GlintyConductorConstant(Stream *stream, InstanceManager *manager)
		: BSDF(stream, manager) {
		m_type = (MicrofacetDistribution::EType) stream->readUInt();
		m_sampleVisible = stream->readBool();
		m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
		m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_eta = Spectrum(stream);
		m_k = Spectrum(stream);

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

		stream->writeUInt((uint32_t)m_type);
		stream->writeBool(m_sampleVisible);
		manager->serialize(stream, m_alphaU.get());
		manager->serialize(stream, m_alphaV.get());
		manager->serialize(stream, m_specularReflectance.get());
		m_eta.serialize(stream);
		m_k.serialize(stream);
	}

	void configure() {
		unsigned int extraFlags = 0;
		if (m_alphaU != m_alphaV)
			extraFlags |= EAnisotropic;

		if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
			!m_specularReflectance->isConstant())
			extraFlags |= ESpatiallyVarying;

		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);

		m_usesRayDifferentials = true;

		BSDF::configure();

	}

	Vector v_unproject(const Vector2& projectedDirection)
	{
		return Vector(
			projectedDirection.x, projectedDirection.y,
			sqrtf(abs(1.0 - projectedDirection.x * projectedDirection.x - projectedDirection.y * projectedDirection.y)));
	}
	bool is_valid(const Vector2& projection) {
		return projection.x * projection.x + projection.y * projection.y < 1.0;
	}


	/// Helper function: reflect \c wi with respect to a given surface normal
	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	Float normrnd(Float mu, Float sigma, Sampler *sampler) const {
		Float U = math::clamp(sampler->next1D(), 1e-20f, 1.0f - 1e-20f);
		Float V = math::clamp(sampler->next1D(), 1e-20f, 1.0f - 1e-20f);
		Float X = sqrt(-2.0f * log(U)) * cos(2.0f * M_PI * V);
		return X * sigma + mu;
	}

	void constrainPRF(Float &sigma_p1, Float &sigma_p2) const {
		const int MAX_QUERY_RANGE = 16;
		if (sigma_p1 * FLAKE_PIXEL_SIGMAS * m_target_resU > MAX_QUERY_RANGE)
			sigma_p1 = MAX_QUERY_RANGE / FLAKE_PIXEL_SIGMAS / m_target_resU;
		if (sigma_p2 * FLAKE_PIXEL_SIGMAS *m_target_resV > MAX_QUERY_RANGE)
			sigma_p2 = MAX_QUERY_RANGE / FLAKE_PIXEL_SIGMAS / m_target_resV;
	}

	Spectrum evalMicrofacet(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		/* Stop if this component was not requested */
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo + bRec.wi);

		/* Construct the microfacet distribution matching the
		roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);

		/* Evaluate the microfacet normal distribution */
		const Float D = distr.eval(H);
		if (D == 0)
			return Spectrum(0.0f);

		/* Fresnel factor */
		const Spectrum F = fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k) *
			m_specularReflectance->eval(bRec.its);

		/* Smith's shadow-masking function */
		const Float G = distr.G(bRec.wi, bRec.wo, H);

		/* Calculate the total amount of reflection */
		Float model = D * G / (4.0f * Frame::cosTheta(bRec.wi));

		return F * model;
	}

	Float pdfMicrofacet(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo + bRec.wi);

		/* Construct the microfacet distribution matching the
		roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);

		if (m_sampleVisible)
			return distr.eval(H) * distr.smithG1(bRec.wi, H)
			/ (4.0f * Frame::cosTheta(bRec.wi));
		else
			return distr.pdf(bRec.wi, H) / (4 * absDot(bRec.wo, H));
	}

	Spectrum sampleMicrofacet(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Construct the microfacet distribution matching the
		roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);

		/* Sample M, the microfacet normal */
		Float pdf;
		Normal m = distr.sample(bRec.wi, sample, pdf);

		if (pdf == 0)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
			m_eta, m_k) * m_specularReflectance->eval(bRec.its);

		Float weight;
		if (m_sampleVisible) {
			weight = distr.smithG1(bRec.wo, m);
		}
		else {
			weight = distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (pdf * Frame::cosTheta(bRec.wi));
		}

		return F * weight;
	}

	Spectrum sampleMicrofacet(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Construct the microfacet distribution matching the
		roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);

		/* Sample M, the microfacet normal */
		Normal m = distr.sample(bRec.wi, sample, pdf);

		if (pdf == 0)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
			m_eta, m_k) * m_specularReflectance->eval(bRec.its);

		Float weight;
		if (m_sampleVisible) {
			weight = distr.smithG1(bRec.wo, m);
		}
		else {
			weight = distr.eval(m) * distr.G(bRec.wi, bRec.wo, m)
				* dot(bRec.wi, m) / (pdf * Frame::cosTheta(bRec.wi));
		}

		/* Jacobian of the half-direction mapping */
		pdf /= 4.0f * dot(bRec.wo, m);

		return F * weight;
	}

	bool intersectAABB(const Vector4f &aa1, const Vector4f &bb1, const Vector4f &aa2, const Vector4f &bb2) const {
		if ((aa1[0] > bb2[0] || bb1[0] < aa2[0]) ||
			(aa1[1] > bb2[1] || bb1[1] < aa2[1]) ||
			(aa1[2] > bb2[2] || bb1[2] < aa2[2]) ||
			(aa1[3] > bb2[3] || bb1[3] < aa2[3]))
			return false;
		return true;
	}
	inline int mod(int x, int y) const{
		return ((x % y) + y) % y;
	}

	Float evalNDF(const Vector h, const BSDFSamplingRecord &bRec) const {

		Float D = 0.0f;

		if (m_blendType == 0)
		{
			bool bound = m_blender->m_gaussianizer->boundByGaussian(Vector2(h.x, h.y));
			if (!bound)
				return D;
		}

		const int resU = m_target_resU;
		const int resV = m_target_resV;

		Float x1 = bRec.its.uv.y * m_tiles;    // In Mitsuba, u and v are along column and row, respectively. So we swap them.
		Float x2 = bRec.its.uv.x * m_tiles;

		Float du = (fabs(bRec.its.dudx) + fabs(bRec.its.dudy));
		Float dv = (fabs(bRec.its.dvdx) + fabs(bRec.its.dvdy));
		Float sigma_p1 = dv / 2.0 * m_tiles;
		Float sigma_p2 = du / 2.0 * m_tiles;

		constrainPRF(sigma_p1, sigma_p2);

		int x1a = (int)floor(x1 - FLAKE_PIXEL_SIGMAS * sigma_p1);
		int x1b = (int)floor(x1 + FLAKE_PIXEL_SIGMAS * sigma_p1);
		int x2a = (int)floor(x2 - FLAKE_PIXEL_SIGMAS * sigma_p2);
		int x2b = (int)floor(x2 + FLAKE_PIXEL_SIGMAS * sigma_p2);

		for (int i = x1a; i <= x1b; i++) {
			for (int j = x2a; j <= x2b; j++) {

				Float x1Prime = x1 - i;
				Float x2Prime = x2 - j;

				Vector4f aa((x1Prime - FLAKE_PIXEL_SIGMAS * sigma_p1) * resU,
					(x2Prime - FLAKE_PIXEL_SIGMAS * sigma_p2) * resV, h[0], h[1]);
				Vector4f bb((x1Prime + FLAKE_PIXEL_SIGMAS * sigma_p1) *resU,
					(x2Prime + FLAKE_PIXEL_SIGMAS * sigma_p2) * resV, h[0], h[1]);

				D += m_blender->queryFlakesBlendingEval(aa, bb,
					Vector2f(x1Prime * resU, x2Prime * resV),
					Vector2f(sigma_p1 * resU, sigma_p2 * resV));

			}
		}

		return D;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (bRec.bounces > 1)
			return evalMicrofacet(bRec, measure);

		//if (!bRec.its.hasUVPartials)
		//	return evalMicrofacet(bRec, measure);

		/* Stop if this component was not requested */
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo + bRec.wi);

		/* Construct the microfacet distribution matching the
		roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);

		/* Evaluate the NDF */
		const Float D = evalNDF(H, bRec);
		if (D == 0)
			return Spectrum(0.0f);

		/* Fresnel factor */
		const Spectrum F = fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k) *
			m_specularReflectance->eval(bRec.its);

		/* Smith's shadow-masking function */
		const Float G = distr.G(bRec.wi, bRec.wo, H);

		/* Calculate the total amount of reflection */
		Float model = D * G / (4.0f * Frame::cosTheta(bRec.wi));

		return F * model;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		//if (!bRec.its.hasUVPartials)
		if (bRec.bounces > 1)
			return pdfMicrofacet(bRec, measure);

		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo + bRec.wi);

		return evalNDF(H, bRec) / (4.0f * fabsf(dot(H, bRec.wo))) * Frame::cosTheta(H);
	}

	Vector sampleNDF(const BSDFSamplingRecord &bRec) const {
		Float x1 = bRec.its.uv.y * m_tiles;    // In Mitsuba, u and v are along column and row, respectively. So we swap them.
		Float x2 = bRec.its.uv.x * m_tiles;

		Float du = (fabs(bRec.its.dudx) + fabs(bRec.its.dudy));
		Float dv = (fabs(bRec.its.dvdx) + fabs(bRec.its.dvdy));
		Float sigma_p1 = dv / 2.0 * m_tiles;
		Float sigma_p2 = du / 2.0 * m_tiles;

		constrainPRF(sigma_p1, sigma_p2);

		Float x1Sample = normrnd(x1, sigma_p1, bRec.sampler);
		Float x2Sample = normrnd(x2, sigma_p2, bRec.sampler);
		x1Sample -= floor(x1Sample);
		x2Sample -= floor(x2Sample);
		Float u1Sample = x1Sample * m_target_resU;
		Float u2Sample = x2Sample * m_target_resV;

		Vector4f aa(u1Sample, u2Sample, -1.0f, -1.0f);
		Vector4f bb(u1Sample, u2Sample, 1.0f, 1.0f);

		vector<LinearFlake> candidateFlakes;
		candidateFlakes.reserve(4);
		m_blender->queryFlakesBlendingSample(aa, bb, candidateFlakes);

		vector<Float> candidateWeights;
		candidateWeights.reserve(candidateFlakes.size());

		Float sumWeight = 0.0f;
		for (int i = 0; i < candidateFlakes.size(); i++) {
			Float weight = g(candidateFlakes[i].u0[0], u1Sample, candidateFlakes[i].shape[0]) *
				g(candidateFlakes[i].u0[1], u2Sample, candidateFlakes[i].shape[1]) *
				candidateFlakes[i].area;

			candidateWeights.push_back(weight);
			sumWeight += weight;
		}
		if (sumWeight == 0.0f) {
			return Vector(0.0f, 0.0f, 0.0f);
		}

		Float randNum = bRec.sampler->next1D();
		LinearFlake selectedFlake;
		for (int i = 0; i < candidateFlakes.size(); i++) {
			randNum -= candidateWeights[i] / sumWeight;
			if (randNum <= 0.0f) {
				selectedFlake = candidateFlakes[i];
				break;
			}
		}
		if (randNum > 0.0f)
		{
			selectedFlake = candidateFlakes.back();			
		}
		Vector2f nSample = selectedFlake.getNormal(Vector2f(u1Sample, u2Sample));

		Float n1Sample = normrnd(nSample[0], m_intrinsic, bRec.sampler);
		Float n2Sample = normrnd(nSample[1], m_intrinsic, bRec.sampler);
		Float n3SampleSqr = 1.0f - n1Sample * n1Sample - n2Sample * n2Sample;

		if (n3SampleSqr < 0.0f) {
			return Vector(0.0f, 0.0f, 0.0f);
		}
		Float n3Sample = sqrtf(n3SampleSqr);

		return Vector(n1Sample, n2Sample, n3Sample);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (bRec.bounces > 1)
			return sampleMicrofacet(bRec, sample);

		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo + bRec.wi);

		/* Construct the microfacet distribution matching the
		roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);

		/* Sample M, the microfacet normal */
		Normal m = sampleNDF(bRec);
		if (m[0] == 0.0f && m[1] == 0.0f && m[2] == 0.0f)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		Float pdf = evalNDF(m, bRec) / (4.0f * fabsf(dot(m, bRec.wo))) * Frame::cosTheta(m);
		if (pdf == 0)
			return Spectrum(0.0f);

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		/* Fresnel factor */
		const Spectrum F = fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k) *
			m_specularReflectance->eval(bRec.its);

		/* Smith's shadow-masking function */
		const Float G = distr.G(bRec.wi, bRec.wo, H);

		/* Calculate the total amount of reflection */
		Float model = G / (4.0f * Frame::cosTheta(bRec.wi));

		return F * model * (4.0f * fabsf(dot(m, bRec.wo))) / Frame::cosTheta(m);
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		//if (!bRec.its.hasUVPartials)
		if (bRec.bounces > 1)
			return sampleMicrofacet(bRec, pdf, sample);

		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo + bRec.wi);

		/* Construct the microfacet distribution matching the
		roughness values at the current surface position. */
		MicrofacetDistribution distr(
			m_type,
			m_alphaU->eval(bRec.its).average(),
			m_alphaV->eval(bRec.its).average(),
			m_sampleVisible
			);

		/* Sample M, the microfacet normal */
		Normal m = sampleNDF(bRec);
		if (m[0] == 0.0f && m[1] == 0.0f && m[2] == 0.0f)
			return Spectrum(0.0f);

		/* Perfect specular reflection based on the microfacet normal */
		bRec.wo = reflect(bRec.wi, m);
		bRec.eta = 1.0f;
		bRec.sampledComponent = 0;
		bRec.sampledType = EGlossyReflection;

		pdf = evalNDF(m, bRec) / (4.0f * fabsf(dot(m, bRec.wo))) * Frame::cosTheta(m);
		if (pdf == 0)
			return Spectrum(0.0f);

		/* Side check */
		if (Frame::cosTheta(bRec.wo) <= 0)
			return Spectrum(0.0f);

		/* Fresnel factor */
		const Spectrum F = fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k) *
			m_specularReflectance->eval(bRec.its);

		/* Smith's shadow-masking function */
		const Float G = distr.G(bRec.wi, bRec.wo, H);

		/* Calculate the total amount of reflection */
		Float model = G / (4.0f * Frame::cosTheta(bRec.wi));

		return F * model * (4.0f * fabsf(dot(m, bRec.wo))) / Frame::cosTheta(m);
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "alpha")
				m_alphaU = m_alphaV = static_cast<Texture *>(child);
			else if (name == "alphaU")
				m_alphaU = static_cast<Texture *>(child);
			else if (name == "alphaV")
				m_alphaV = static_cast<Texture *>(child);
			else if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		}
		else {
			BSDF::addChild(name, child);
		}
	}

	Float getRoughness(const Intersection &its, int component) const {
		return 0.5f * (m_alphaU->eval(its).average()
			+ m_alphaV->eval(its).average());
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "GlintyConductor[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  distribution = " << MicrofacetDistribution::distributionName(m_type) << "," << endl
			<< "  sampleVisible = " << m_sampleVisible << "," << endl
			<< "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
			<< "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  eta = " << m_eta.toString() << "," << endl
			<< "  k = " << m_k.toString() << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:
	MicrofacetDistribution::EType m_type;
	ref<Texture> m_specularReflectance;
	ref<Texture> m_alphaU, m_alphaV;
	bool m_sampleVisible;
	Spectrum m_eta, m_k;
	Float m_intrinsic;
	Float m_tiles;

	int m_visTexels;

	//	ref<Texture> m_normalMap_example;
	int m_exam_resU;
	int m_exam_resV;

	int m_target_resU;
	int m_target_resV;

	int m_texture_scale; 

	string m_name_nor_map;

	BlendCPU *m_blender;

	int m_blendType; //0 for gaussian, 1 for variance, 2 for linear.
	Float m_exampleScale; //to adjust the roughness of the normal map



};

/**
* GLSL port of the glinty conductor shader. This version is much more
* approximate -- it only supports the Ashikhmin-Shirley distribution,
* does everything in RGB, and it uses the Schlick approximation to the
* Fresnel reflectance of conductors. When the roughness is lower than
* \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
* reasonably well in a VPL-based preview.
*/
class GlintyConductorConstantShader : public Shader {
public:
	GlintyConductorConstantShader(Renderer *renderer, const Texture *specularReflectance,
		const Texture *alphaU, const Texture *alphaV, const Spectrum &eta,
		const Spectrum &k) : Shader(renderer, EBSDFShader),
		m_specularReflectance(specularReflectance), m_alphaU(alphaU), m_alphaV(alphaV){
		m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
		m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
		m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

		/* Compute the reflectance at perpendicular incidence */
		m_R0 = fresnelConductorExact(1.0f, eta, k);
	}

	bool isComplete() const {
		return m_specularReflectanceShader.get() != NULL &&
			m_alphaUShader.get() != NULL &&
			m_alphaVShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_specularReflectanceShader.get());
		deps.push_back(m_alphaUShader.get());
		deps.push_back(m_alphaVShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_specularReflectance.get());
		renderer->unregisterShaderForResource(m_alphaU.get());
		renderer->unregisterShaderForResource(m_alphaV.get());
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_R0);
	}

	void generateCode(std::ostringstream &oss,
		const std::string &evalName,
		const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_R0;" << endl
			<< endl
			<< "float " << evalName << "_D(vec3 m, float alphaU, float alphaV) {" << endl
			<< "    float ct = cosTheta(m), ds = 1-ct*ct;" << endl
			<< "    if (ds <= 0.0)" << endl
			<< "        return 0.0f;" << endl
			<< "    alphaU = 2 / (alphaU * alphaU) - 2;" << endl
			<< "    alphaV = 2 / (alphaV * alphaV) - 2;" << endl
			<< "    float exponent = (alphaU*m.x*m.x + alphaV*m.y*m.y)/ds;" << endl
			<< "    return sqrt((alphaU+2) * (alphaV+2)) * 0.15915 * pow(ct, exponent);" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
			<< "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
			<< "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
			<< "        return 0.0;" << endl
			<< "    float nDotM = cosTheta(m);" << endl
			<< "    return min(1.0, min(" << endl
			<< "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
			<< "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_schlick(float ct) {" << endl
			<< "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
			<< "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "   vec3 H = normalize(wi + wo);" << endl
			<< "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
			<< "   float alphaU = max(0.2, " << depNames[1] << "(uv).r);" << endl
			<< "   float alphaV = max(0.2, " << depNames[2] << "(uv).r);" << endl
			<< "   float D = " << evalName << "_D(H, alphaU, alphaV)" << ";" << endl
			<< "   float G = " << evalName << "_G(H, wi, wo);" << endl
			<< "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
			<< "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);" << endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_specularReflectance;
	ref<const Texture> m_alphaU;
	ref<const Texture> m_alphaV;
	ref<Shader> m_specularReflectanceShader;
	ref<Shader> m_alphaUShader;
	ref<Shader> m_alphaVShader;
	Spectrum m_R0;
};

Shader *GlintyConductorConstant::createShader(Renderer *renderer) const {
	return new GlintyConductorConstantShader(renderer,
		m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(GlintyConductorConstantShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(GlintyConductorConstant, false, BSDF)
MTS_EXPORT_PLUGIN(GlintyConductorConstant, "Glinty conductor BRDF");
MTS_NAMESPACE_END
