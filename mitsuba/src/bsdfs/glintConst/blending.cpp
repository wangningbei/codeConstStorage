
#include "flake.h"
#include "normalmap.h"
#include "gaussianization.h"
#include "blending.h"

#include "mitsuba/core/matrix.h"
#include "mitsuba/core/fwd.h"
#include "mitsuba/core/random.h"
#include <omp.h>

#include <time.h>
#include <fstream>

using namespace OPENEXR_IMF_NAMESPACE;
using namespace mitsuba;

const int MAX_QUERY_RANGE = 128;

BlendCPU::BlendCPU(string _path, int blendType, float exampleScale, int &width, int &height){
	path = _path;

	//example scale here 
	m_example_normal = new NormalMap(path);
	m_example_normal->scaleNormalMap(exampleScale);

	m_blend_type = blendType;
	// preprocess
	cout << "Start Precomputation." << endl;

	m_gaussianizer = new Gaussianizaition();

	//this is a bit redundency, but to fit the Gaussian structure.
	TextureDataFloat input;
	ReadEXRtoTextureData(path, input, exampleScale);

	m_gaussianizer->setInputData(input);
	m_gaussianizer->Precomputations();

	m_res_example_width = m_example_normal->getWidth();
	m_res_example_height = m_example_normal->getHeight();

	width = m_res_example_width;
	height = m_res_example_height;

	m_res_target_p_width = m_res_example_width / 4;
	m_res_target_p_height = m_res_example_height / 4;

	m_invres_target_p_width = 1.0f / m_res_target_p_width;
	m_invres_target_p_height = 1.0f / m_res_target_p_height;

	m_bit_target_p_width = log2(m_res_target_p_width);
	m_bit_target_p_height = log2(m_res_target_p_height);

	m_example_to_target_patch[0] = Vector2i(m_res_target_p_width, m_res_target_p_height);
	m_example_to_target_patch[1] = Vector2i(m_res_target_p_width, 0);
	m_example_to_target_patch[2] = Vector2i(0, m_res_target_p_height);
	m_example_to_target_patch[3] = Vector2i(0, 0);

	m_example_normal_Gaussian = new NormalMap(m_res_example_width, m_res_example_height);

	//compute the normlmap of m_example_normal_Gaussian
	Array2D<Rgba> &normalData = m_example_normal_Gaussian->m_normal_map;

	for (int u = 0; u < m_res_example_width; u++)
	{
		for (int v = 0; v < m_res_example_height; v++)
		{
			Vector2 G = m_gaussianizer->getG(Vector2((u + 0.5f) / m_res_example_width,
				(v + 0.5f) / m_res_example_height));

			normalData[v][u].r = G.x;
			normalData[v][u].g = G.y;
			normalData[v][u].b = 1.0;
			normalData[v][u].a = 1;
		}
	}

	m_step = 1;
	//this is from Yan et al. 2016.
	m_sigmaH = 1.5*m_step / sqrt(12.0);
}


BlendCPU::~BlendCPU()
{
	delete m_example_normal;
	delete m_target_normal;
	delete m_example_normal_Gaussian;

	delete m_gaussianizer;
	for (int i = 0; i < m_hash.size(); i++)
	{
		vector<Vector2i>().swap(m_hash[i]);
	}
	vector<vector<Vector2i>>().swap(m_hash);
}

Vector2 BlendCPU::hash1(Vector2i p)
{
	Vector2f temp = mul(Vector2(p), Matrix2x2(127.1, 269.5, 311.7, 183.3));
	temp = Vector2f(sin(temp.x), sin(temp.y)) * 43758.5453f;

	return  Vector2(temp.x - floor(temp.x), temp.y - floor(temp.y));
}


void BlendCPU::buildFlakes()
{
	if (m_blend_type == EGaussian)
		m_example_normal_Gaussian->constructFlakes(u_intrinsic);
	else
		m_example_normal->constructFlakes(u_intrinsic);
	
}

void BlendCPU::setTargetInform(const int &_witdth, const int &_height)
{
	m_res_target_width = _witdth;
	m_res_target_height = _height;

	m_inv_target_width = 1.0f / m_res_target_width;
	m_inv_target_height = 1.0f / m_res_target_height;

	int patchCount = 2 * PATCHSIZE;

	//when the target is extramely large, use the hash directly, rather than caching it.
	m_hash.resize(patchCount);
	//m_res_example - m_res_patch_res
	int x1 = m_res_example_width *0.5;
	int y1 = m_res_example_height *0.5;
#pragma omp parallel for
	for (int i = -patchCount/2; i < patchCount/2; i++)
	{
		int iidex = i + patchCount / 2;
		m_hash[i + patchCount/2].resize(patchCount);
		for (int j = -patchCount / 2; j < patchCount / 2; j++)
		{
			Vector2 temp = hash1(Vector2i(i, j)) * m_res_example_width;
			m_hash[iidex][j+patchCount / 2].x = int(temp.x) % x1;
			m_hash[iidex][j+patchCount / 2].y = int(temp.y) % y1;
		}
	}
}


void BlendCPU::generateEntireTargetImage( string outputFile)
{
	m_target_normal = new NormalMap(m_res_target_width, m_res_target_height);
	Array2D<Rgba> &normalData = m_target_normal->m_normal_map;

	for (int u = 0; u < m_res_target_width; u++)
	{
		for (int v = 0; v < m_res_target_height; v++)
		{
			Matrix2x2 J;
			Vector2 color = ByExampleProceduralNoise(Vector2i(u, v), J);
			normalData[v][u].r = color.x;
			normalData[v][u].g = color.y;
			normalData[v][u].b = 1;
			normalData[v][u].a = 1;
		}
	}

	RgbaOutputFile file(outputFile.c_str(), m_res_target_width, m_res_target_height, WRITE_RGBA); // 1
	file.setFrameBuffer(&normalData[0][0], 1, m_res_target_width); // 2
	file.writePixels(m_res_target_height); // 3
	printf("Done blended normal map writing\n");

}


inline Matrix2x2 BlendCPU::computeJExample(float x, float y)  {
	float x_plus_1 = x + delta;
	float x_miunus_1 = x - delta;

	Vector2f NX_plus_1_Y = Vector2f(m_example_normal->p(x_plus_1, y), m_example_normal->q(x_plus_1, y));
	Vector2f NX_minus_1_Y = Vector2f(m_example_normal->p(x_miunus_1, y), m_example_normal->q(x_miunus_1, y));

	float y_plus_1 = y + delta;
	float y_miunus_1 = y - delta;

	Vector2f NX_Y_PLUS_1 = Vector2f(m_example_normal->p(x, y_plus_1), m_example_normal->q(x, y_plus_1));
	Vector2f NX_Y_MINUS_1 = Vector2f(m_example_normal->p(x, y_miunus_1), m_example_normal->q(x, y_miunus_1));

	Vector2f J1 = (NX_plus_1_Y - NX_minus_1_Y) / (2.0f * delta);
	Vector2f J2 = (NX_Y_PLUS_1 - NX_Y_MINUS_1) / (2.0f * delta);

	Matrix2x2 J(J1.x, J2.x, J1.y, J2.y);
	return J;
}


void BlendCPU::divideLargeFootprint(Vector2i uv0, Vector2i uv1,
	std::vector<Vector2i> &uv0List, std::vector<Vector2i> &uv1List)
{

	Vector2i patchIndex0 = Vector2i(floor((float)uv0.x / m_res_target_p_width), floor((float)uv0.y / m_res_target_p_height));
	Vector2i patchIndex1 = Vector2i(floor((float)uv1.x / m_res_target_p_width), floor((float)uv1.y / m_res_target_p_height));

	if (patchIndex0 == patchIndex1)
	{
		uv0List.push_back(uv0);
		uv1List.push_back(uv1);
		return;
	}

	Vector2i patchCount = patchIndex1 - patchIndex0;

	std::vector<float> widthValue(patchCount.x + 2);
	widthValue[0] = uv0.x;  widthValue[patchCount.x + 1] = uv1.x;
	
	for (int i = 1; i <= patchCount.x; i++)
	{
		widthValue[i] = (patchIndex0.x + i) * m_res_target_p_width - 1; //from 0.5
	}
	
	std::vector<float> heightValue(patchCount.y + 2);
	heightValue[0] = uv0.y; heightValue[patchCount.y + 1] = uv1.y;

	for (int i = 1; i <= patchCount.y; i++)
	{
		heightValue[i] = (patchIndex0.y + i) * m_res_target_p_height - 1; //from 0.5
	}
	
	for (int i = 0; i < patchCount.x + 1; i++)
	{
		for (int j = 0; j < patchCount.y + 1; j++)
		{
			Vector2i temp0 = Vector2i(i == 0 ? widthValue[i] : widthValue[i] + 1,
				j == 0 ? heightValue[j] : heightValue[j] + 1);
			Vector2i temp1 = Vector2i(widthValue[i + 1], heightValue[j + 1]);

			if (temp1.x < temp0.x || temp1.y < temp0.y) 
			{
				SLog(EDebug, "false");
				continue;
			}
			uv0List.push_back(temp0);
			uv1List.push_back(temp1);
		}
	}
}


inline void BlendCPU::QuadGrid(const Vector2i &uv, Vector4 &w, Vector2i &vertex1)
{
	vertex1.x = (uv.x >> m_bit_target_p_width);
	vertex1.y = (uv.y >> m_bit_target_p_height);

	float weight1 = (mod1(uv.x, m_res_target_p_width)) * m_invres_target_p_width;
	float weight2 = (mod1(uv.y, m_res_target_p_height)) * m_invres_target_p_height;

	w.x = (1.0f- weight1) * (1-weight2);
	w.y = (1.0f - weight1) * weight2;
	w.z = weight1* (1 - weight2);
	w.w = weight1 * weight2;

}

inline void BlendCPU::QuadWeight(const Vector2i &uv,Vector4 &w)
{
	float weight1 = (mod1(uv.x, m_res_target_p_width)) * m_invres_target_p_width;
	float weight2 = (mod1(uv.y, m_res_target_p_height)) * m_invres_target_p_height;

	w.x = (1.0f - weight1) * (1 - weight2);
	w.y = (1.0f - weight1) * weight2;
	w.z = weight1* (1 - weight2);
	w.w = weight1 * weight2;
}


Vector2 BlendCPU::ByExampleProceduralNoise(const Vector2i &_uv, Matrix2x2 &Joutput)
{
	int u = _uv.x;
	int v = _uv.y;

	float w[4];
	//the ID of target patch
	Vector2i vertex[4];

	vertex[0].x = (u >> m_bit_target_p_width);
	vertex[0].y = (v >> m_bit_target_p_height);
	vertex[1] = vertex[0] + Vector2i(0, 1);
	vertex[2] = vertex[0] + Vector2i(1, 0);
	vertex[3] = vertex[0] + Vector2i(1, 1);

	Vector2i uvInTarget = Vector2i(mod1(u, m_res_target_p_width), mod1(v, m_res_target_p_height));

	float weight1 = uvInTarget.x * m_invres_target_p_width;
	float weight2 = uvInTarget.y * m_invres_target_p_height;

	w[0] = (1.0f - weight1) * (1 - weight2);
	w[1] = (1.0f - weight1) * weight2;
	w[2] = weight1* (1 - weight2);
	w[3] = weight1 * weight2;

	Vector2 G(0.0);
	Matrix2x2 J(0.0);
	for (int i = 0; i < 4; i++)
	{
		Vector2i exampleP = getRandom(vertex[i]);
		Vector2i uv1 = exampleP + m_example_to_target_patch[i] + uvInTarget;
		Vector2 G1 = m_gaussianizer->getGInt(uv1);

		Vector4 G1Prime = m_gaussianizer->getGIntPrime(uv1);
		Matrix2x2 GP = Matrix2x2(G1Prime.x, G1Prime.y, G1Prime.z, G1Prime.w);
		G += G1 * w[i];
		J += GP * w[i];
	}

	// Variance-preserving blending	
	float weight = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
	G -= Vector2(0.5f);
	G /= weight;
	G += Vector2(0.5f);

	J /= weight;

	Vector2 GinvPrime = m_gaussianizer->FromG2OrgPrime(G);
	Joutput = Matrix2x2(J(0, 0) * GinvPrime.x, J(0, 1)* GinvPrime.x,
		J(1, 0) * GinvPrime.y, J(1, 1)* GinvPrime.y);
	//J <-G^{ -1 }’(n) * J
	return  m_gaussianizer->FromG2Org(G);

}

Vector2 BlendCPU::ByExampleProceduralNoiseVariance(const Vector2i &_uv, Matrix2x2 &Joutput)
{
	int u = _uv.x;
	int v = _uv.y;

	float w[4];
	Vector2i vertex[4];

	vertex[0].x = (u >> m_bit_target_p_width);
	vertex[0].y = (v >> m_bit_target_p_height);
	vertex[1] = vertex[0] + Vector2i(0, 1);
	vertex[2] = vertex[0] + Vector2i(1, 0);
	vertex[3] = vertex[0] + Vector2i(1, 1);

	Vector2i uvInTarget = Vector2i(mod1(u, m_res_target_p_width), mod1(v, m_res_target_p_height));

	float weight1 = uvInTarget.x * m_invres_target_p_width;
	float weight2 = uvInTarget.y * m_invres_target_p_height;

	w[0] = (1.0f - weight1) * (1 - weight2);
	w[1] = (1.0f - weight1) * weight2;
	w[2] = weight1* (1 - weight2);
	w[3] = weight1 * weight2;

	Vector2 N(0.0f);
	Matrix2x2 J(0.0f);
	for (int i = 0; i < 4; i++)
	{
		Vector2i exampleP = getRandom(vertex[i]);
		Vector2i uv1 = exampleP + m_example_to_target_patch[i] + uvInTarget;
		Vector2 G1 = getExampleNormal(uv1);
		Matrix2x2 J1 = computeJExample(uv1.x, uv1.y);
		N += G1 * w[i];
		J += J1 * w[i];
	}
	Joutput = J;

	// Variance-preserving blending	
	if (m_blend_type == EVariance)
	{
		//the center is 0, so there is no +/- center
		float weight = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
		N /= weight;
		Joutput /= weight;
	}
	return  N;

}

void BlendCPU::constrainPRF(Float &sigma_p1, Float &sigma_p2) {
	
	if (sigma_p1 * FLAKE_PIXEL_SIGMAS * m_res_target_width > MAX_QUERY_RANGE)
		sigma_p1 = MAX_QUERY_RANGE / FLAKE_PIXEL_SIGMAS / m_res_target_width;
	if (sigma_p2 * FLAKE_PIXEL_SIGMAS * m_res_target_height > MAX_QUERY_RANGE)
		sigma_p2 = MAX_QUERY_RANGE / FLAKE_PIXEL_SIGMAS / m_res_target_height;
}


void BlendCPU::processNodeVariance(const Vector4i &node, const Vector4i &footprint, const Vector2 &halfV,
	const Vector4 &x1, float &contribution, const Vector2i *examplePos)
{
	//if the current node has intersection with the footprint, then perform further process
	if (node.x > footprint.z || node.y > footprint.w || node.z < footprint.x || node.w < footprint.y)
		return;

	bool sameX = (node[0] == node[2]);
	bool inside = sameX ?
		pointQueryVaraince(Vector2i(node[0], node[1]), halfV, examplePos) :
		rangeQueryRMQVariance(node, Vector4(halfV.x,halfV.y,halfV.x,halfV.y), examplePos);

	if (!inside)
		return;

	if (sameX)
	{
		int u = node.x; int  v = node.y;
		Matrix2x2 J;
		Vector2 n0 = ByExampleProceduralNoiseVariance(Vector2i(u, v), J);

		LinearFlake flake = LinearFlake(Vector2(u, v), n0, Vector2(m_sigmaH, m_sigmaH), m_step*m_step, J);

		Vector2f aa, bb;
		flake.getNormalBoundingBox(u_intrinsic, aa, bb);
		if (halfV.x < aa.x || halfV.x > bb.x ||
			halfV.y < aa.y || halfV.y > bb.y)
		{
			return;
		}

		contribution += flake.contributionToNdf(Vector2(x1.x, x1.y), Vector2(x1.z, x1.w), halfV, u_intrinsic);
		return;
	}
	else
	{
		int midX = ceil((node.x + node.z) * 0.5);
		int midY = ceil((node.y + node.w) * 0.5);

		processNodeVariance(Vector4i(node.x, node.y, midX - 1, midY - 1), footprint, halfV, x1, contribution, examplePos); //, file);
		processNodeVariance(Vector4i(node.x, midY, midX - 1, node.w), footprint, halfV, x1, contribution, examplePos); //, file);
		processNodeVariance(Vector4i(midX, node.y, node.z, midY - 1), footprint, halfV, x1, contribution, examplePos); //, file);
		processNodeVariance(Vector4i(midX, midY, node.z, node.w), footprint, halfV, x1, contribution, examplePos); //, file);
		return;
	}
}

void BlendCPU::processNode(const Vector4i &node, const Vector4i &footprint, const Vector2 &halfV,
	const Vector4 &hG, const Vector4 &x1, float &contribution, const Vector2i *examplePos)
{
	//if the current node has intersection with the footprint, then perform further process
	if (node.x > footprint.z || node.y > footprint.w || node.z < footprint.x || node.w < footprint.y)
		return;

	bool sameX = (node[0] == node[2]);
	bool inside = (sameX) ?
		pointQuery(Vector2i(node[0], node[1]), hG, examplePos) :
		rangeQueryRMQ(node, hG, examplePos);


	if (!inside)
		return;

	if (sameX)
	{
		int u = node.x; int  v = node.y;

		Matrix2x2 J;
		Vector2	n0 = ByExampleProceduralNoise(Vector2i(u, v), J);

		LinearFlake flake = LinearFlake(Vector2(u, v), n0,
			Vector2(m_sigmaH, m_sigmaH), m_step*m_step, J);

		Vector2f aa, bb;
		flake.getNormalBoundingBox(u_intrinsic, aa, bb);
		if (halfV.x < aa.x || halfV.x > bb.x ||
			halfV.y < aa.y || halfV.y > bb.y)
		{
			return;
		}

		contribution += flake.contributionToNdf(Vector2(x1.x, x1.y),
			Vector2(x1.z, x1.w), halfV, u_intrinsic);
		return;
	}
	else
	{
		int midX = ceil((node.x + node.z) * 0.5);
		int midY = ceil((node.y + node.w) * 0.5);

		processNode(Vector4i(node.x, node.y, midX - 1, midY - 1), footprint, halfV, hG, x1, contribution, examplePos); //, file);
		processNode(Vector4i(node.x, midY, midX - 1, node.w), footprint, halfV, hG, x1, contribution, examplePos); //, file);
		processNode(Vector4i(midX, node.y, node.z, midY - 1), footprint, halfV, hG, x1, contribution, examplePos); //, file);
		processNode(Vector4i(midX, midY, node.z, node.w), footprint, halfV, hG, x1, contribution, examplePos); //, file);
		return;
	}
}
	
float BlendCPU::queryFlakesBlendingEval(Vector4f uvMin, Vector4f uvMax, Vector2f x1, Vector2f x2)
{
	float result(0.0f);

	if (!(uvMin[2] == uvMax[2] && uvMin[3] == uvMax[3]))
		return result;

	Vector2i iuvMin = Vector2i(floor(uvMin.x), floor(uvMin.y));
	Vector2i iuvMax = Vector2i(ceil(uvMax.x), ceil(uvMax.y));

	std::vector<Vector2i>uv0List, uv1List;
	uv0List.reserve(1);
	uv1List.reserve(1);

	divideLargeFootprint(iuvMin, iuvMax,uv0List, uv1List);

	const float halfFac = u_intrinsic * FLAKE_NORMAL_SIGMAS;
	Vector2 hGMin = m_gaussianizer->getGFromValue(Vector2(uvMin[2] - halfFac, uvMin[3] - halfFac));
	Vector2 hGMax = m_gaussianizer->getGFromValue(Vector2(uvMin[2] + halfFac, uvMin[3] + halfFac));

	for (int i = 0; i < uv0List.size(); i++)
	{
		Vector2i vertex[4];

		int u0x = uv0List[i].x;
		int u0y = uv0List[i].y;
		int u1x = uv1List[i].x;
		int u1y = uv1List[i].y;

		vertex[0].x = (u0x >> m_bit_target_p_width);
		vertex[0].y = (u0y >> m_bit_target_p_height);

		vertex[1] = vertex[0] + Vector2i(0, 1);
		vertex[2] = vertex[0] + Vector2i(1, 0);
		vertex[3] = vertex[0] + Vector2i(1, 1);

		Vector2i exampleP[4];
		for (int j = 0; j < 4; j++)
		{
			exampleP[j] = getRandom(vertex[j]) + m_example_to_target_patch[j];
		}

		//the quad of the target patch		
		int tp0x = vertex[0].x * m_res_target_p_width;
		int tp0y = vertex[0].y * m_res_target_p_height;

		Vector2i patchMaximum = Vector2i(tp0x + m_res_target_p_width - 1, tp0y + m_res_target_p_height - 1);

		//find a proper size for the starting level
		int length = max(u1x - u0x, u1y - u0y) + 1;
		int level = min(m_bit_target_p_width, (int)ceil(log2(length)));
		int gridSize = (1 << level);

		Vector2i tp1 = Vector2i(u0x + gridSize - 1, u0y + gridSize - 1);
		Vector2i tp0 = Vector2i(u0x, u0y);

		//if tp1 extends the patchMaximum, move both the tp0 and tp1 back the offset
		Vector2i offset = tp1 - patchMaximum;
		offset.x = max(offset.x, 0);
		offset.y = max(offset.y, 0);
		tp0 -= offset;
		tp1 -= offset;
		
		Vector4i footprint(u0x, u0y, u1x, u1y);

		if (m_blend_type == EVariance || m_blend_type == ELinear)
			processNodeVariance(
			Vector4i(tp0.x, tp0.y, tp1.x, tp1.y),
			footprint,
			Vector2(uvMin[2], uvMin[3]),
			Vector4f(x1.x, x1.y, x2.x, x2.y), result, exampleP);
		else
			processNode(Vector4i(tp0.x, tp0.y, tp1.x, tp1.y), 
			footprint,
			Vector2(uvMin[2], uvMin[3]),
			Vector4(hGMin.x, hGMin.y, hGMax.x, hGMax.y),
			Vector4f(x1.x, x1.y, x2.x, x2.y), result, exampleP);

	}

	return result;

}
void BlendCPU::queryFlakesBlendingSample(Vector4f uvMin, Vector4f uvMax, vector<LinearFlake> &candidateFlakes)
{
	if (!(uvMin[0] == uvMax[0] && uvMin[1] == uvMax[1]))
		return;

	int u0 = floor(uvMin.x);
	int u1 = ceil(uvMin.x);
	int v0 = floor(uvMin.y);
	int v1 = ceil(uvMin.y);

	for (int u = u0; u <= u1; u++)
	for (int v = v0; v <= v1; v++)
	{
		Matrix2x2 J;
		Vector2 n0 = (m_blend_type == EVariance || m_blend_type == ELinear) ? ByExampleProceduralNoiseVariance(Vector2i(u, v), J) :
			ByExampleProceduralNoise(Vector2i(u, v), J);

		LinearFlake flake(Vector2(u, v), n0, Vector2(m_sigmaH), m_step*m_step, J);

		Vector4f aa, bb;
		flake.getBoundingBox(u_intrinsic, aa, bb);
		if (!intersectAABB(aa, bb, uvMin, uvMax))
		{
			continue;
		}
		candidateFlakes.push_back(flake);
	}
}


float BlendCPU::evalDirNDF(const NormalMap *texture, Vector2i uvQuery, float texels, Vector2f hQuery)
{

	int x1 = uvQuery.x;    // In Mitsuba, u and v are along column and row, respectively. So we swap them.
	int x2 = uvQuery.y;

	int sigma_p1 = texels / 2.0;
	int sigma_p2 = texels / 2.0;

	//constrainPRF(sigma_p1, sigma_p2);

	float D = 0.0f;
	
	Vector4f aa( x1 - FLAKE_PIXEL_SIGMAS * sigma_p1, x2 - FLAKE_PIXEL_SIGMAS * sigma_p2, hQuery[0], hQuery[1]);
	Vector4f bb(x1 + FLAKE_PIXEL_SIGMAS * sigma_p1, x2 + FLAKE_PIXEL_SIGMAS * sigma_p2, hQuery[0], hQuery[1]);

	D += queryFlakesBlendingEval(aa, bb,  Vector2f(x1, x2), Vector2f(sigma_p1, sigma_p2));

	return D;
}


void BlendCPU::visPNDF(const NormalMap *texture, const Vector2i center, const float visTexels)
{
	const int res = 512;
	Array2D<Rgba> ndfData;
	ndfData.resizeErase(res, res);
	const float halfFac = u_intrinsic * FLAKE_NORMAL_SIGMAS;

#pragma omp parallel for
	for (int i = 0; i < res; i++)
	{
		for (int j = 0; j < res; j++) 
		{
			const Vector2 h = Vector2(2.0f * i / res - 1.0, 1.0 - 2.0f * j / res);
			if (!is_valid(h))
			{
				ndfData[j][i].r = 1.0f;
				ndfData[j][i].g = 1.0f;
				ndfData[j][i].b = 1.0f;
				ndfData[j][i].a = 1.0f;
				continue;
			}

			Vector2 halfVmin = h - Vector2(halfFac);
			Vector2 halfVmax = h + Vector2(halfFac);

			bool bound = m_blend_type == EGaussian ? m_gaussianizer->boundByGaussian(Vector2f(h.x, h.y)) : true;
			if (!bound)
			{
				ndfData[j][i].r = 0.0f;
				ndfData[j][i].g = 0.0f;
				ndfData[j][i].b = 0.0f;
				ndfData[j][i].a = 1.0f;
				continue;
			}

			float result = evalDirNDF(texture, center, visTexels, Vector2f(h[0], h[1]));

			ndfData[j][i].r = result;
			ndfData[j][i].g = result;
			ndfData[j][i].b = result;
			ndfData[j][i].a = 1.0f;
		}
	}

	std::string path = "PNDF.exr";

	cout << "Start ouput " << endl;
	RgbaOutputFile file(path.c_str(), res, res, WRITE_RGBA); // 1
	file.setFrameBuffer(&ndfData[0][0], 1, res); // 2
	file.writePixels(res); // 3
	cout << "Output Done " << endl;
}


bool BlendCPU::rangeQueryRMQ(const Vector4i &uvMinMAX, const Vector4 &hG, const Vector2i *examplePoses)
{
	Vector2i uv00 = Vector2i(uvMinMAX.x, uvMinMAX.y);
	Vector2i uv01 = Vector2i(uvMinMAX.x, uvMinMAX.w);
	Vector2i uv10 = Vector2i(uvMinMAX.z, uvMinMAX.y);
	Vector2i uv11 = Vector2i(uvMinMAX.z, uvMinMAX.w);

	Vector4 wT[4];
	QuadWeight(uv00, wT[0]);
	QuadWeight(uv01, wT[1]);
	QuadWeight(uv10, wT[2]);
	QuadWeight(uv11, wT[3]);

	Vector2i uvInTarget = Vector2i(mod1(uv00.x, m_res_target_p_width), mod1(uv00.y, m_res_target_p_height));

	Vector4 GminMaxT[4] = { Vector4(0.0), Vector4(0.0), Vector4(0.0), Vector4(0.0) };
	
	int length = uv11.x - uv00.x;
	int h = log2_32(length + 1);

	const vector<vector<int64_t>> &minMaxTable = m_example_normal_Gaussian->getMinMaxTable(h);
	for (int i = 0; i < 4; i++)
	{
		Vector2i exampleP = examplePoses[i];
		Vector2i uv0 = exampleP + uvInTarget;
		Vector2i uv1 = uv0 + uv11 - uv00;

		int64_t data = minMaxTable[uv0.x][uv0.y];
		Vector4f min_max;
		unpack(data, min_max);

		for (int j = 0; j < 4; j++)
		{
			GminMaxT[j] += min_max*wT[j][i];
		}
	}

	for (int j = 0; j < 4; j++)
	{
		GminMaxT[j] -= Vector4(0.5f);
		GminMaxT[j] /= wT[j].length();
		GminMaxT[j] += Vector4(0.5f);
	}

	float minXFinal = min(min(GminMaxT[0].x, GminMaxT[1].x), min(GminMaxT[2].x, GminMaxT[3].x));
	float minYFinal = min(min(GminMaxT[0].y, GminMaxT[1].y), min(GminMaxT[2].y, GminMaxT[3].y));

	if (hG.z < minXFinal || hG.w < minYFinal)
		return false;


	float maxXFinal = max(max(GminMaxT[0].z, GminMaxT[1].z), max(GminMaxT[2].z, GminMaxT[3].z));
	float maxYFinal = max(max(GminMaxT[0].w, GminMaxT[1].w), max(GminMaxT[2].w, GminMaxT[3].w));

	if (hG.x > maxXFinal || hG.y > maxYFinal)
		return false;

	return true;
}

bool BlendCPU::rangeQueryRMQVariance(const Vector4i &uvMinMAX, const Vector4 &hG, const Vector2i *examplePoses)
{
	Vector2i uv00 = Vector2i(uvMinMAX.x, uvMinMAX.y);
	Vector2i uv01 = Vector2i(uvMinMAX.x, uvMinMAX.w);
	Vector2i uv10 = Vector2i(uvMinMAX.z, uvMinMAX.y);
	Vector2i uv11 = Vector2i(uvMinMAX.z, uvMinMAX.w);

	Vector4 wT[4];
	QuadWeight(uv00, wT[0]);
	QuadWeight(uv01, wT[1]);
	QuadWeight(uv10, wT[2]);
	QuadWeight(uv11, wT[3]);

	Vector2i uvInTarget = Vector2i(mod1(uv00.x, m_res_target_p_width), mod1(uv00.y, m_res_target_p_height));

	Vector4 GminMaxT[4] = { Vector4(0.0), Vector4(0.0), Vector4(0.0), Vector4(0.0) };

	int length = uv11.x - uv00.x;
	int h = log2_32(length + 1);

	const vector<vector<int64_t>> &minMaxTable = m_example_normal->getMinMaxTable(h);
	for (int i = 0; i < 4; i++)
	{
		Vector2i exampleP = examplePoses[i];
		Vector2i uv0 = exampleP + uvInTarget;
		Vector2i uv1 = uv0 + uv11 - uv00;

		int64_t data = minMaxTable[uv0.x][uv0.y];
		Vector4f min_max;
		unpack(data, min_max);

		//inside the four vertex j
		for (int j = 0; j < 4; j++)
		{
			GminMaxT[j] += min_max*wT[j][i];
		}
	}

	if (m_blend_type == EVariance)
	{
		for (int j = 0; j < 4; j++)
		{
			GminMaxT[j] /= wT[j].length();
		}
	}

	float minXFinal = min(min(GminMaxT[0].x, GminMaxT[1].x), min(GminMaxT[2].x, GminMaxT[3].x));
	float minYFinal = min(min(GminMaxT[0].y, GminMaxT[1].y), min(GminMaxT[2].y, GminMaxT[3].y));

	if (hG.z < minXFinal || hG.w < minYFinal)
		return false;


	float maxXFinal = max(max(GminMaxT[0].z, GminMaxT[1].z), max(GminMaxT[2].z, GminMaxT[3].z));
	float maxYFinal = max(max(GminMaxT[0].w, GminMaxT[1].w), max(GminMaxT[2].w, GminMaxT[3].w));

	if (hG.x > maxXFinal || hG.y > maxYFinal)
		return false;

	return true;

}

bool BlendCPU::pointQuery(const Vector2i& _uv, const Vector4& hG, const Vector2i *examplePoses)
{

	float w[4];
	Vector2i uvInTarget = Vector2i(mod1(_uv.x, m_res_target_p_width), mod1(_uv.y, m_res_target_p_height));
	float weight1 = uvInTarget.x * m_invres_target_p_width;
	float weight2 = uvInTarget.y * m_invres_target_p_height;

	w[0] = (1.0f - weight1) * (1 - weight2);
	w[1] = (1.0f - weight1) * weight2;
	w[2] = weight1* (1 - weight2);
	w[3] = weight1 * weight2;

	Vector4 GminMaxT1(0.0f, 0.0f, 0.0f, 0.0f);

	const Flakes * flakes = m_example_normal_Gaussian->m_flakes;

	for (int i = 0; i < 4; i++)
	{
		Vector2i uv1 = examplePoses[i] + uvInTarget;
		Vector4f min_max = flakes->getFlakeMinMax(uv1.x * m_res_example_height + uv1.y);
		GminMaxT1 += min_max*w[i];

	}

	float weight = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);

	GminMaxT1 -= Vector4(0.5f);
	GminMaxT1 /= weight;
	GminMaxT1 += Vector4(0.5f);

	if (hG.z < GminMaxT1.x || hG.w < GminMaxT1.y || hG.x > GminMaxT1.z || hG.y > GminMaxT1.w)
		return false;

	return true;

}

bool BlendCPU::pointQueryVaraince(const Vector2i& _uv, const Vector2& hG, const Vector2i *examplePoses)
{
	float w[4];
	Vector2i uvInTarget = Vector2i(mod1(_uv.x, m_res_target_p_width), mod1(_uv.y, m_res_target_p_height));
	float weight1 = uvInTarget.x * m_invres_target_p_width;
	float weight2 = uvInTarget.y * m_invres_target_p_height;

	w[0] = (1.0f - weight1) * (1 - weight2);
	w[1] = (1.0f - weight1) * weight2;
	w[2] = weight1* (1 - weight2);
	w[3] = weight1 * weight2;

	Vector4 GminMaxT1(0.0f, 0.0f, 0.0f, 0.0f);

	const Flakes * flakes = m_example_normal->m_flakes;

	for (int i = 0; i < 4; i++)
	{
		Vector2i uv1 = examplePoses[i] + uvInTarget;
		Vector4f min_max = flakes->getFlakeMinMax(uv1.x * m_res_example_height + uv1.y); 
		GminMaxT1 += min_max*w[i];
	}

	float weight;

	if (m_blend_type == EVariance)
	{
		weight = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
		GminMaxT1 /= weight;
	}

	if (hG.x < GminMaxT1.x || hG.y < GminMaxT1.y || hG.x > GminMaxT1.z || hG.y > GminMaxT1.w)
		return false;
	return true;

}