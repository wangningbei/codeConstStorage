#pragma once
#if !defined(__FLAKE_H_)
#define __FLAKE_H_

#include <vector>
#include <string>

#include "mitsuba/core/platform.h"
#include "mitsuba/core/fwd.h"
#include "mitsuba/core/aabb.h"
#include "mitsuba/core/matrix.h"

#include "normalmap.h"

using namespace std;
using namespace mitsuba;

#define SQRT_2 1.4142135623730951
#define SQRT_PI 1.7724538509055159
#define SQRT_2PI 2.5066282746310002

#define FLAKE_SHAPE_SIGMAS 3.0f
#define FLAKE_NORMAL_SIGMAS 4.0f
#define FLAKE_PIXEL_SIGMAS 3.0f

inline float g(const float x, const float sigma) {
    return exp(-0.5f * x * x / (sigma * sigma)) / (sigma * SQRT_2PI);
}

inline float g(const float x, const float mu, const float sigma) {
    return g(x - mu, sigma);
}

class Flake {
public:
    Vector2f u0;
    Vector2f n0;
    Vector2f shape;         // Standard deviation.
    Vector4f aa, bb;        // Bounding box.
    float area;             // Actual area it covers.
public:
    Flake() {}
    Flake(Vector2f u0, Vector2f n0, Vector2f shape, float area): u0(u0), n0(n0), shape(shape), area(area) {}
    virtual ~Flake() {}

    virtual Vector2f getNormal(const Vector2f &u) const = 0;
    virtual void getBoundingBox(float intrinsicRoughness, Vector4f &aa, Vector4f &bb) const = 0;
	virtual void getTightBoundingBox(float intrinsicRoughness, Vector4f &aa, Vector4f &bb) const = 0;
    virtual float contributionToNdf(Vector2f uQuery, Vector2f sigma_p, Vector2f nQuery, float sigma_r) const = 0;
};



class LinearFlake: public Flake {
public:
    Matrix2x2 J;
private:
	float gaussian(Vector2 u, Vector2 mu, float c, Matrix2x2 invCov) const
	{
		Vector2 dU = u - mu;
		Vector2 s_T_C = Vector2((dU.x * invCov(0, 0) + dU.y * invCov(1, 0)), (dU.x * invCov(0, 1) + dU.y * invCov(1, 1)));
		float t1 = s_T_C.x * dU.x + s_T_C.y * dU.y;
	//	cout << s_T_C.x * dU.x <<  "," << s_T_C.y * dU.y << "," << exp(-0.5*t1) << endl;
		return c * exp(-0.5*t1);
	}
	float gaussian(Vector4 u, Vector4 mu, float c, Matrix4x4 invCov) const
	{
		Vector4 dU = u - mu;
		Vector4 s_T_C = Vector4((dU.x * invCov(0, 0) + dU.y * invCov(1, 0) + dU.z * invCov(2, 0) + dU.w*invCov(3, 0)),
			(dU.x * invCov(0, 1) + dU.y * invCov(1, 1) + dU.z * invCov(2, 1) + dU.w*invCov(3, 1)),
			(dU.x * invCov(0, 2) + dU.y * invCov(1, 2) + dU.z * invCov(2, 2) + dU.w*invCov(3, 2)),
			(dU.x * invCov(0, 3) + dU.y * invCov(1, 3) + dU.z * invCov(2, 3) + dU.w*invCov(3, 3)));

		float t1 = s_T_C.x * dU.x + s_T_C.y * dU.y + s_T_C.z * dU.z + s_T_C.w * dU.w;
		//	cout << s_T_C.x * dU.x <<  "," << s_T_C.y * dU.y << "," << exp(-0.5*t1) << endl;
		return c * exp(-0.5*t1);
	}
    inline double intExp(double a, double b, double c, double d, double e, double f) const {
        double t1 = 4.0 * a * c - b * b;
        if (t1 <= 0.0) return 0.0;
        double t2 = f + (b * d * e - c * d * d - a * e * e) / t1;
        return 2.0 * M_PI / sqrt(t1) * exp(t2);
    }

    inline float intExp(float *coeffs) const {
        return intExp(coeffs[0], coeffs[1], coeffs[2], coeffs[3], coeffs[4], coeffs[5]);
    }

    inline float getC(float sigma) const {
        return 1.0f / (SQRT_2PI * sigma);
    }

    inline void addCoeffs(float k1, float k2, float b, float sigma, float *coeffs) const {
        float denominatorInv = 1.0f / (-2.0f * sigma * sigma);
        coeffs[0] += k1 * k1 * denominatorInv;
        coeffs[1] += 2.0f * k1 * k2 * denominatorInv;
        coeffs[2] += k2 * k2 * denominatorInv;
        coeffs[3] += 2.0f * b * k1 * denominatorInv;
        coeffs[4] += 2.0f * b * k2 * denominatorInv;
        coeffs[5] += b * b * denominatorInv;
    }

public:
    LinearFlake() {}
    LinearFlake(Vector2f u0, Vector2f n0, Vector2f shape, float area, Matrix2x2 J): Flake(u0, n0, shape, area), J(J) {}

    ~LinearFlake() {}

    virtual Vector2f getNormal(const Vector2f &u) const {
        float nx = n0[0] + J(0, 0) * (u[0] - u0[0]) + J(0, 1) * (u[1] - u0[1]);
        float ny = n0[1] + J(1, 0) * (u[0] - u0[0]) + J(1, 1) * (u[1] - u0[1]);
        return Vector2f(nx, ny);
    }

    virtual void getBoundingBox(float intrinsicRoughness, Vector4f &aa, Vector4f &bb) const {
        float nxMin = n0[0] - fabsf(J(0, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] - fabsf(J(0, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];
        float nxMax = n0[0] + fabsf(J(0, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] + fabsf(J(0, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];
        float nyMin = n0[1] - fabsf(J(1, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] - fabsf(J(1, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];
        float nyMax = n0[1] + fabsf(J(1, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] + fabsf(J(1, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];

       aa = Vector4f(u0[0] - FLAKE_SHAPE_SIGMAS * shape[0], u0[1] - FLAKE_SHAPE_SIGMAS * shape[1], nxMin - FLAKE_NORMAL_SIGMAS * intrinsicRoughness, nyMin - FLAKE_NORMAL_SIGMAS * intrinsicRoughness);
       bb = Vector4f(u0[0] + FLAKE_SHAPE_SIGMAS * shape[0], u0[1] + FLAKE_SHAPE_SIGMAS * shape[1], nxMax + FLAKE_NORMAL_SIGMAS * intrinsicRoughness, nyMax + FLAKE_NORMAL_SIGMAS * intrinsicRoughness);
	}

	virtual inline void getNormalBoundingBox(float intrinsicRoughness, Vector2f &aa, Vector2f &bb) const {

		float flakeFac = FLAKE_SHAPE_SIGMAS * shape[0];

		float fac1 = (fabsf(J(0, 0)) + fabsf(J(0, 1))) * flakeFac;
		float nxMin = n0[0] - fac1;
		float nxMax = n0[0] + fac1;
		float fac2 = (fabsf(J(1, 0)) + fabsf(J(1, 1))) * flakeFac;
		float nyMin = n0[1] - fac2;
		float nyMax = n0[1] + fac2;

		float roughFac = FLAKE_NORMAL_SIGMAS * intrinsicRoughness;
		aa = Vector2f(nxMin - roughFac, nyMin - roughFac);
		bb = Vector2f(nxMax + roughFac, nyMax + roughFac);
	}

	virtual void getTightBoundingBox(float intrinsicRoughness, Vector4f &aa, Vector4f &bb) const {
		float scaleJ = 1.0;
		float nxMin = n0[0] - scaleJ * fabsf(J(0, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] - scaleJ * fabsf(J(0, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];
		float nxMax = n0[0] + scaleJ * fabsf(J(0, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] + scaleJ * fabsf(J(0, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];
		float nyMin = n0[1] - scaleJ * fabsf(J(1, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] - scaleJ * fabsf(J(1, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];
		float nyMax = n0[1] + scaleJ * fabsf(J(1, 0)) * FLAKE_SHAPE_SIGMAS * shape[0] + scaleJ * fabsf(J(1, 1)) * FLAKE_SHAPE_SIGMAS * shape[1];

		aa = Vector4f(u0[0], u0[1], nxMin, nyMin);
		bb = Vector4f(u0[0], u0[1], nxMax, nyMax);
	}

    virtual float contributionToNdf(Vector2f uQuery, Vector2f sigma_p, Vector2f nQuery, float sigma_r) const {

        float c = area;
        c *= getC(shape[0]);
        c *= getC(shape[1]);
        c *= getC(sigma_p[0]);
        c *= getC(sigma_p[1]);
        c *= getC(sigma_r);
        c *= getC(sigma_r);
		
        float coeffs[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        addCoeffs(1.0f, 0.0f, 0.0f, shape[0], coeffs);
        addCoeffs(0.0f, 1.0f, 0.0f, shape[1], coeffs);
        addCoeffs(1.0f, 0.0f, u0[0] - uQuery[0], sigma_p[0], coeffs);
        addCoeffs(0.0f, 1.0f, u0[1] - uQuery[1], sigma_p[1], coeffs);
        addCoeffs(J(0, 0), J(0, 1), n0[0] - nQuery[0], sigma_r, coeffs);
        addCoeffs(J(1, 0), J(1, 1), n0[1] - nQuery[1], sigma_r, coeffs);

		return c * intExp(coeffs);

    }
};

class Flakes {

public:

	Flakes() {}

	Flakes(NormalMap *_normalMap, const float _intrinsicRoughness) {

		m_resolutionU = _normalMap->getWidth();
		m_resolutionV = _normalMap->getHeight();
		int numFlakes = m_resolutionU * m_resolutionU;

		float step = 1.0f;
		for (float i = 0.0f; i < m_resolutionU; i += step) {
			for (float j = 0.0f; j < m_resolutionV; j += step) {
				Vector2f u0(i, j);
				Vector2f n0 = _normalMap->getNormal(i * 1.0f / m_resolutionU, j * 1.0f / m_resolutionV);

				const float k = 1.5; //from 1.5
				Vector2f shape(step * k / sqrt(12.0f), step * k / sqrt(12.0f));
				Matrix2x2 J = _normalMap->getJ(i, j); //computeSimpleJ(i, j);// v
				float area = step * step;
				Flake *currentFlake = new LinearFlake(u0, n0, shape, area, J);
				m_flakes.push_back(currentFlake);
			}
		}

		configure(_intrinsicRoughness);
	}


	~Flakes() {
		for (int i = 0; i < m_flakes.size(); i++)
		{
			delete m_flakes[i];
		}
		std::vector<Flake*>().swap(m_flakes);
	}

	Vector4 getFlakeMinMax(int index)const
	{
		const Flake *flake = m_flakes[index];
		return Vector4(flake->aa.z, flake->aa.w, flake->bb.z, flake->bb.w);
	}
	size_t getFlakeSize()
	{
		return m_flakes.size() * sizeof(Flake);
	}

private:

	void configure(float intrinsicRoughness) {

		printf("Building bounding boxes for flakes...");
		for (int i = 0; i < m_flakes.size(); i++) {
			Flake *flake = m_flakes[i];
			m_flakes[i]->getTightBoundingBox(intrinsicRoughness, flake->aa, flake->bb);
		}
		printf("OK!\n");
	}

    vector<Flake*> m_flakes;
    int m_resolutionU, m_resolutionV;
};

#endif