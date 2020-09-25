#pragma once
#if !defined(__BLENDCPU_H_)
#define __BLENDCPU_H_

#include "normalmap.h"
#include "gaussianization.h"
#include "flake.h"

const int PATCHSIZE = 1000;


enum BlendType
{
	EGaussian = 0,
	EVariance,
	ELinear,
};

class BlendCPU{

public:
	BlendCPU(string _path, int blendType, float exampleScale, int &width, int &height);
	~BlendCPU();

	size_t getSize()
	{
		if (m_blend_type == 0)
		{
			return m_gaussianizer->getSize() + m_example_normal_Gaussian->getSize() + m_example_normal_Gaussian->m_flakes->getFlakeSize();
		}
		else
			return m_example_normal->getSize() + m_example_normal->m_flakes->getFlakeSize();
		
	}

	void setTargetInform(const int &_witdth, const int &_height);
	void generateEntireTargetImage(string outputFile);

	void setIntrinsic(float intrinsic){ u_intrinsic = intrinsic; };

	void buildFlakes();

	void visPNDF(const NormalMap *texture, const Vector2i center, const float visTexels);

	//Perform point query for Gaussian Preserving Blending case
	bool pointQuery(const Vector2i& _uv, const Vector4& hG, const Vector2i *examplePoses);
	//Perform point query for non-Gaussian Preserving Blending case
	bool pointQueryVaraince(const Vector2i& _uv, const Vector2& hG, const Vector2i *examplePoses);

	//Perform range query with RMQ for Gaussian Preserving Blending case
	bool rangeQueryRMQ(const Vector4i& uvMinMAX, const Vector4 &h, const Vector2i *examplePos);
	//Perform range query with RMQ for non- Gaussian Preserving Blending case
	bool rangeQueryRMQVariance(const Vector4i &uvMinMAX, const Vector4 &hG, const Vector2i *examplePoses);

	const NormalMap *getExampleTexture(){ return m_example_normal; }
	const NormalMap *getExampleGaussianTexture(){ return m_example_normal_Gaussian; }
	const NormalMap *getTargetTexture(){ return m_target_normal; }

	//This is the evaluation and sampling
	float queryFlakesBlendingEval(Vector4f uvMin, Vector4f uvMax, Vector2f x1, Vector2f x2);
	void queryFlakesBlendingSample(Vector4f uvMin, Vector4f uvMax, vector<LinearFlake> &candidateFlakes);

	inline bool intersectAABB(const Vector4f &aa1, const Vector4f &bb1, const Vector4f &aa2, const Vector4f &bb2) {
		if ((aa1[0] > bb2[0] || bb1[0] < aa2[0]) ||
			(aa1[1] > bb2[1] || bb1[1] < aa2[1]) ||
			(aa1[2] > bb2[2] || bb1[2] < aa2[2]) ||
			(aa1[3] > bb2[3] || bb1[3] < aa2[3]))
			return false;
		return true;
	}

	NormalMap *m_example_normal;
	NormalMap *m_target_normal; //this is only for debug
	NormalMap *m_example_normal_Gaussian;
	Gaussianizaition *m_gaussianizer;

private:

	//if the range in the target texture comes from different patches, then subdivide then to small subranges.
	void divideLargeFootprint(Vector2i uv0, Vector2i uv1, std::vector<Vector2i> &uv0List, std::vector<Vector2i> &uv1List);

	inline Vector2i getRandom(const Vector2i& input)
	{
		return m_hash[input.x + PATCHSIZE][input.y + PATCHSIZE];
	}

	/******Helper funtions for RMQ unpacking*****/
	inline int log2_32(uint32_t value)
	{
		value |= value >> 1;
		value |= value >> 2;
		value |= value >> 4;
		value |= value >> 8;
		value |= value >> 16;
		return tab32[(uint32_t)(value * 0x07C4ACDD) >> 27];
	}

	inline void unpack(int64_t value, Vector4 &minmax)
	{
		uint32_t aUnpacked = (value >> 32);
		uint32_t bUnpacked = (value & 0xFFFFFFFF);
		minmax.x = (aUnpacked >> 16);
		minmax.y = uint16_t(aUnpacked & 0xFFFF);
		minmax.z = (bUnpacked >> 16);
		minmax.w = uint16_t(bUnpacked & 0xFFFF);
		minmax /= SCALEPACKING;
		minmax = minmax * 6.0f - Vector4(3.0f);
	}
	/******Helper funtions for RMQ unpacking end*****/


	inline bool inside(const Vector2 &point)
	{
		return (point.x <= 1.0 && point.y <= 1.0 && point.x >= 0.0f && point.y >= 0.0f);
	}

	float evalDirNDF(const NormalMap *texture, Vector2i uvQuery, float texels, Vector2f hQuery);

	//For implicit hiearhcy traversal. Gaussian case.
	void processNode(const Vector4i &node, const Vector4i &footprint, const Vector2 &halfV,
		const Vector4 &hG, const Vector4 &x1, float &contribution, const Vector2i *examplePos);

	//For implicit hiearhcy traversal. non-Gaussian case.
	void processNodeVariance(const Vector4i &node, const Vector4i &footprint, const Vector2 &halfV,
		const Vector4 &x1, float &contribution, const Vector2i *examplePos);

	inline Matrix2x2 computeJExample(float x, float y);

	inline int mod1(int x, int y) {
		return x - floor((float)x / y) * y;
	}

	Vector2 ByExampleProceduralNoise(const Vector2i &_uv, Matrix2x2 &Joutput);
	Vector2 ByExampleProceduralNoiseVariance(const Vector2i &_uv, Matrix2x2 &J);

	inline Vector2 getExampleNormal(Vector2i uv)
	{
		return m_example_normal->getPixel(uv.x, uv.y);
	}

	void constrainPRF(Float &sigma_p1, Float &sigma_p2);

	inline Vector2f mul(Vector2f vecA, Matrix2x2 matB)
	{
		Vector2f result;
		result.x = vecA[0] * matB.m[0][0] + vecA[1] * matB.m[1][0];
		result.y = vecA[0] * matB.m[0][1] + vecA[1] * matB.m[1][1];
		return result;
	}

	bool is_valid(const Vector2& projection) {
		return projection.x * projection.x + projection.y * projection.y < 1.0;
	}

	Vector v_unproject(const Vector2& projectedDirection)
	{
		return Vector(
			projectedDirection.x, projectedDirection.y,
			sqrtf(abs(1.0 - projectedDirection.x * projectedDirection.x - projectedDirection.y * projectedDirection.y)));
	}

	Vector2 hash1(Vector2i p);

	//compute the ID / weight of the target patch
	inline void QuadWeight(const Vector2i &uv, Vector4 &w);
	inline void QuadGrid(const Vector2i &uv, Vector4 &w, Vector2i &vertex1);

	float u_intrinsic;

	float delta = 1.0f;  //change from 0.5f to 1.0f to avoid interpolation
	string path;

	int m_res_target_width;
	int m_res_target_height;

	float m_inv_target_width;
	float m_inv_target_height;

	int m_res_example_width;
	int m_res_example_height;

	float m_sigmaH;
	int m_step;

	int m_res_target_p_width;
	int m_res_target_p_height;

	float m_invres_target_p_width;
	float m_invres_target_p_height;

	int m_bit_target_p_width;;
	int m_bit_target_p_height;

	Vector2i m_example_to_target_patch[4];

	std::vector<vector<Vector2i>> m_hash;

	int m_blend_type;

};

#endif