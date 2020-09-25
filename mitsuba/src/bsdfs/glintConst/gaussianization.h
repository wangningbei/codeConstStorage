/*****************************************************************************/
//This is from Eric Heitz paper: Histogram preserving blending.
/*****************************************************************************/
#pragma once
#if !defined(__GAUSSIANIZATION_H_)
#define __GAUSSIANIZATION_H_

/*****************************************************************************/
/*****************************************************************************/
/********************************* Includes **********************************/
/*****************************************************************************/
/*****************************************************************************/

#include <vector>
#include <algorithm>
#include <cmath>

#include "datastructure.h"

/*****************************************************************************/
/*****************************************************************************/
/******************************** Parameters *********************************/
/*****************************************************************************/
/*****************************************************************************/

#define USE_DXT_COMPRESSION false // Use DXT1 (true) or GL_RGB8 (false) (Section 1.6)
#define GAUSSIAN_AVERAGE 0.5f // Expectation of the Gaussian distribution
#define GAUSSIAN_SCALE 1.0f
#define GAUSSIAN_STD 0.16666f 
//0.16666f // Std of the Gaussian distribution
#define LUT_WIDTH 128 // Size of the look-up table

using namespace std;

class Gaussianizaition {

public:
	Gaussianizaition(){};
	~Gaussianizaition(){

		for (int i = 0; i < 3; i++)
		{
			vector<PixelSortStruct>().swap(m_sortedInputValues[3]);
			vector<float>().swap(m_sortedValues[3]);
		}
	
	}; 
	void setInputData(const TextureDataFloat &_input);

	void Precomputations();
	Vector2 getG(Vector2 uv) const;

	inline Vector2 getGInt(const Vector2i &uv) const
	{
		return m_Tinput.Get2ColorAt(uv.x, uv.y);
	}

	inline Vector4 getGIntPrime(const Vector2i &uv) const
	{
		return m_T_prime.Get4ColorAt(uv.x, uv.y);
	}

	inline Vector2 getJInt(const Vector2i &uv) const
	{
		return m_Tinput.Get2ColorAt(uv.x, uv.y);
	}

	bool boundByGaussian(Vector2 valueN)const;

	Vector2 getGFromValue(Vector2 valueN) const;

	inline Vector2 FromG2Org(const Vector2 &_G) const
	{
		int GX = _G.x * m_Tinv.width;
		int GY = _G.y * m_Tinv.width;
		GX = max(0, min(m_Tinv.width - 1, GX));
		GY = max(0, min(m_Tinv.width - 1, GY));

		GX *= 3; GY = GY * 3 + 1;

		Vector2 color;
		color.x = m_Tinv.data[GX];
		color.y = m_Tinv.data[GY];

		return color;
	}

	inline Vector2 FromG2OrgPrime(const Vector2 &_G) const
	{
		int GX = _G.x * m_Tinv.width;
		int GY = _G.y * m_Tinv.width;
		GX = max(0, min(m_Tinv.width - 1, GX));
		GY = max(0, min(m_Tinv.width - 1, GY));

		GX *= 2; GY = GY * 2 + 1;

		Vector2 color;
		color.x = m_Tinv_prime.data[GX];
		color.y = m_Tinv_prime.data[GY];

		return color;
	}

	const TextureDataFloat& getT() const{ return m_Tinput; };
	const TextureDataFloat& getInvT() const{ return m_Tinv; };

	TextureDataFloat& getT(){ return m_Tinput; };
	TextureDataFloat& getInvT(){ return m_Tinv; };

	size_t getSize(){
		return (m_Tinput.data.size() + m_Tinv.data.size()) * 4;//B	
	}
private:

	/*****************************************************************************/
	/*****************************************************************************/
	/**************** Section 1.3.1 Target Gaussian distribution *****************/
	/*****************************************************************************/
	/*****************************************************************************/

	float Erf(float x);

	float ErfInv(float x);

	float CDF(float x, float mu, float sigma);

	float invCDF(float U, float mu, float sigma);

	/*****************************************************************************/
	/*****************************************************************************/
	/**** Section 1.3.2 Applying the histogram transformation T on the input *****/
	/*****************************************************************************/
	/*****************************************************************************/

	void ComputeTinput(TextureDataFloat& input, TextureDataFloat& T_input, int channel);

	/*****************************************************************************/
	/*****************************************************************************/
	/*  Section 1.3.3 Precomputing the inverse histogram transformation T^{-1}   */
	/*****************************************************************************/
	/*****************************************************************************/

	void ComputeinvT(TextureDataFloat& input, TextureDataFloat& Tinv, int channel);
	void computeTPrime();
	void computeInvTPrime();

	TextureDataFloat m_input;
	TextureDataFloat m_Tinput;
	TextureDataFloat m_Tinv;

	TextureDataFloat m_T_prime;
	TextureDataFloat m_Tinv_prime;

	vector<PixelSortStruct> m_sortedInputValues[3];
	vector<float> m_sortedValues[3];
};


#endif