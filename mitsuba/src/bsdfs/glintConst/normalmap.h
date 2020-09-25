#pragma once
#if !defined(__NORMALMAP_H_)
#define __NORMALMAP_H_

#include <ImfRgbaFile.h>
#include <ImfArray.h>
#include <ImfNamespace.h>
#include <OpenEXRConfig.h>

#include "mitsuba/core/platform.h"
#include "mitsuba/core/fwd.h"
#include "mitsuba/core/aabb.h"
#include "mitsuba/core/matrix.h"


//the sample rate MUST BE 1 in our case. 
#define SAMPLE_RATE  1
class Flakes;


using namespace std;
using namespace mitsuba;

const float SCALEPACKING = 65535.0f;

const int tab32[32] = {
	0, 9, 1, 10, 13, 21, 2, 29,
	11, 14, 16, 18, 22, 25, 3, 30,
	8, 12, 20, 28, 15, 17, 24, 7,
	19, 27, 23, 6, 26, 5, 4, 31 };

class NormalMap{

public:
	NormalMap(){
		m_flakes = NULL;
	}

	NormalMap(int _width, int _height){
		m_flakes = NULL;
		m_width = _width;
		m_height = _height;
		m_delta = 1.0; //this is for generating the J
		m_normal_map.resizeErase(_height, _width);
	}

	NormalMap(string path);

	~NormalMap(){ 

		if (m_flakes)
			delete m_flakes;

		for (int i = 0; i < m_min_max_pretable_packed_3D.size(); i++)
		{
			for (int j = 0; j < m_min_max_pretable_packed_3D[0].size(); j++)
			{
				vector<int64_t>().swap(m_min_max_pretable_packed_3D[i][j]);
			}
			vector<vector<int64_t>>().swap(m_min_max_pretable_packed_3D[i]);
		}
		vector<vector<vector<int64_t>>>().swap(m_min_max_pretable_packed_3D);

	}

	size_t getSize()
	{
		const int n = m_width;
		const int m = m_height;
		int patchSize = m_width / 4;
		const int k = log2(patchSize); 
		return n * m * k * 8;
	}
	inline int mod(int x, int y) {
		return ((x % y) + y) % y;
	}

	inline float p(int x, int y) {
		return m_normal_map[mod(y, m_height)][mod(x, m_width)].r;
	}
	inline float px(int x, int y) {
		return (m_normal_map[mod(y, m_height)][mod(x + 1, m_width)].r - m_normal_map[mod(y, m_height)][mod(x - 1, m_width)].r) / 2.0f;
	}

	inline float py(int x, int y) {
		return (m_normal_map[mod(y + 1, m_height)][mod(x, m_width)].r - m_normal_map[mod(y - 1, m_height)][mod(x, m_width)].r) / 2.0f;
	}

	inline float pxy(int x, int y) {
		return (p(x + 1, y + 1) - p(x + 1, y) - p(x, y + 1) + 2.0f * p(x, y) - p(x - 1, y) - p(x, y - 1) + p(x - 1, y - 1)) / 2.0f;
	}

	inline float q(int x, int y) {
		return m_normal_map[mod(y, m_height)][mod(x, m_width)].g;
	}

	inline float qx(int x, int y) {
		return (m_normal_map[mod(y, m_height)][mod(x + 1, m_width)].g - m_normal_map[mod(y, m_height)][mod(x - 1, m_width)].g) / 2.0f;
	}

	inline float qy(int x, int y) {
		return (m_normal_map[mod(y + 1, m_height)][mod(x, m_width)].g - m_normal_map[mod(y - 1, m_height)][mod(x, m_width)].g) / 2.0f;
	}

	inline float qxy(int x, int y) {
		return (q(x + 1, y + 1) - q(x + 1, y) - q(x, y + 1) + 2.0f * q(x, y) - q(x - 1, y) - q(x, y - 1) + q(x - 1, y - 1)) / 2.0f;
	}

	inline float getNx(float u, float v) {
		return p(mod(u, m_width), mod(v, m_height));
	}

	inline float getNy(float u, float v) {
		return q(mod(u, m_width), mod(v, m_height));
	}

	inline float getNxDx(float x, float y) {
		return (getNx(x + m_delta, y) - getNx(x - m_delta, y)) / (2.0f * m_delta);
	}

	inline float getNxDy(float x, float y) {
		return (getNx(x, y + m_delta) - getNx(x, y - m_delta)) / (2.0f * m_delta);
	}

	inline float getNyDx(float x, float y) {
		return (getNy(x + m_delta, y) - getNy(x - m_delta, y)) / (2.0f * m_delta);
	}

	inline float getNyDy(float x, float y) {
		return (getNy(x, y + m_delta) - getNy(x, y - m_delta)) / (2.0f * m_delta);
	}
	void scaleNormalMap(float scale);

	void setWidth(int _width){ m_width = _width; }
	void setHeight(int _height){ m_height = _height; }

	int getWidth(){ return m_width; } 
	int getHeight(){ return m_height; }
	int getWidth() const{ return m_width; }
	int getHeight() const{ return m_height; }

	void constructFlakes(float _intrinsic);

	//precompute the RMQ Table
	void precomputeMinMaxTable();
	Vector4f queryRangeNorMinMax(const Vector2i&start, const Vector2i&end, int h);

	// Get Jacobian matrix
	inline Matrix2x2 getJ(float x, float y)  {
		Matrix2x2 J(getNxDx(x, y), getNxDy(x, y), getNyDx(x, y), getNyDy(x, y));
		return J;
	}

	Vector2f getNormal(float u, float v)  {
		return Vector2f(p(mod(u* m_width, m_width), mod(v* m_height, m_height)), q(mod(u* m_width, m_width), mod(v* m_height, m_height)));
	}

	const vector<vector<int64_t>> &getMinMaxTable(int h) const
	{
		return m_min_max_pretable_packed_3D[h];
	}

	inline Vector2 getPixel(int u, int v) const{
		Imf::Rgba pixel = m_normal_map[v][u];
		return Vector2(pixel.r, pixel.g);
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


	Imf::Array2D<Imf::Rgba> m_normal_map;
	Flakes* m_flakes;

private:

	int m_width;
	int m_height;
	float m_delta;

	vector<vector<vector<int64_t>>>  m_min_max_pretable_packed_3D;
	int m_res1, m_res2, m_res3, m_res4;
};


#endif