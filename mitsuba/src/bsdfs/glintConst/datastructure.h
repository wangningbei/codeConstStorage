
#pragma once
#if !defined(__DATASTRUCTURE_H_)
#define __DATASTRUCTURE_H_

/*****************************************************************************/
/*****************************************************************************/
/********************************* Includes **********************************/
/*****************************************************************************/
/*****************************************************************************/

#include <vector>
#include <algorithm>
#include <cmath>

#include "mitsuba/mitsuba.h"


using namespace std;
using namespace mitsuba;


/*****************************************************************************/
/*****************************************************************************/
/******************************** Parameters *********************************/
/*****************************************************************************/
/*****************************************************************************/

#define USE_DXT_COMPRESSION false // Use DXT1 (true) or GL_RGB8 (false) (Section 1.6)


/*****************************************************************************/
/*****************************************************************************/
/*************************** Types and Structures ****************************/
/*****************************************************************************/
/*****************************************************************************/
struct vec2i;

struct vec2
{
	vec2() : x(0), y(0) {}
	vec2(float a, float b) : x(a), y(b) {}
	vec2(int a, int b) : x((float)a), y((float)b) {}

	float x;
	float y;

	vec2& operator*=(const float& value)
	{
		x *= value;
		y *= value;
		return *this;
	}

	vec2& operator+=(const vec2& other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}

	vec2& operator-=(const vec2& other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}

	vec2 operator-(const vec2& other) const
	{
		return vec2(x - other.x, y - other.y);
	}

	vec2 operator*(const float& value) const
	{
		return vec2(x * value, y * value);
	}

	vec2 operator+(const vec2& other) const
	{
		return vec2(x + other.x, y + other.y);
	}
	vec2& operator-=(const float& value)
	{
		x -= value;
		y -= value;
		return *this;
	}
	vec2& operator/=(const float& value)
	{
		x /= value;
		y /= value;
		
		return *this;
	}
	vec2& operator+=(const float& value)
	{
		x += value;
		y += value;
		return *this;
	}
};

struct vec2i
{
	vec2i() : x(0), y(0) {}
	vec2i(int a, int b) : x(a), y(b) {}
	int x;
	int y;

	vec2i& operator*=(const float& value)
	{
		x *= value;
		y *= value;
		return *this;
	}

	vec2i& operator+=(const vec2i& other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}

	vec2i& operator-=(const vec2i& other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}

	vec2i operator-(const vec2i& other) const
	{
		return vec2i(x - other.x, y - other.y);
	}

	vec2i operator*(const float& value) const
	{
		return vec2i(x * value, y * value);
	}

};


struct vec3
{
	vec3() : x(0), y(0), z(0) {}
	vec3(float a, float b, float c) : x(a), y(b), z(c) {}
	float x;
	float y;
	float z;

	float& operator[](const int i)
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		else
			return z;
	}

	float Dot(const vec3& other)
	{
		return x * other.x + y * other.y + z * other.z;
	}

	vec3& operator*=(const float& value)
	{
		x *= value;
		y *= value;
		z *= value;
		return *this;
	}

	vec3& operator/=(const float& value)
	{
		x /= value;
		y /= value;
		z /= value;
		return *this;
	}
	vec3& operator+=(const float& value)
	{
		x += value;
		y += value;
		z += value;
		return *this;
	}
	vec3& operator-=(const float& value)
	{
		x -= value;
		y -= value;
		z -= value;
		return *this;
	}
	vec3 operator+(const vec3& other) const
	{
		return vec3(x + other.x, y + other.y, z + other.z);
	}

	vec3 operator*(const float& value) const
	{
		return vec3(x * value, y * value, z * value);
	}
};

struct vec4
{
	vec4() : x(0), y(0), z(0), w(0) {}
	vec4(float a, float b, float c, float d) : x(a), y(b), z(c),w(d) {}
	float x;
	float y;
	float z;
	float w;

	float& operator[](const int i)
	{
		if (i == 0)
			return x;
		else if (i == 1)
			return y;
		else if(i == 2)
			return z;
		else return w;
	}

	float Dot(const vec4& other)
	{
		return x * other.x + y * other.y + z * other.z + w* other.w;
	}

	vec4& operator*=(const float& value)
	{
		x *= value;
		y *= value;
		z *= value;
		w *= value;
		return *this;
	}

	vec4& operator/=(const float& value)
	{
		x /= value;
		y /= value;
		z /= value;
		w /= value;
		return *this;
	}
	vec4& operator+=(const float& value)
	{
		x += value;
		y += value;
		z += value;
		w += value;
		return *this;
	}
	vec4& operator-=(const float& value)
	{
		x -= value;
		y -= value;
		z -= value;
		w -= value;
		return *this;
	}
	vec4 operator+(const vec4& other) const
	{
		return vec4(x + other.x, y + other.y, z + other.z, w + other.w);
	}

	vec4 operator*(const float& value) const
	{
		return vec4(x * value, y * value, z * value, w * value);
	}
};


class TextureDataFloat
{
public:
	TextureDataFloat() : data(), width(0), height(0), channels(0) {}
	TextureDataFloat(const int w, const int h, const int c) :
		data(w * h * c), width(w), height(h), channels(c)
	{
	}

	~TextureDataFloat(){
		vector<float>().swap(data);
	}

	float GetPixel(int w, int h, int c)
	{
		return data[h * width * channels + w * channels + c];
	}

	inline float GetPixel(int w, int h, int c)const
	{
		return data[h * width * channels + w * channels + c];
	}
	vec3 GetColorAt(int w, int h)
	{
		return vec3(
			data[h * width * channels + w * channels + 0],
			data[h * width * channels + w * channels + 1],
			data[h * width * channels + w * channels + 2]);
	}
	vec3 GetColorAt(int w, int h) const
	{
		return vec3(
			data[h * width * channels + w * channels + 0],
			data[h * width * channels + w * channels + 1],
			data[h * width * channels + w * channels + 2]);
	}

	inline Vector2 Get2ColorAt(int w, int h) const
	{
		int idx = h * width * channels + w * channels;
		return Vector2(data[idx], data[idx + 1]);
	}

	inline Vector4 Get4ColorAt(int w, int h) const
	{
		int idx = h * width * channels + w * channels;
		return Vector4(data[idx], data[idx + 1], data[idx + 2], data[idx + 3]);
	}

	void SetPixel(int w, int h, int c, float value)
	{
		data[h * width * channels + w * channels + c] = value;
	}
	void SetColorAt(int w, int h, vec3 value)
	{
		data[h * width * channels + w * channels + 0] = value.x;
		data[h * width * channels + w * channels + 1] = value.y;
		data[h * width * channels + w * channels + 2] = value.z;
	}
	void SetColorAt(int w, int h, vec4 value)
	{
		data[h * width * channels + w * channels + 0] = value.x;
		data[h * width * channels + w * channels + 1] = value.y;
		data[h * width * channels + w * channels + 2] = value.z;
		data[h * width * channels + w * channels + 3] = value.w;
	}
	vector<float> data;

	int width;
	int height;
	int channels;
};

struct PixelSortStruct
{
	int x;
	int y;
	float value;
	bool operator < (const PixelSortStruct& other) const
	{
		return (value < other.value);
	}
};

int ReadEXRtoTextureData(string path, TextureDataFloat& out, float scale);

#endif