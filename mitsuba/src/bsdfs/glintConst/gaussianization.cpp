
#include "gaussianization.h"

void Gaussianizaition::setInputData(const TextureDataFloat &_input)
{
	m_input = TextureDataFloat(_input.width, _input.height, _input.channels);
	for (int i = 0; i < m_input.data.size(); i++)
		m_input.data[i] = _input.data[i];

	m_Tinput = TextureDataFloat(m_input.width, m_input.height, _input.channels);
	m_Tinv = TextureDataFloat(LUT_WIDTH, 1, 3);

	m_T_prime = TextureDataFloat(m_input.width, m_input.height, 4);
	m_Tinv_prime = TextureDataFloat(LUT_WIDTH, 1, 2);

	// Sort pixels of example image
	for (int i = 0; i < 3; i++)
	{
		m_sortedValues[i].resize(_input.width * _input.height);
		m_sortedInputValues[i].resize(_input.width * _input.height);
		for (int y = 0; y < _input.height; y++)
		for (int x = 0; x < _input.width; x++)
		{
			m_sortedInputValues[i][y * _input.width + x].x = x;
			m_sortedInputValues[i][y * _input.width + x].y = y;
			m_sortedInputValues[i][y * _input.width + x].value = _input.GetPixel(x, y, i);
		}
		sort(m_sortedInputValues[i].begin(), m_sortedInputValues[i].end());

		for (int k = 0; k < _input.width * _input.height; k++)
		{
			m_sortedValues[i][k] = m_sortedInputValues[i][k].value;

		}
	}
}

/*********************************************************************/
/*********************************************************************/
/*************************** Main Function ***************************/
/*********************************************************************/
/*********************************************************************/

void Gaussianizaition::Precomputations()
{
	// Section 1.3.2 Applying the histogram transformation T on the input
	for (int channel = 0; channel < 3; channel++)
	{
		ComputeTinput(m_input, m_Tinput, channel);
	}

	// Section 1.3.3 Precomputing the inverse histogram transformation T^{-1}
	for (int channel = 0; channel < 3; channel++)
	{
		ComputeinvT(m_input, m_Tinv, channel);
	}

	computeTPrime();
	computeInvTPrime();
}

/*****************************************************************************/
/*****************************************************************************/
/**************** Section 1.3.1 Target Gaussian distribution *****************/
/*****************************************************************************/
/*****************************************************************************/

float Gaussianizaition::Erf(float x)
{
	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = abs(x);

	// A&S formula 7.1.26
	float t = 1.0f / (1.0f + 0.3275911f * x);
	float y = 1.0f - (((((1.061405429f * t + -1.453152027f) * t) + 1.421413741f)
		* t + -0.284496736f) * t + 0.254829592f) * t * exp(-x * x);

	return sign * y;
}

float Gaussianizaition::ErfInv(float x)
{
	float w, p;
	w = -log((1.0f - x) * (1.0f + x));
	if (w < 5.000000f)
	{
		w = w - 2.500000f;
		p = 2.81022636e-08f;
		p = 3.43273939e-07f + p * w;
		p = -3.5233877e-06f + p * w;
		p = -4.39150654e-06f + p * w;
		p = 0.00021858087f + p * w;
		p = -0.00125372503f + p * w;
		p = -0.00417768164f + p * w;
		p = 0.246640727f + p * w;
		p = 1.50140941f + p * w;
	}
	else
	{
		w = sqrt(w) - 3.000000f;
		p = -0.000200214257f;
		p = 0.000100950558f + p * w;
		p = 0.00134934322f + p * w;
		p = -0.00367342844f + p * w;
		p = 0.00573950773f + p * w;
		p = -0.0076224613f + p * w;
		p = 0.00943887047f + p * w;
		p = 1.00167406f + p * w;
		p = 2.83297682f + p * w;
	}
	return p * x;
}

float Gaussianizaition::CDF(float x, float mu, float sigma)
{
	float U = 0.5f * (1 + Erf((x - mu) / (sigma*sqrtf(2.0f))));
	return U;
}

float Gaussianizaition::invCDF(float U, float mu, float sigma)
{
	float x = sigma*sqrtf(2.0f) * ErfInv(2.0f*U - 1.0f) + mu;
	return x;
}

/*****************************************************************************/
/*****************************************************************************/
/**** Section 1.3.2 Applying the histogram transformation T on the input *****/
/*****************************************************************************/
/*****************************************************************************/

void Gaussianizaition::ComputeTinput(TextureDataFloat& input, TextureDataFloat& T_input, int channel)
{
	// Assign Gaussian value to each pixel
	for (unsigned int i = 0; i < m_sortedInputValues[channel].size(); i++)
	{
		// Pixel coordinates
		int x = m_sortedInputValues[channel][i].x;
		int y = m_sortedInputValues[channel][i].y;
		// Input quantile (given by its order in the sorting)
		float U = (i + 0.5f) / (m_sortedInputValues[channel].size());
		// Gaussian quantile
		float G = GAUSSIAN_SCALE * invCDF(U, GAUSSIAN_AVERAGE, GAUSSIAN_STD);
		// Store
		T_input.SetPixel(x, y, channel, G);
	}
}

void Gaussianizaition::computeTPrime()
{
	const int delta = 1;
	for (int x = 0; x < m_Tinput.width; x++)
	{
		int x_plus_1 = std::min(m_Tinput.width - 1, x + delta);
		int x_miunus_1 = std::max(0, x - delta);

		for (int y = 0; y < m_Tinput.height; y++)
		{
			Vector2f NX_plus_1_Y = m_Tinput.Get2ColorAt(x_plus_1, y);
			Vector2f NX_minus_1_Y = m_Tinput.Get2ColorAt(x_miunus_1, y);

			int y_plus_1 = std::min(m_Tinput.height - 1, y + delta);
			int y_miunus_1 = std::max(0, y - delta);

			Vector2f NX_Y_PLUS_1 = m_Tinput.Get2ColorAt(x, y_plus_1);
			Vector2f NX_Y_MINUS_1 = m_Tinput.Get2ColorAt(x, y_miunus_1);

			Vector2f J1 = (NX_plus_1_Y - NX_minus_1_Y) / (x_plus_1 - x_miunus_1);
			Vector2f J2 = (NX_Y_PLUS_1 - NX_Y_MINUS_1) / (y_plus_1 - y_miunus_1);

	//		Matrix2x2 J(J1.x, J2.x, J1.y, J2.y);
			m_T_prime.SetColorAt(x, y, vec4(J1.x, J2.x, J1.y, J2.y));
		}
	}
}

void Gaussianizaition::computeInvTPrime()
{
	const int delta = 1;
	int y = 0;
	for (int x = 0; x <m_Tinv.width; x++)
	{
		int x_plus_1 = std::min(m_Tinv.width - 1, x + delta);
		int x_miunus_1 = std::max(0, x - delta);

		float NX_plus_1_Y_x = m_Tinv.data[x_plus_1 * 3];
		float NX_minus_1_Y_x = m_Tinv.data[x_miunus_1 * 3];

		float NX_plus_1_Y_y = m_Tinv.data[x_plus_1 * 3 + 1];
		float NX_minus_1_Y_y = m_Tinv.data[x_miunus_1 * 3 + 1];

		float J1 = (NX_plus_1_Y_x - NX_minus_1_Y_x) * m_Tinv.width / (x_plus_1 - x_miunus_1);
		float J2 = (NX_plus_1_Y_y - NX_minus_1_Y_y) * m_Tinv.width / (x_plus_1 - x_miunus_1);

		m_Tinv_prime.SetPixel(x, 0, 0, J1);
		m_Tinv_prime.SetPixel(x, 0, 1, J2);

	}
}
/*****************************************************************************/
/*****************************************************************************/
/*  Section 1.3.3 Precomputing the inverse histogram transformation T^{-1}   */
/*****************************************************************************/
/*****************************************************************************/

void Gaussianizaition::ComputeinvT(TextureDataFloat& input, TextureDataFloat& Tinv, int channel)
{
	// Sort pixels of example image
	vector<float> sortedInputValues;
	sortedInputValues.resize(input.width * input.height);
	for (int y = 0; y < input.height; y++)
		for (int x = 0; x < input.width; x++)
		{
			sortedInputValues[y * input.width + x] = input.GetPixel(x, y, channel);
		}
	sort(sortedInputValues.begin(), sortedInputValues.end());

	// Generate Tinv look-up table 
	for (int i = 0; i < Tinv.width; i++)
	{
		// Gaussian value in [0, 1]
		float G = GAUSSIAN_SCALE * (i + 0.5f) / (Tinv.width);
		// Quantile value 
		float U = CDF(G, GAUSSIAN_AVERAGE, GAUSSIAN_STD);
		// Find quantile in sorted pixel values
		int index = (int)floor(U * sortedInputValues.size());
		// Get input value 
		float I = sortedInputValues[index];
		// Store in LUT
		Tinv.SetPixel(i, 0, channel, I);
	}
}

Vector2 Gaussianizaition::getG(Vector2 uv) const
{
	Vector2 uv1 = uv;

	uv1.x -= floor(uv.x);
	uv1.y -= floor(uv.y);

	vec2i uv10 = vec2i(floor(uv1.x * m_Tinput.width), floor(uv1.y * m_Tinput.height));

	return m_Tinput.Get2ColorAt(uv10.x, uv10.y);
}


Vector2 Gaussianizaition::getGFromValue(Vector2 valueN) const
{
	float result[2];
	
	for (int i = 0; i < 2; i++)
	{
		//from the value to the index for each channel seperately
		const vector<PixelSortStruct> &channelFullData = m_sortedInputValues[i];
		const vector<float> &channelValueData = m_sortedValues[i];

		int index = std::min((int)channelFullData.size() - 1, (int)(std::lower_bound(channelValueData.begin(),
			channelValueData.end(), valueN[i]) - channelValueData.begin()));
	
		const PixelSortStruct &texelData = channelFullData[index];

		if (texelData.value == valueN[i] || index == 0 || index == channelFullData.size() - 1)
		{			
			int u = texelData.x;
			int v = texelData.y;
			result[i] = m_Tinput.GetPixel(u, v, i);
			continue;
		}

		float preValue = channelFullData[index - 1].value;
		float nextValue = texelData.value;

		float G0 = m_Tinput.GetPixel(channelFullData[index - 1].x, channelFullData[index - 1].y, i);
		float G1 = m_Tinput.GetPixel(texelData.x, texelData.y, i);

		float weight = (preValue != nextValue) ? (nextValue - valueN[i]) / (nextValue - preValue) : 1.0f;
		result[i] = G0 * weight + G1 * (1.0f - weight);
	}

	return Vector2(result[0], result[1]);
}


bool Gaussianizaition::boundByGaussian(Vector2 input)const
{
	for (int i = 0; i < 2; i++)
	{
		//from the value to the index for each channel seperately
		const vector<PixelSortStruct> &channelFullData = m_sortedInputValues[i];

		if (input[i] < channelFullData[0].value || input[i] > channelFullData[channelFullData.size() - 1].value)
		{
			return false;
		}
	}

	return true;
	
}
