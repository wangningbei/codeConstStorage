#include <fstream>
#include "normalmap.h"
#include "flake.h"


NormalMap::NormalMap(string path)
{
	m_flakes = NULL;
	m_delta = 1.0; //this is for generating the J

	int width, height;
	try {
		Imf::RgbaInputFile file(path.c_str());
		Imath::Box2i dw = file.dataWindow();
		m_width = dw.max.x - dw.min.x + 1;
		m_height = dw.max.y - dw.min.y + 1;
		m_normal_map.resizeErase(m_height, m_width);
		printf("%d %d\n", m_width, m_height);

		file.setFrameBuffer(&m_normal_map[0][0] - dw.min.x - dw.min.y * m_width, 1, m_width);
		file.readPixels(dw.min.y, dw.max.y);
		printf("Done image reading\n");
	}
	catch (Iex::BaseExc &e) {
		std::cerr << e.what() << std::endl;
	}
}


void NormalMap::constructFlakes(float _intrinsic)
{
	if (m_flakes == NULL)
	{
		m_flakes = new Flakes(this, _intrinsic);
	}
}

Vector4f NormalMap::queryRangeNorMinMax(const Vector2i&start, const Vector2i& end, int h) //the index in the array)
{
	Vector4 Nmin_max;
	if (m_min_max_pretable_packed_3D[h][start.x][start.y] != -1)
	{
		unpack(m_min_max_pretable_packed_3D[h][start.x][start.y], Nmin_max);
		return Nmin_max;

	}
	Vector2 Nmin(1e7), Nmax(-1e7);
	
	if (h == 0)
	{
		//there is only one texel, read it.
		int index = start.x * m_height + start.y;
		Vector4f Nmin_max = m_flakes->getFlakeMinMax(index);
		return Nmin_max;
	}
	else
	{
		//subdivide the square with size power of 2
		int newH = h - 1;
		int newLength = (1 << newH);
		int midX = start.x + newLength;
		int midY = start.y + newLength;

		Vector4f min_max[4];
		min_max[0] = queryRangeNorMinMax(start, Vector2i(midX - 1, midY - 1), h - 1);
		min_max[1] = queryRangeNorMinMax(Vector2i(start.x, midY), Vector2i(midX - 1, end.y), h - 1);
		min_max[2] = queryRangeNorMinMax(Vector2i(midX, start.y), Vector2i(end.x, midY - 1), h - 1);
		min_max[3] = queryRangeNorMinMax(Vector2i(midX, midY), end, h - 1);

		for (int i = 0; i < 4; i++)
		{
			Nmin[0] = min(Nmin[0], min_max[i][0]);
			Nmin[1] = min(Nmin[1], min_max[i][1]);
			Nmax[0] = max(Nmax[0], min_max[i][2]);
			Nmax[1] = max(Nmax[1], min_max[i][3]);

		}

		Nmin_max = Vector4f(Nmin[0], Nmin[1], Nmax[0], Nmax[1]);
		if (h < m_min_max_pretable_packed_3D.size())
		{
			//some scaling for packing.
			Nmin = (Nmin + Vector2(3.0)) / 6.0f;
			Nmax = (Nmax + Vector2(3.0)) / 6.0f;

			uint32_t aScaled = Nmin.x * SCALEPACKING;
			uint16_t bScaled = Nmin.y * SCALEPACKING;
			uint32_t cScaled = Nmax.x * SCALEPACKING;
			uint16_t dScaled = Nmax.y * SCALEPACKING;

			uint32_t abPacked = (aScaled << 16) | (bScaled & 0xFFFF);
			uint32_t cdPacked = (cScaled << 16) | (dScaled & 0xFFFF);

			int64_t abcdPacked = ((int64_t)abPacked << 32) | (cdPacked);

			m_min_max_pretable_packed_3D[h][start.x][start.y] = abcdPacked;
		}
	}

	return Nmin_max;

}

void NormalMap::precomputeMinMaxTable(){

	//the precomputation table with 3D
	const int n = m_width;
	const int m = m_height;
	int patchSize = m_width / 4;
	const int k = log2(patchSize) + 1; 

	m_min_max_pretable_packed_3D.resize(k);

#pragma omp parallel for
	for (int i = 0; i < k; i++)
	{
		m_min_max_pretable_packed_3D[i].resize(n);
		for (int j = 0; j < n; j++)
		{
			m_min_max_pretable_packed_3D[i][j].resize(m, -1);
		}
	}
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			int h = k - 1;
			int length = (1 << h);

			while (i + length >= m_width || j + length >= m_width)
			{
				h--;
				length = (1 << h);
			}
			if (h < 0)
				continue;

			Vector4 Nmin_max = queryRangeNorMinMax(Vector2i(i, j), Vector2i(i + length - 1, j + length - 1), h);
		}
	}
}

void NormalMap::scaleNormalMap(float scale)
{
#pragma omp parallel for
	for (int u = 0; u < m_width; u++)
	{
		for (int v = 0; v < m_height; v++)
		{
			Vector3 normal = Vector3(m_normal_map[v][u].r, m_normal_map[v][u].g, m_normal_map[v][u].b);
			normal /= normal.z;

			normal.x *= scale;
			normal.y *= scale;

			normal = normalize(normal);

			m_normal_map[v][u].r = normal.x;
			m_normal_map[v][u].g = normal.y;
			m_normal_map[v][u].b = normal.z;
		}
	}
}

