#include "datastructure.h"
#include <iostream>
#include <sstream>
#include <fstream>

#include <OpenEXRConfig.h>
#include <ImfRgbaFile.h>
#include <ImfArray.h>
#include <ImfNamespace.h>
#include <half.h>

int ReadEXRtoTextureData(string path, TextureDataFloat& out, float scale) {
	Imf::Array2D<Imf::Rgba> exrImage;
	int width, height;
	try {
		Imf::RgbaInputFile file(path.c_str());
		Imath::Box2i dw = file.dataWindow();
		width = dw.max.x - dw.min.x + 1;
		height = dw.max.y - dw.min.y + 1;

		printf("%d %d\n", width, height);
		exrImage.resizeErase(height, width);
		file.setFrameBuffer(&exrImage[0][0] - dw.min.x - dw.min.y * width, 1, width);
		file.readPixels(dw.min.y, dw.max.y);
		printf("Done image reading\n");
	}
	catch (Iex::BaseExc &e) {
		std::cerr << e.what() << std::endl;
	}

	out.channels = 3;
	out.width = width;
	out.height = height;
	// Store in floating point structure
	out.data.resize((out.width) * (out.height) * 3);

	for (int v = 0; v < height; v += 1) {
		for (int u = 0; u < width; u += 1) {
			int index = (v * width + u) * 3;

			//scale normal map
			Vector3 normal = Vector3(exrImage[v][u].r, exrImage[v][u].g, exrImage[v][u].b);
			normal /= normal.z;

			normal.x *= scale;
			normal.y *= scale;

			normal = normalize(normal);

			out.data[index] = normal.x;
			out.data[index + 1] = normal.y;
			out.data[index + 2] = normal.z;
		}
	}

	printf("Done texture reading\n");
	return 0;
}

