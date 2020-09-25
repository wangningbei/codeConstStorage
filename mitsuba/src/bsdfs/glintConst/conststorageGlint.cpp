#include "blending.h"
#include <time.h>

/* =============================== MAIN ==================================== */
int main(int argc, char** argv)
{
	// Ask user for input texture file name
	string inputSampleFile = "";
	string outputFile = "";
	clock_t t;

	cout << argc;
	if (argc != 5) {
		cout << "example_filename_exr target_filename_exr target_res_width target_res_height" << endl;
		return 0;
	}

	inputSampleFile = argv[1];
	outputFile = argv[2];
	const int res_target_width = std::stoi(argv[3]);
	const int res_target_height = std::stoi(argv[4]);

	//set the blend type: EGaussian = 0, EVariance, ELinear,
	BlendType blendtype = EGaussian; 

	//scale the normals: larger value means more rough
	float scale = 0.5;

	BlendCPU *blender = new BlendCPU(inputSampleFile, blendtype, scale);

	cout << "Generate the entire target image for validation. \n";
	blender->setTargetInform(res_target_width, res_target_height);

	//this is for debug, we do not have to generate the entire target
	blender->generateEntireTargetImage(outputFile);

	const float intrinsic = 0.001;
	blender->setIntrinsic(intrinsic);

	printf("Start a NDF Visulization.\n");

	//This is the pixel's footprint: center and size, taking the resolution into consideration.
	float visTexels = 10.0f;
	Vector2i center(256, 256);

	printf("Construct the flakes.\n");
	blender->buildFlakes();

	printf("Start to precompute MinMax Table.\n");
	//use gaussian space or not
	if (blendtype == 0)
		blender->m_example_normal_Gaussian->precomputeMinMaxTable();
	else
		blender->m_example_normal->precomputeMinMaxTable();

	t = clock();

	printf("Start the PNDF visualization.\n");

	blender->visPNDF(blender->getTargetTexture(), center, visTexels);

	t = clock() - t;
	printf("It cost %f seconds).\n",((float)t) / CLOCKS_PER_SEC);
  
	delete blender;
	return 0;
}






