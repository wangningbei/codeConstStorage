# Constant Storage Glints

Description

This code implements the key ideas of the ACM Transactions on Graphics 2020 paper:

Example-Based Microstructure Rendering with Constant Storage", by

Beibei Wang, Miloš Hašan, Nicolas Holzschuch, and Ling-Qi Yan.

It's under the framework of Mitsuba and a new BSDF plug-in named glintyconductor_constant is added. Given a normal map as the input example, the blending type and targeting size, it implicitly generate the targeting normal map and hierarchy during rendering. 

Usage:

The code is under the framework of Mitsuba Renderer 0.5. It's only tested on Windows, so not sure about Mac OS. To compile this code, please use the following steps:

1. Find a dependencies for Mitsuba in the website of Mitsuba renderer, and then make it under the folder mitsuba.

2. Use CMake 2.8 to generate a solution of VS 2013. (We found CMake 3.0+ does not work in our case). 

3. Open the solution, and then build it. 


Example:

We also provided a testing scene TEAPOT. With the provided script, the scene could be rendered directly.

The XML can be found [here](scenes/teapot/teapot.xml), an input normal map can be found [here](scenes/teapot/normalmap/brush512.exr) and a binary (Windows) can be found [here](scenes/Release). 

The command is also shown here:

Release\mitsuba.exe teapot\teapot.xml -o teapot.exr


The settings of glintyconductor_constant in the XML is as follows:

	<bsdf type="glintyconductor_constant">
		
		<string name="distribution" value="beckmann"/>
		<float name="alpha" value="0.1"/>
		<string name="material" value="Cu"/>

		<!--intrinsic roughness-->
		<float name="intrinsic" value="0.005"/>

		<!--this is the size of the target normal map. For example, if the size of example normal map is 1K x 1K, then the size of target normal map is 32K x 32 K-->
		<float name="textureScale" value="32"/>

		<!--Blend type: 0 for histogram preserving, 1 for variance preserving and 2 for linear preserving-->
		<integer name="blendType" value="0"/>  
		
		<!--input example normal map-->
		<string name="nameNormalMap" value="normalmap\brush512.exr"/>
	</bsdf> 

Some other tips:

The flakes per texel should be set as one, as this is required by the Jacobian computation. 

The input normal map should have a range [-1,1]. One example of the normal map is included in the KETTLE folder.


