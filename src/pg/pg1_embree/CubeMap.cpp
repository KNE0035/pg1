#include "stdafx.h"
#include "CubeMap.h"


CubeMap::CubeMap(string posxFile, string negxFile, string posyFile, string negyFile, string poszFile, string negzFile)
{
	this->textures[0] = new Texture(posxFile.c_str());
	this->textures[1] = new Texture(negxFile.c_str());
	
	this->textures[2] = new Texture(posyFile.c_str());
	this->textures[3] = new Texture(negyFile.c_str());

	this->textures[4] = new Texture(poszFile.c_str());
	this->textures[5] = new Texture(negzFile.c_str());
}


Color4f CubeMap::getTexel(Vector3 viewToIntersectionVector) {
	float u, v;
	int mapPossition;
	char axis = viewToIntersectionVector.LargestComponent(true);
	int orientationOffset = viewToIntersectionVector.data[axis] < 0 ? 1 : 0;



	mapPossition = 2 * axis + orientationOffset;
	CubePossition cubePossition = CubePossition(mapPossition);
	
	if (cubePossition == POS_X) {
		printf("test");
	}

	switch (cubePossition)
	{
	case CubeMap::POS_X:
	case CubeMap::NEG_X:
		setUVForLargestComponent(viewToIntersectionVector.x, viewToIntersectionVector.y, viewToIntersectionVector.z, &u, &v);
		break;
	case CubeMap::POS_Y:
	case CubeMap::NEG_Y:
		setUVForLargestComponent(viewToIntersectionVector.y, viewToIntersectionVector.z, viewToIntersectionVector.x, &u, &v);
		break;
	case CubeMap::POS_Z:
	case CubeMap::NEG_Z:
		setUVForLargestComponent(viewToIntersectionVector.z, viewToIntersectionVector.y, viewToIntersectionVector.x, &u, &v);
		break;
	}

	/*if (cubePossition == NEG_X || cubePossition == NEG_Y || cubePossition == NEG_Z) {
		u = 1 - u;
		v = 1 - v;
	}*/

	//u = 1 - u;
	//v = 1 - v;

	Color3f texel = this->textures[mapPossition]->get_texel(u, v);
	return Color4f {texel.r, texel.g, texel.b, 1};
}

void CubeMap::setUVForLargestComponent(float largestComponent, float componentU, float componentV, float* u, float *v) {
	float axixNormalizingFactor = 1.0f / abs(largestComponent);
	*u = (componentU * axixNormalizingFactor + 1) * 0.5f;
	*v = (componentV * axixNormalizingFactor + 1) * 0.5f;
}

CubeMap::~CubeMap()
{
	delete[] this->textures;
}


/*	switch (cubePossition)
	{
	case CubeMap::POS_X:
		setUVForLargestComponent(viewToIntersectionVector.x, viewToIntersectionVector.y, viewToIntersectionVector.z, &u, &v);
		break;
	case CubeMap::NEG_X:
		setUVForLargestComponent(viewToIntersectionVector.x, viewToIntersectionVector.y, viewToIntersectionVector.z, &u, &v);
		break;
	case CubeMap::POS_Y:
		setUVForLargestComponent(viewToIntersectionVector.y, viewToIntersectionVector.x, viewToIntersectionVector.z, &u, &v);
		break;
	case CubeMap::NEG_Y:
		setUVForLargestComponent(viewToIntersectionVector.y, viewToIntersectionVector.x, viewToIntersectionVector.z, &u, &v);
		break;
	case CubeMap::POS_Z:
		setUVForLargestComponent(viewToIntersectionVector.z, viewToIntersectionVector.y, viewToIntersectionVector.x, &u, &v);
		break;
	case CubeMap::NEG_Z:
		setUVForLargestComponent(viewToIntersectionVector.z, viewToIntersectionVector.y, viewToIntersectionVector.x, &u, &v);
		break;
	}*/