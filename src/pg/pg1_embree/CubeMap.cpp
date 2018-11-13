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


Color4f CubeMap::getTexel(Vector3 toIntersectionVector) {
	//viewToIntersectionVector.z = -viewToIntersectionVector.z;
	float u, v;
	int mapPossition;
	char axis = toIntersectionVector.LargestComponent(true);
	int orientationOffset = toIntersectionVector.data[axis] < 0 ? 1 : 0;



	mapPossition = 2 * axis + orientationOffset;
	CubePossition cubePossition = CubePossition(mapPossition);

	switch (cubePossition)
	{
	case CubeMap::POS_X:
		setUVForLargestComponent(toIntersectionVector.x, -toIntersectionVector.z, toIntersectionVector.y, &u, &v);
		break;
	case CubeMap::NEG_X:
		setUVForLargestComponent(toIntersectionVector.x, toIntersectionVector.z, toIntersectionVector.y, &u, &v);
		break;
	case CubeMap::POS_Y:
		setUVForLargestComponent(toIntersectionVector.y, toIntersectionVector.x, -toIntersectionVector.z, &u, &v);
		break;
	case CubeMap::NEG_Y:
		setUVForLargestComponent(toIntersectionVector.y, toIntersectionVector.x, toIntersectionVector.z, &u, &v);
		break;
	case CubeMap::POS_Z:
		setUVForLargestComponent(toIntersectionVector.z, toIntersectionVector.x, toIntersectionVector.y, &u, &v);
		break;
	case CubeMap::NEG_Z:
		setUVForLargestComponent(toIntersectionVector.z, -toIntersectionVector.x, toIntersectionVector.y, &u, &v);
		break;
	}
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
	for (int i = 0; i < 6; i++) {
		delete this->textures[i];
	}
}