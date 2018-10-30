#pragma once

#ifndef CUBE_MAP
#define CUBE_MAP

#include "texture.h"
#include "vector3.h"

using namespace std;

class CubeMap
{
public:
	Texture* textures[6];

	CubeMap(string posxFile, string negxFile, string posyFile, string negyFile, string poszFile, string negzFile);

	Color4f getTexel(Vector3 viewToIntersectionVector);

	~CubeMap();

private:
	enum CubePossition {
		POS_X = 0,
		NEG_X = 1,
		POS_Y = 2,
		NEG_Y = 3,
		POS_Z = 4,
		NEG_Z = 5
	};

	void setUVForLargestComponent(float largestComponent, float componentU, float componentV, float* u, float *v);
};
#endif // !CUBE_MAP