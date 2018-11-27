
#ifndef SPHERICAL_MAP
#define SPHERICAL_MAP

#include "stdafx.h"
#include "texture.h"

#define _USE_MATH_DEFINES
#include <math.h>

class SphericalMap
{
public:
	SphericalMap() { 
		texture = new Texture("../../../data/green_point_park.jpg");
	}

	Color4f getTexel(Vector3 viewToIntersectionVector) {
		float theta = acosf(viewToIntersectionVector.z) / float(M_PI);
		float phi = (atan2f(viewToIntersectionVector.y, viewToIntersectionVector.x) + float(M_PI)) / (2 * float(M_PI));

		Color3f color =  texture->get_texel(1 - phi, theta);
		return Color4f {color.r, color.g, color.b, 1};
	}


private:
	Texture* texture;
};
#endif
