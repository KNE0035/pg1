#ifndef MY_MATH_H_
#define MY_MATH_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include "structs.h"
#include "vector3.h"

template <class T> inline T sqr( const T x )
{
	return x * x;
}

inline Normal3f normalize( const Normal3f & n )
{
	float tmp = sqr( n.x ) + sqr( n.y ) + sqr( n.z );

	if ( fabsf( tmp ) > FLT_EPSILON )
	{
		tmp = 1.0f / tmp;
		return Normal3f{ n.x * tmp, n.y * tmp, n.z * tmp };
	}

	return n;
}

inline float deg2rad( const float x )
{
	return x * float( M_PI ) / 180.0f;
}

inline Vector3 orthogonal(const Vector3 & v)
{
	return (abs(v.x) > abs(v.z)) ? Vector3(-v.y, v.x, 0.0f) : Vector3(0.0f, -v.z, v.y);
}

#endif
