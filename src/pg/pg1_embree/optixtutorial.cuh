#include <optix_world.h>
//#include <optix.h>
//#include <optix_math.h>

struct PerRayData_radiance
{
	float3 result;
	float  importance;
	int depth;
};

struct PerRayData_shadow
{
	float3 attenuation;
};