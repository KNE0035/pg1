#pragma once
#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"
#include "sphericalMap.h"
#include "utils.h"

/*! \class Raytracer
\brief General ray tracer class.

\author Tomáš Fabián
\version 0.1
\date 2018
*/
class Raytracer : public SimpleGuiDX11
{
public:
	Raytracer( const int width, const int height, 
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=0,verbose=3" );
	~Raytracer();

	void LoadScene( const std::string file_name );

	Color4f get_pixel( const int x, const int y, const float t = 0.0f ) override;

	int Ui();

	Vector3 getInterpolatedPoint(RTCRay ray);
private:
	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_; 
	
	float castShadowRay(const Vector3 origin, Vector3 vectorToLight, const float dist, RTCIntersectContext context);
	
	Color4f applyShader(const int x, const int y, const float t);
	Color4f applyShaderInternal(RTCRayHitWithIor rtcRayHitWithIor, float t, int depth);

	void getIntersectionInfo(RTCRayHitWithIor rtcRayHitWithIor, Vector3* vectorToLight, Vector3* normal, Vector3 viewVector, Vector3* intersectionPoint, Vector3 lightPossition, float* dstToLight, Material** material);


	int InitDeviceAndScene(const char * config);
	int ReleaseDeviceAndScene();

	SphericalMap* sphericalMap;
	RTCDevice device_;
	RTCScene scene_;
	Camera camera_;
};