#pragma once
#ifndef RAY_TRACER
#define RAY_TRACER

#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"
#include "CubeMap.h"
#include "sphericalMap.h"

/*! \class Raytracer
\brief General ray tracer class.

\author Tomáš Fabián
\version 0.1
\date 2018
*/


struct IntersectionInfo { Vector3 intersectionPoint, vectorToLight, normal, viewToIntersectionVector; float dstToLight, enlighted; Material* material; Coord2f tex_coord; };

class Raytracer : public SimpleGuiDX11
{
public:
	Raytracer(const int width, const int height,
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=0,verbose=3");
	~Raytracer();

	void LoadScene(const std::string file_name);

	Color4f get_pixel(const int x, const int y, const float t = 0.0f) override;

	int Ui();

	Vector3 getInterpolatedPoint(RTCRay ray);
private:
	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_;

	int const antiAliasingSubSamplingConst = 80;
	float const antialiasingNormalizingCoef = 1 / float(antiAliasingSubSamplingConst);
	float const sRGBToLinearPower = 1 / 2.4f;

	float castShadowRay(IntersectionInfo intersectionInfo, RTCIntersectContext context);

	Color4f applyShader(const int x, const int y, const float t);
	Color4f applyShaderInternal(RTCRayHitWithIor rtcRayHitWithIor, float t, int depth, bool* isObjectSource, IntersectionInfo* nextIntersectionInfo = NULL);

	IntersectionInfo getIntersectionInfo(RTCRayHitWithIor rtcRayHitWithIor, Vector3 vectorFromCamera, Vector3 lightPossition, RTCIntersectContext context);

	Color4f applyPhongShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t, int depth);
	Color4f applyGlassShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t, int depth);
	Color4f applyWhittedShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t, int depth);
	Color4f applyNormalShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t);
	Color4f applyLambertShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t);
	Color4f applyPhysicallyBasedShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t, int depth);

	Color4f getAttenuationOfRay(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float actualIor, Material* material);

	Vector3 sampleHemisphere(Vector3 normal);

	int InitDeviceAndScene(const char * config);
	int ReleaseDeviceAndScene();

	CubeMap* cubeMap;
	SphericalMap* sphericalMap;
	RTCDevice device_;
	RTCScene scene_;
	Camera camera_;
};

#endif // !RAY_TRACER