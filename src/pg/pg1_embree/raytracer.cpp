#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include "mymath.h" 

Raytracer::Raytracer( const int width, const int height,
	const float fov_y, const Vector3 view_from, const Vector3 view_at,
	const char * config ) : SimpleGuiDX11( width, height )
{
	InitDeviceAndScene( config );

	camera_ = Camera( width, height, fov_y, view_from, view_at );
}

Raytracer::~Raytracer()
{
	ReleaseDeviceAndScene();
}

int Raytracer::InitDeviceAndScene( const char * config )
{
	device_ = rtcNewDevice( config );
	error_handler( nullptr, rtcGetDeviceError( device_ ), "Unable to create a new device.\n" );
	rtcSetDeviceErrorFunction( device_, error_handler, nullptr );

	ssize_t triangle_supported = rtcGetDeviceProperty( device_, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED );

	// create a new scene bound to the specified device
	scene_ = rtcNewScene( device_ );

	return S_OK;
}

int Raytracer::ReleaseDeviceAndScene()
{
	rtcReleaseScene( scene_ );
	rtcReleaseDevice( device_ );

	return S_OK;
}

void Raytracer::LoadScene( const std::string file_name )
{
	const int no_surfaces = LoadOBJ( file_name.c_str(), surfaces_, materials_ );
	sphericalMap = new SphericalMap();
	
	// surfaces loop
	for ( auto surface : surfaces_ )
	{
		RTCGeometry mesh = rtcNewGeometry( device_, RTC_GEOMETRY_TYPE_TRIANGLE );

		Vertex3f * vertices = ( Vertex3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
			sizeof( Vertex3f ), 3 * surface->no_triangles() );

		Triangle3ui * triangles = ( Triangle3ui * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			sizeof( Triangle3ui ), surface->no_triangles() );

		rtcSetGeometryUserData( mesh, ( void* )( surface->get_material() ) );

		rtcSetGeometryVertexAttributeCount( mesh, 2 );

		Normal3f * normals = ( Normal3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3,
			sizeof( Normal3f ), 3 * surface->no_triangles() );

		Coord2f * tex_coords = ( Coord2f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2,
			sizeof( Coord2f ), 3 * surface->no_triangles() );		

		// triangles loop
		for ( int i = 0, k = 0; i < surface->no_triangles(); ++i )
		{
			Triangle & triangle = surface->get_triangle( i );

			// vertices loop
			for ( int j = 0; j < 3; ++j, ++k )
			{
				const Vertex & vertex = triangle.vertex( j );

				vertices[k].x = vertex.position.x;
				vertices[k].y = vertex.position.y;
				vertices[k].z = vertex.position.z;

				normals[k].x = vertex.normal.x;
				normals[k].y = vertex.normal.y;
				normals[k].z = vertex.normal.z;

				tex_coords[k].u = vertex.texture_coords[0].u;
				tex_coords[k].v = vertex.texture_coords[0].v;
			} // end of vertices loop

			triangles[i].v0 = k - 3;
			triangles[i].v1 = k - 2;
			triangles[i].v2 = k - 1;
		} // end of triangles loop

		rtcCommitGeometry( mesh );
		unsigned int geom_id = rtcAttachGeometry( scene_, mesh );
		rtcReleaseGeometry( mesh );
	} // end of surfaces loop

	rtcCommitScene( scene_ );
}

Color4f Raytracer::applyShader(const int x, const int y, const float t = 0.0f) {
	RTCRayHitWithIor rtcRayHitWithIor;
	Color4f resultColor = {0 ,0 ,0, 1};

	for (int i = 0; i < this->antiAliasingSubSamplingConst; i++)
	{
		float randomX = Random();
		float randomY = Random();

		rtcRayHitWithIor.rtcRayHit.ray = camera_.GenerateRay(x + randomX, y + randomY);
		rtcRayHitWithIor.rtcRayHit.hit = createEmptyHit();
		rtcRayHitWithIor.ior = IOR_AIR;
		resultColor = resultColor + applyShaderInternal(rtcRayHitWithIor, t, 0);
	}

	resultColor = Color4f{ resultColor.r * this->antialiasingNormalizingCoef, resultColor.g * this->antialiasingNormalizingCoef, resultColor.b * this->antialiasingNormalizingCoef, 1 };
	return resultColor;
}

Color4f Raytracer::applyShaderInternal(RTCRayHitWithIor rtcRayHitWithIor, float t, int depth)
{
	Vector3 vectorFromCamera = Vector3{ rtcRayHitWithIor.rtcRayHit.ray.dir_x, rtcRayHitWithIor.rtcRayHit.ray.dir_y, rtcRayHitWithIor.rtcRayHit.ray.dir_z };
	if(depth > 5) return sphericalMap->getTexel(vectorFromCamera);

	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcIntersect1(scene_, &context, &(rtcRayHitWithIor.rtcRayHit));

	Vector3 lightPossition = Vector3{ 250, 250, 400 };
	Color4f resultColor = Color4f{ 0, 0, 0, 1};
	if (rtcRayHitWithIor.rtcRayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		IntersectionInfo intersectionInfo;

		intersectionInfo = getIntersectionInfo(rtcRayHitWithIor, vectorFromCamera, lightPossition, context);

		switch (intersectionInfo.material->shader) {
			case PHONG_SHADER: {
				resultColor = applyPhondShader(rtcRayHitWithIor, intersectionInfo, t, depth);
				break;
			}
			case GLASS_SHADER:				
				resultColor = applyGlassShader(rtcRayHitWithIor, intersectionInfo, t, depth);
				break;
			case NORMAL_SHADER:
				resultColor = applyNormalShader(rtcRayHitWithIor, intersectionInfo, t);
				break;
		}
	}
	else {
		return sphericalMap->getTexel(vectorFromCamera);
	}
	return resultColor;
}

Color4f Raytracer::applyPhondShader(RTCRayHitWithIor rtcRayHitWithIor,IntersectionInfo intersectionInfo, float t, int depth) {
	Color4f resultColor;
	float normalLigthScalarProduct = intersectionInfo.normal.DotProduct(intersectionInfo.vectorToLight);
	Vector3 lr = 2 * (normalLigthScalarProduct)* intersectionInfo.normal - intersectionInfo.vectorToLight;

	RTCRayHitWithIor phongReflectedRay = createRayWithEmptyHitAndIor(intersectionInfo.intersectionPoint, lr, FLT_MAX, 0.1f, rtcRayHitWithIor.ior);

	resultColor = Color4f{
		(intersectionInfo.material->ambient.x + intersectionInfo.enlighted * ((intersectionInfo.material->diffuse.x * normalLigthScalarProduct) + pow(intersectionInfo.material->specular.x * (-intersectionInfo.viewToIntersectionVector).DotProduct(lr), intersectionInfo.material->shininess))),
		(intersectionInfo.material->ambient.y + intersectionInfo.enlighted * ((intersectionInfo.material->diffuse.y * normalLigthScalarProduct) + pow(intersectionInfo.material->specular.y * (-intersectionInfo.viewToIntersectionVector).DotProduct(lr), intersectionInfo.material->shininess))),
		(intersectionInfo.material->ambient.z + intersectionInfo.enlighted * ((intersectionInfo.material->diffuse.z * normalLigthScalarProduct) + pow(intersectionInfo.material->specular.z * (-intersectionInfo.viewToIntersectionVector).DotProduct(lr), intersectionInfo.material->shininess))),
		1 } +(applyShaderInternal(phongReflectedRay, t, ++depth) * intersectionInfo.material->reflectivity);
	return resultColor;
}

Color4f Raytracer::applyGlassShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t, int depth) {
	Vector3 dirTowardsObs = -(intersectionInfo.viewToIntersectionVector), dirOfTransmittedRay, dirOfReflectedRay;
	Color4f resultColor = Color4f {0,0,0,1};
	float ior1, ior2;
	float cosAngle1 = intersectionInfo.normal.DotProduct(dirTowardsObs), cosAngle2 = 0;
	float Rs, Rp, R, T;

	ior1 = rtcRayHitWithIor.ior;

	ior2 = ior1 != intersectionInfo.material->ior ? intersectionInfo.material->ior : IOR_AIR;

	float sqrCos2 = 1 - sqr(ior1 / ior2) * (1 - sqr(cosAngle1));
	if (acosf(cosAngle1) > 0 && sqrCos2 > 0) {
		cosAngle2 = sqrt(sqrCos2);
		
		dirOfTransmittedRay = (ior1 / ior2) * intersectionInfo.viewToIntersectionVector + ((ior1 / ior2) * cosAngle1 - cosAngle2) * intersectionInfo.normal;
		dirOfReflectedRay = (2 * (intersectionInfo.normal.DotProduct(dirTowardsObs))) * intersectionInfo.normal - dirTowardsObs;

		Rs = sqr((ior2 * cosAngle2 - ior1 * cosAngle1) / (ior2 * cosAngle2 + ior1 * cosAngle1));
		Rp = sqr((ior2 * cosAngle1 - ior1 * cosAngle2) / (ior2 * cosAngle1 + ior1 * cosAngle2));
		R = (Rs + Rp) * 0.5f;

		T = 1 - R;

		RTCRayHitWithIor transmittedRayHitWithIor, reflectedRayHitWithIor;

		transmittedRayHitWithIor = createRayWithEmptyHitAndIor(intersectionInfo.intersectionPoint, dirOfTransmittedRay, FLT_MAX, 0.1f, ior2);
		reflectedRayHitWithIor = createRayWithEmptyHitAndIor(intersectionInfo.intersectionPoint, dirOfReflectedRay, FLT_MAX, 0.1f, ior2);

		Color4f attenuation = getAttenuationOfReflectedRay(rtcRayHitWithIor, intersectionInfo.intersectionPoint, ior2, intersectionInfo.material);

		depth++;
		resultColor = (resultColor + applyShaderInternal(transmittedRayHitWithIor, t, depth) * T) * attenuation;
		resultColor = (resultColor + applyShaderInternal(reflectedRayHitWithIor, t, depth)  * R) * attenuation;
	}
	return resultColor;
}

Color4f Raytracer::getAttenuationOfReflectedRay(RTCRayHitWithIor rtcRayHitWithIor, Vector3 intersectionPoint, float ior2, Material* material) {
	Color4f attenuation = { 0,0,0,1 };
	if (ior2 != IOR_AIR) {
		Vector3 vectorToIntersection = (intersectionPoint - Vector3{ rtcRayHitWithIor.rtcRayHit.ray.org_x, rtcRayHitWithIor.rtcRayHit.ray.org_y, rtcRayHitWithIor.rtcRayHit.ray.org_z });
		float dstToIntersection = vectorToIntersection.L2Norm();
		attenuation = Color4f{ exp(-material->diffuse.x *dstToIntersection), exp(-material->diffuse.y *dstToIntersection), exp(-material->diffuse.z *dstToIntersection) };
	}
	return attenuation;
}

Color4f Raytracer::applyNormalShader(RTCRayHitWithIor rtcRayHitWithIor, IntersectionInfo intersectionInfo, float t) {
	Color4f resultColor;
	float normalLigthScalarProduct = intersectionInfo.normal.DotProduct(intersectionInfo.vectorToLight);

	resultColor = Color4f{
		(intersectionInfo.material->diffuse.x * normalLigthScalarProduct),
		(intersectionInfo.material->diffuse.y * normalLigthScalarProduct),
		(intersectionInfo.material->diffuse.z * normalLigthScalarProduct),
		1 };
	return resultColor;
}

Color4f Raytracer::get_pixel(const int x, const int y, const float t)
{
	return applyShader(x, y, t);;
}

int Raytracer::Ui()
{
	static float f = 0.0f;
	static int counter = 0;

	// we use a Begin/End pair to created a named window
	ImGui::Begin( "Ray Tracer Params" );
	
	ImGui::Text( "Surfaces = %d", surfaces_.size() );
	ImGui::Text( "Materials = %d", materials_.size() );
	ImGui::Separator();
	ImGui::Checkbox( "Vsync", &vsync_ );
	
	//ImGui::Checkbox( "Demo Window", &show_demo_window );      // Edit bools storing our window open/close state
	//ImGui::Checkbox( "Another Window", &show_another_window );

	ImGui::SliderFloat( "float", &f, 0.0f, 1.0f );            // Edit 1 float using a slider from 0.0f to 1.0f    
	//ImGui::ColorEdit3( "clear color", ( float* )&clear_color ); // Edit 3 floats representing a color

	if ( ImGui::Button( "Button" ) )                            // Buttons return true when clicked (most widgets return true when edited/activated)
		counter++;
	ImGui::SameLine();
	ImGui::Text( "counter = %d", counter );

	ImGui::Text( "Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate );
	ImGui::End();

	// 3. Show another simple window.
	/*if ( show_another_window )
	{
	ImGui::Begin( "Another Window", &show_another_window );   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
	ImGui::Text( "Hello from another window!" );
	if ( ImGui::Button( "Close Me" ) )
	show_another_window = false;
	ImGui::End();
	}*/

	return 0;
}

Vector3 Raytracer::getInterpolatedPoint(RTCRay ray) {
	return Vector3{
			ray.org_x + ray.tfar * ray.dir_x,
			ray.org_y + ray.tfar * ray.dir_y,
			ray.org_z + ray.tfar * ray.dir_z
	};
}

float Raytracer::castShadowRay(IntersectionInfo intersectionInfo, RTCIntersectContext context) {
	RTCRay rayFromIntersectPointToLight = createRay(intersectionInfo.intersectionPoint, intersectionInfo.vectorToLight, intersectionInfo.dstToLight, 0.1f);
	rtcOccluded1(scene_, &context, &rayFromIntersectPointToLight);
	return rayFromIntersectPointToLight.tfar < intersectionInfo.dstToLight ? 0.0f : 1.0f;
}

IntersectionInfo Raytracer::getIntersectionInfo(RTCRayHitWithIor rtcRayHitWithIor, Vector3 vectorFromCamera, Vector3 lightPossition, RTCIntersectContext context) {
	IntersectionInfo intersectionInfo;
	
	intersectionInfo.intersectionPoint = Raytracer::getInterpolatedPoint(rtcRayHitWithIor.rtcRayHit.ray);

	intersectionInfo.vectorToLight = (lightPossition - intersectionInfo.intersectionPoint);
	intersectionInfo.dstToLight = (intersectionInfo.vectorToLight).L2Norm();
	(intersectionInfo.vectorToLight).Normalize();

	RTCGeometry geometry = rtcGetGeometry(scene_, rtcRayHitWithIor.rtcRayHit.hit.geomID);
	
	rtcInterpolate0(geometry, rtcRayHitWithIor.rtcRayHit.hit.primID, rtcRayHitWithIor.rtcRayHit.hit.u, rtcRayHitWithIor.rtcRayHit.hit.v,
		RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &(intersectionInfo.normal.x), 3);

	intersectionInfo.material = ((Material *)rtcGetGeometryUserData(geometry));

	if (vectorFromCamera.DotProduct(intersectionInfo.normal) > 0) {
		intersectionInfo.normal *= -1;
	}

	intersectionInfo.enlighted = castShadowRay(intersectionInfo, context);
	intersectionInfo.viewToIntersectionVector = vectorFromCamera;

	return intersectionInfo;
}