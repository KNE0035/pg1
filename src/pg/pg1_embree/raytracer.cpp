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

	/*if ((x == 306) && (y == 91)) {
		printf("test");
	}*/

	rtcRayHitWithIor.rtcRayHit.ray = camera_.GenerateRay(x + 0.5f, y + 0.5f);
	rtcRayHitWithIor.rtcRayHit.hit = Raytracer::createEmptyHit();
	rtcRayHitWithIor.ior = IOR_AIR;
	Color4f test = applyShaderInternal(rtcRayHitWithIor, t, 0);
	return applyShaderInternal(rtcRayHitWithIor, t, 0);
}

Color4f Raytracer::applyShaderInternal(RTCRayHitWithIor rtcRayHitWithIor, float t, int depth)
{
	Vector3 viewVector = Vector3{ rtcRayHitWithIor.rtcRayHit.ray.dir_x, rtcRayHitWithIor.rtcRayHit.ray.dir_y, rtcRayHitWithIor.rtcRayHit.ray.dir_z };
	if(depth > 5) return sphericalMap->getTexel(viewVector);

	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcIntersect1(scene_, &context, &(rtcRayHitWithIor.rtcRayHit));

	Vector3 lightPossition = Vector3{ 250, 250, 400 };
	Color4f resultColor = Color4f{ 0, 0, 0, 1};
	if (rtcRayHitWithIor.rtcRayHit.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		Vector3 vectorToLight, normal, intersectionPoint;
		float dstToLight;
		Material* material = NULL;

		getIntersectionInfo(rtcRayHitWithIor, &vectorToLight, &normal, viewVector, &intersectionPoint, lightPossition, &dstToLight, &material);

		switch (material->shader) {
			case PHONG_SHADER: {
				float normalLigthScalarProduct = normal.DotProduct(vectorToLight);
				Vector3 lr = 2 * (normalLigthScalarProduct)* normal - vectorToLight;
				float enlighted = castShadowRay(intersectionPoint, vectorToLight, dstToLight, context);

				resultColor = Color4f{
					(material->ambient.x + enlighted * ((material->diffuse.x * normalLigthScalarProduct) + pow(material->specular.x * (-viewVector).DotProduct(lr),material->shininess))),
					(material->ambient.y + enlighted * ((material->diffuse.y * normalLigthScalarProduct) + pow(material->specular.y * (-viewVector).DotProduct(lr),material->shininess))),
					(material->ambient.z + enlighted * ((material->diffuse.z * normalLigthScalarProduct) + pow(material->specular.z * (-viewVector).DotProduct(lr),material->shininess))),
					1 };
				break;
			}
			case GLASS_SHADER:				
				Vector3 dirTowardsObs = -(viewVector), dirOfTransmittedRay, dirOfReflectedRay;
				float ior1, ior2;
				float cosAngle1 = normal.DotProduct(dirTowardsObs), cosAngle2 = 0;
				float Rs, Rp, R, T;

				ior1 = rtcRayHitWithIor.ior;
			
				ior2 = ior1 != material->ior ? material->ior : IOR_AIR;

				float sqrCos2 = 1 - sqr(ior1 / ior2) * (1 - sqr(cosAngle1));
				if (acosf(cosAngle1) > 0 && sqrCos2 > 0) {
					cosAngle2 = sqrt(sqrCos2);
					dirOfTransmittedRay = (ior1 / ior2) * viewVector + ((ior1 / ior2) * cosAngle1 - cosAngle2) * normal;

					dirOfReflectedRay = (2 * (normal.DotProduct(dirTowardsObs))) * normal - dirTowardsObs;

					Rs = sqr((ior2 * cosAngle2 - ior1 * cosAngle1) / (ior2 * cosAngle2 + ior1 * cosAngle1));
					Rp = sqr((ior2 * cosAngle1 - ior1 * cosAngle2) / (ior2 * cosAngle1 + ior1 * cosAngle2));
					R = (Rs + Rp) * 0.5f;

					T = 1 - R;

					RTCRayHitWithIor transmittedRayHitWithIor, reflectedRayHitWithIor;

					transmittedRayHitWithIor.rtcRayHit.ray = Raytracer::createRay(intersectionPoint, dirOfTransmittedRay, FLT_MAX, 0.1f);
					transmittedRayHitWithIor.rtcRayHit.hit = Raytracer::createEmptyHit();
					transmittedRayHitWithIor.ior = ior2;

					reflectedRayHitWithIor.rtcRayHit.ray = Raytracer::createRay(intersectionPoint, dirOfReflectedRay, FLT_MAX, 0.1f);
					reflectedRayHitWithIor.rtcRayHit.hit = Raytracer::createEmptyHit();
					reflectedRayHitWithIor.ior = ior2;

					resultColor = resultColor + Color4f{pow(material->specular.x, material->shininess), pow(material->specular.y, material->shininess),pow(material->specular.z, material->shininess), 1 };

					float T1 = 1;
					/*if (ior2 != IOR_AIR) {
						Vector3 vectorToIntersection = (intersectionPoint - Vector3{ rtcRayHitWithIor.rtcRayHit.ray.org_x, rtcRayHitWithIor.rtcRayHit.ray.org_y, rtcRayHitWithIor.rtcRayHit.ray.org_z });
						float dstToIntersection = vectorToIntersection.L2Norm();
						T1 = exp(-0.0000001*dstToIntersection);
					}*/
					
					depth++;
					resultColor = resultColor + applyShaderInternal(transmittedRayHitWithIor, t, depth) * T * T1;
					resultColor = resultColor + applyShaderInternal(reflectedRayHitWithIor, t, depth)  * R * T1;


				}
			}
	}
	else {
		return sphericalMap->getTexel(viewVector);
	}
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

RTCRay Raytracer::createRay(Vector3 origin, Vector3 dir, float tfar, float tnear) {
	RTCRay ray;

	ray.tnear = tnear; // start of ray segment

	ray.org_x = origin.x;
	ray.org_y = origin.y;
	ray.org_z = origin.z;

	ray.dir_x = dir.x; // ray direction
	ray.dir_y = dir.y;
	ray.dir_z = dir.z;
	ray.time = 0.0f; // time of this ray for motion blur
	ray.tfar = tfar; // end of ray segment (set to hit distance)

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved
	return ray;
}

RTCHit Raytracer::createEmptyHit() {
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f; // geometry normal
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;
	return hit;
}

Vector3 Raytracer::getInterpolatedPoint(RTCRay ray) {
	return Vector3{
			ray.org_x + ray.tfar * ray.dir_x,
			ray.org_y + ray.tfar * ray.dir_y,
			ray.org_z + ray.tfar * ray.dir_z
	};
}

float Raytracer::castShadowRay(const Vector3 intersectionPoint, Vector3 vectorToLight, const float dist, RTCIntersectContext context) {
	RTCRay rayFromIntersectPointToLight = Raytracer::createRay(intersectionPoint, vectorToLight, dist, 0.1f);
	rtcOccluded1(scene_, &context, &rayFromIntersectPointToLight);
	return rayFromIntersectPointToLight.tfar < dist ? 0.0f : 1.0f;
}

void Raytracer::getIntersectionInfo(RTCRayHitWithIor rtcRayHitWithIor, Vector3* vectorToLight, Vector3* normal, Vector3 viewVector, Vector3* intersectionPoint, Vector3 lightPossition, float* dstToLight, Material** material) {
	*intersectionPoint = Raytracer::getInterpolatedPoint(rtcRayHitWithIor.rtcRayHit.ray);

	*vectorToLight = (lightPossition - (*intersectionPoint));
	*dstToLight = (*vectorToLight).L2Norm();
	(*vectorToLight).Normalize();

	RTCGeometry geometry = rtcGetGeometry(scene_, rtcRayHitWithIor.rtcRayHit.hit.geomID);
	
	rtcInterpolate0(geometry, rtcRayHitWithIor.rtcRayHit.hit.primID, rtcRayHitWithIor.rtcRayHit.hit.u, rtcRayHitWithIor.rtcRayHit.hit.v,
		RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &(normal->x), 3);

	*material = ((Material *)rtcGetGeometryUserData(geometry));

	if (viewVector.DotProduct(*normal) > 0) {
		(*normal) *= -1;
	}
}