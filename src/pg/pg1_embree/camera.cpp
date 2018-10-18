#include "stdafx.h"
#include "camera.h"
#include "raytracer.h"

Camera::Camera( const int width, const int height, const float fov_y,
	const Vector3 view_from, const Vector3 view_at )
{
	width_ = width;
	height_ = height;
	fov_y_ = fov_y;

	view_from_ = view_from;
	view_at_ = view_at;

	// TODO compute focal lenght based on the vertical field of view and the camera resolution
	f_y_ = height  / (2 * tanf(fov_y * 0.5f));
	
	// TODO build M_c_w_ matrix	
	
	Vector3 basis_x, basis_y, basis_z;

	basis_z = view_from - view_at;
	basis_z.Normalize();

	basis_x = up_.CrossProduct(basis_z);
	basis_x.Normalize();
	
	basis_y = basis_z.CrossProduct(basis_x);
	basis_y.Normalize();


	M_c_w_ = Matrix3x3( basis_x, basis_y, basis_z );
}

RTCRay Camera::GenerateRay( const float x_i, const float y_i ) const
{
	Vector3 d_c = { x_i - width_ * 0.5f, height_ * 0.5f - y_i, - f_y_  };
	Vector3 d_w = M_c_w_ * d_c; 
	d_w.Normalize();

	return createRay(view_from_, d_w);
}
