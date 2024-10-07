#ifndef Vector3DH
#define Vector3DH
#include "Vector.h"
using namespace std;
class Vector3D : public Vector
{
protected:
	void CopyBody(const bMatrix& A);
	Vector3D CrossProduct(const Vector3D& Q) const;
public:
	Vector3D(double x = 0, double y = 0, double z = 0);
	Vector3D(const Vector& U);
	Vector3D(const bMatrix& U);
	Vector3D& operator=(const Vector3D& P);
	Vector3D operator&&(const Vector3D& Q) const;
	void SetCoords(double x = 0, double y = 0, double z = 0);
	Vector3D rotate_x(double tx) const;
	Vector3D rotate_y(double ty) const;
	Vector3D rotate_z(double tz) const;
	double CrossProduct2D(const Vector3D& Q) const;
	double GetAngleAboutXAxis() const;// a number between 0 and 2*pi
	double GetOrientationAboutXAxis() const;//a number between 0 to pi that is independent of vector direction
	static double RotationAngle_z(const Vector3D& A, const Vector3D& B, const Vector3D& C);//theta of rotation AB about B to be oriented along BC
	static double GetAngle(const Vector3D& A, const Vector3D& B, const Vector3D& C);//the ABC angle < pi
};
#endif