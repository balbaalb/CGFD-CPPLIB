#include "Vector3D.h"
#include "MathUtils.h"
void Vector3D::CopyBody(const bMatrix& A)
{
	if (A.GetDim1() != 3)
		throw "Vector3D::CopyBody(): Input dim1 != 3";
	this->Vector::CopyBody(A);
}
Vector3D::Vector3D(double x, double y, double z) : Vector(3)
{
	(*this)(0) = x;
	(*this)(1) = y;
	(*this)(2) = z;
}
Vector3D::Vector3D(const Vector& U)
{
	this->CopyBody(U);
}
Vector3D::Vector3D(const bMatrix& U)
{
	this->CopyBody(U);
}
Vector3D& Vector3D::operator=(const Vector3D& P)
{
	this->DeleteBody();
	this->CopyBody(P);
	return *this;
}
Vector3D Vector3D::CrossProduct(const Vector3D& Q) const
{
	double ax = (*this)(1) * Q(2) - (*this)(2) * Q(1);
	double ay = (*this)(2) * Q(0) - (*this)(0) * Q(2);
	double az = (*this)(0) * Q(1) - (*this)(1) * Q(0);
	Vector3D A(ax, ay, az);
	return A;
}
Vector3D Vector3D::operator&&(const Vector3D& Q) const
{
	return this->CrossProduct(Q);
}
void Vector3D::SetCoords(double x, double y, double z)
{
	(*this)(0) = x;
	(*this)(1) = y;
	(*this)(2) = z;
}
Vector3D Vector3D::rotate_x(double tx) const
{
	bMatrix R(3, 3);
	R(0, 0) = 1;
	R(1, 1) = cos(tx); R(1, 2) = -sin(tx);
	R(2, 1) = sin(tx); R(2, 2) = cos(tx);
	Vector3D A = R * (*this);
	return A;
}
Vector3D Vector3D::rotate_y(double ty) const
{
	bMatrix R(3, 3);
	R(1, 1) = 1;
	R(0, 0) = cos(ty); R(0, 2) = -sin(ty);
	R(2, 0) = sin(ty); R(2, 2) = cos(ty);
	Vector3D A = R * (*this);
	return A;
}
Vector3D Vector3D::rotate_z(double tz) const
{
	bMatrix R(3, 3);
	R(2, 2) = 1;
	R(0, 0) = cos(tz); R(0, 1) = -sin(tz);
	R(1, 0) = sin(tz); R(1, 1) = cos(tz);
	Vector3D A = R * (*this);
	return A;
}
double Vector3D::CrossProduct2D(const Vector3D& Q) const
{
	return ((*this)(0) * Q(1) - (*this)(1) * Q(0));
}
double Vector3D::GetAngleAboutXAxis() const
{
	double x = (*this)(0);
	double y = (*this)(1);
	if (fabs(x) < 1.0e-10)
		return (y >= 0 ? pi / 2.0 : 3.0 * pi / 2.0);
	double theta = fabs(atan(fabs(y) / fabs(x)));
	if (x >= 0 && y >= 0)
		return theta;
	if (x < 0 && y >= 0)
		return pi - theta;
	if (x < 0 && y < 0)
		return pi + theta;
	return 2.0 * pi - theta;
}
double Vector3D::GetOrientationAboutXAxis() const
{
	double theta = this->GetAngleAboutXAxis();
	return theta <= pi - 1.0e-10 ? theta : theta - pi;
}
double Vector3D::RotationAngle_z(const Vector3D& A, const Vector3D& B, const Vector3D& C)//static
{//theta of rotation AB about B to be oriented along BC
	Vector3D AB = B - A;
	Vector3D BC = C - B;
	if (AB.abs() < 1.0e-10 || BC.abs() < 1.0e-10)
		throw "Vector3D::RotationAngle(): inputs are coincidental";
	AB = AB / AB.abs();
	BC = BC / BC.abs();
	double sinTheta = AB.CrossProduct2D(BC);
	double cosTheta = AB || BC;
	double theta = asin(sinTheta);
	if (cosTheta < 0){
		if (theta > 0)
			theta = pi - theta;
		else
			theta = -pi - theta;
	}
	return theta;
}
double Vector3D::GetAngle(const Vector3D& A, const Vector3D& B, const Vector3D& C)
{
	Vector3D BA = A - B;
	Vector3D BC = C - B;
	double ba = BA.abs();
	double bc = BC.abs();
	if (ba < 1.0e-15 || bc < 1.0e-15)
		throw "Vector3D::RotationAngle(): inputs are coincidental";
	double dotProd = (BA || BC);
	double cosTheta =  dotProd /ba / bc;
	double theta = acos(cosTheta);
	if (theta < 0)
		theta += pi;
	return theta;
}