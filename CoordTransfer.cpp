#include "CoordTransfer.h"
#include "math.h"
void CoordTransfer::CopyBody(const CoordTransfer& rhs)
{
	this->M = new bMatrix(3, 3);
	(*this->M)(0, 0) = (*this->M)(1, 1) = (*this->M)(2, 2) = 1;
	if (rhs.M)
		*this->M = *rhs.M;
	this->RevM = new bMatrix(3, 3);
	(*this->RevM)(0, 0) = (*this->RevM)(1, 1) = (*this->RevM)(2, 2) = 1;
	if (rhs.RevM)
		*(this->RevM) = *rhs.RevM;
	this->X0 = rhs.X0;
}
void CoordTransfer::DeleteBody()
{
	delete this->M;
	this->M = 0;
}
void CoordTransfer::UpdateReverseMatrix()
{
	double a = (*this->M)(0, 0); double b = (*this->M)(0, 1); double c = (*this->M)(0, 2);
	double d = (*this->M)(1, 0); double e = (*this->M)(1, 1); double f = (*this->M)(1, 2);
	double g = (*this->M)(2, 0); double h = (*this->M)(2, 1); double k = (*this->M)(2, 2);
	double det = a * (e*k - f*h) - b * (d*k - f*g) + c * (d*h - e*g);
	if (fabs(det) < 1e-10)
		throw "CoordTransfer::UpdateReverseMatrix(): determinant = 0";
	(*this->RevM)(0, 0) = (e*k - h*f) / det;  (*this->RevM)(0, 1) = -(b*k - h*c) / det; (*this->RevM)(0, 2) = (b*f - e*c) / det;
	(*this->RevM)(1, 0) = -(d*k - g*f) / det;  (*this->RevM)(1, 1) = (a*k - g*c) / det; (*this->RevM)(1, 2) = -(a*f - d*c) / det;
	(*this->RevM)(2, 0) = (d*h - g*e) / det;  (*this->RevM)(2, 1) = -(a*h - g*b) / det; (*this->RevM)(2, 2) = (a*e - d*b) / det;
}
CoordTransfer::CoordTransfer()
{
	this->M = new bMatrix(3, 3);
	(*this->M)(0, 0) = (*this->M)(1, 1) = (*this->M)(2, 2) = 1;
	this->RevM = new bMatrix(3, 3);
	(*this->RevM)(0, 0) = (*this->RevM)(1, 1) = (*this->RevM)(2, 2) = 1;
}
CoordTransfer::CoordTransfer(const CoordTransfer& rhs)
{
	this->CopyBody(rhs);
}
CoordTransfer::~CoordTransfer()
{
	this->DeleteBody();
}
void CoordTransfer::operator=(const CoordTransfer& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
}
void CoordTransfer::Translate(double x, double y, double z)
{
	this->X0(0) = x;
	this->X0(1) = y;
	this->X0(2) = z;
}
void CoordTransfer::Translate(const Vector3D& C)
{
	this->X0 = C;
}
void CoordTransfer::Rotate(double theta)
{
	(*this->M)(0, 0) = cos(theta);  (*this->M)(0, 1) = sin(theta);  (*this->M)(0, 2) = 0;
	(*this->M)(1, 0) = -sin(theta); (*this->M)(1, 1) = cos(theta);  (*this->M)(1, 2) = 0;
	(*this->M)(2, 0) = 0; 		   (*this->M)(2, 1) = 0;  		  (*this->M)(2, 2) = 1;
}
void CoordTransfer::Rotate(double dx, double dy)
{
	double cosTheta = dx / sqrt(dx*dx + dy*dy);
	double sinTheta = dy / sqrt(dx*dx + dy*dy);
	(*this->M)(0, 0) = cosTheta;  (*this->M)(0, 1) = sinTheta;  (*this->M)(0, 2) = 0;
	(*this->M)(1, 0) = -sinTheta; (*this->M)(1, 1) = cosTheta;  (*this->M)(1, 2) = 0;
	(*this->M)(2, 0) = 0; 		 (*this->M)(2, 1) = 0;  		  (*this->M)(2, 2) = 1;
}
void CoordTransfer::SetRotationMatrix(const bMatrix& A)
{
	if (!this->M)
		this->M = new bMatrix(3, 3);
	*(this->M) = A;
	this->UpdateReverseMatrix();
}
void CoordTransfer::OriginalToNew(double x, double y, double& x1, double& y1) const
{
	Vector3D p0(x, y);
	Vector3D p1 = this->OriginalToNew(p0);
	x1 = p1(0);
	y1 = p1(1);
}
void CoordTransfer::NewToOriginal(double x1, double y1, double& x, double& y) const
{
	Vector3D p1(x1, y1);
	Vector3D p0 = this->NewToOriginal(p1);
	x = p0(0);
	y = p0(1);
}
Vector3D CoordTransfer::OriginalToNew(const Vector3D& p0) const
{
	return (*this->M) * p0 + this->X0;
}
Vector3D CoordTransfer::NewToOriginal(const Vector3D& p1) const
{
	return (*this->RevM) * (p1 - this->X0);
}
bool tester_CoordTransfer(int& NumTests)
{
	bMatrix M(3, 3);
	M(0, 0) = 11.2; M(0, 1) = 6.8;
	M(1, 0) = -14.5; M(1, 1) = -5.9;
	M(2, 2) = 1.0;
	Vector3D C(12.5, -9.8);
	CoordTransfer CTr;
	CTr.SetRotationMatrix(M);
	CTr.Translate(C);
	CoordTransfer CTr1(CTr), CTr2;
	CTr2 = CTr1;
	Vector3D p0(12.5, -101);
	Vector3D p1 = M * p0 + C;
	Vector3D q1 = CTr2.OriginalToNew(p0);
	if (p1 != q1)
		return false;
	Vector3D q0 = CTr2.NewToOriginal(q1);
	if (q0 != p0)
		return false;
	NumTests += 1;
	return true;
}