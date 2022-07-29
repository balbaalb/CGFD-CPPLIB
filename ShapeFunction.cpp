#include <math.h>
#include "ShapeFunction.h"
void ShapeFunction::CopyBody(const ShapeFunction& rhs)
{
	this->type = rhs.type;
	this->N = 0;
	if (rhs.N)
		this->N = new SquareMatrix(*rhs.N);
	this->Xmax = rhs.Xmax;
	this->uh = rhs.uh;
	this->nh = rhs.nh;
	this->Centeriod = rhs.Centeriod;
	this->Pe = rhs.Pe;
}
void ShapeFunction::DeleteBody()
{
	delete this->N;
	this->N = 0;
}
ShapeFunction::ShapeFunction(const Triangle& T)
{
	this->type = SHAPE_FUNCTION_TYPE::LINEAR;
	this->Pe = 0;
	this->Xmax = 0;
	this->N = new SquareMatrix(3);
	SquareMatrix K(3);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			K(i, j) = (j != 2) ? T.GetPoint(i)(j) : 1;
		}
	}
	for (int i = 0; i < 3; ++i){
		Vector C(3);
		for (int j = 0; j < 3; ++j){
			C(j) = (i == j) ? 1 : 0;
		}
		Vector coeffs = K.Solve(C);
		for (int j = 0; j < 3; ++j){
			(*this->N)(i, j) = coeffs(j);
		}
	}
}
ShapeFunction::ShapeFunction(const Triangle& T, double Pe_, const Vector3D& uh_)
{
	this->type = SHAPE_FUNCTION_TYPE::EXPONENTIAL;
	this->Pe = Pe_;
	this->uh = uh_ / uh_.abs();
	this->Xmax = 0;
	this->N = new SquareMatrix(3);
	SquareMatrix K(3);
	double X[3], Y[3];
	Vector3D ez(0, 0, 1);
	this->nh = ez && uh;
	nh(2) = 0;
	this->Centeriod = (T.GetPoint(0) + T.GetPoint(1) + T.GetPoint(2)) / 3.0;
	for (int i = 0; i < 3; ++i)
	{
		Vector3D r = T.GetPoint(i) - this->Centeriod;
		X[i] = r || uh;
		if (X[i] > Xmax)
			Xmax = X[i];
		Y[i] = r || nh;
	}
	for (int i = 0; i < 3; ++i){
		K(i, 0) = exp(Pe*(X[i] - Xmax));
		K(i, 1) = Y[i];
		K(i, 2) = 1;
	}
	for (int i = 0; i < 3; ++i){
		Vector C(3);
		for (int j = 0; j < 3; ++j){
			C(j) = (i == j) ? 1 : 0;
		}
		Vector coeffs = K.Solve(C);
		for (int j = 0; j < 3; ++j){
			(*N)(i, j) = coeffs(j);
		}
	}
}
ShapeFunction::ShapeFunction(const Triangle& T, const Vector3D& uh_)
{
	this->type = SHAPE_FUNCTION_TYPE::EXPONENTIAL;
	this->Pe = 0;
	this->Xmax = 0;
	this->uh = uh_ / uh_.abs();
	this->N = new SquareMatrix(3);
	SquareMatrix K(2);
	double X[3], Y[3];
	Vector3D ez(0, 0, 1);
	this->nh = ez && uh;
	nh(2) = 0;
	this->Centeriod = (T.GetPoint(0) + T.GetPoint(1) + T.GetPoint(2)) / 3.0;
	for (int i = 0; i < 3; ++i){
		Vector3D r = T.GetPoint(i) - this->Centeriod;
		Y[i] = r || nh;
	}
	int i = 0;//the node that will not participate in the shape function
	int j = 1;
	int k = 2;
	if ((Y[1] - Y[0]) * (Y[1] - Y[2]) <= 0){
		i = 1; 		j = 2; 		k = 0;
	}
	else if ((Y[2] - Y[0]) * (Y[2] - Y[1]) <= 0){
		i = 2;		j = 0;		k = 1;
	}
	K(0, 0) = Y[j];	K(0, 1) = 1.0;
	K(1, 0) = Y[k];	K(1, 1) = 1.0;
	Vector C(2);
	C(0) = 1;	C(1) = 0;
	Vector coeffs_j = K.Solve(C);
	(*N)(j, 1) = coeffs_j(0);
	(*N)(j, 2) = coeffs_j(1);
	C(0) = 0;	C(1) = 1;
	Vector coeffs_k = K.Solve(C);
	(*N)(k, 1) = coeffs_k(0);
	(*N)(k, 2) = coeffs_k(1);
}
ShapeFunction::ShapeFunction(const ShapeFunction& rhs)
{
	this->CopyBody(rhs);
}
ShapeFunction::~ShapeFunction()
{
	this->DeleteBody();
}
void ShapeFunction::operator=(const ShapeFunction& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
}
double ShapeFunction::GetValue(int i, const Vector3D& P) const
{
	if (this->type == SHAPE_FUNCTION_TYPE::EXPONENTIAL)
	{
		Vector3D r = P - this->Centeriod;
		double X = r || this->uh;
		double Y = r || this->nh;
		double Z = exp(this->Pe * (X - this->Xmax));//different than the Z defined in Baliga & Patankar (1985)
		return ((*this->N)(i, 0) * Z + (*this->N)(i, 1) * Y + (*this->N)(i, 2));
	}
	return ((*this->N)(i, 0) * P(0) + (*this->N)(i, 1) * P(1) + (*this->N)(i, 2));
}
Vector3D ShapeFunction::Grad(int i, const Vector3D& P) const
{
	if (this->type == SHAPE_FUNCTION_TYPE::EXPONENTIAL)
	{
		Vector3D r = P - this->Centeriod;
		double X = r || this->uh;
		double Y = r || this->nh;
		double Z = exp(this->Pe * (X - this->Xmax));//different than the Z defined in Baliga & Patankar (1985)
		double Nx = (*this->N)(i, 0) * Z * Pe * uh(0) + (*this->N)(i, 1) * nh(0);
		double Ny = (*this->N)(i, 0) * Z * Pe * uh(1) + (*this->N)(i, 1) * nh(1);
		Vector3D DeltaN(Nx, Ny);
		return DeltaN;
	}
	Vector3D DeltaN((*this->N)(i, 0), (*this->N)(i, 1));
	return DeltaN;
}
bool tester_ShapeFunction(int& NumTests)
{
	if (!tester_ShapeFunction_1(NumTests))
		return false;
	if (!tester_ShapeFunction_2(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_ShapeFunction_1(int& NumTests)
{
	Vector3D A, B(2, 1), C(0, 3);
	Triangle T(A, B, C);
	ShapeFunction N(T);
	if (fabs(N.GetValue(0, A) - 1) > 1.0e-10 && fabs(N.GetValue(0, B)) > 1.0e-10 && fabs(N.GetValue(0, C)) > 1.0e-10)
		return false;
	if (fabs(N.GetValue(1, A)) > 1.0e-10 && fabs(N.GetValue(1, B) - 1) > 1.0e-10 && fabs(N.GetValue(1, C)) > 1.0e-10)
		return false;
	if (fabs(N.GetValue(2, A)) > 1.0e-10 && fabs(N.GetValue(2, B)) > 1.0e-10 && fabs(N.GetValue(2, C) - 1) > 1.0e-10)
		return false;
	++NumTests;
	return true;
}
bool tester_ShapeFunction_2(int& NumTests)
{
	Vector3D A, B(1), C(0, 1);
	Triangle T(A, B, C);
	Vector3D uh(1,1);//perpendicular to BC
	double Pe = 2;
	ShapeFunction N(T,Pe,uh);
	if (fabs(N.GetValue(0, A) - 1) > 1.0e-10 && fabs(N.GetValue(0, B)) > 1.0e-10 && fabs(N.GetValue(0, C)) > 1.0e-10)
		return false;
	if (fabs(N.GetValue(1, A)) > 1.0e-10 && fabs(N.GetValue(1, B) - 1) > 1.0e-10 && fabs(N.GetValue(1, C)) > 1.0e-10)
		return false;
	if (fabs(N.GetValue(2, A)) > 1.0e-10 && fabs(N.GetValue(2, B)) > 1.0e-10 && fabs(N.GetValue(2, C) - 1) > 1.0e-10)
		return false;
	Vector3D M = (A + B + C) / 3.0;
	double N0 = N.GetValue(0, M);
	double N1 = N.GetValue(1, M);
	double N2 = N.GetValue(2, M);
	if (N.GetValue(0, M) < 1.0 / 3.0 || N.GetValue(1, M) > 1.0 / 3.0 || N.GetValue(2, M) > 1.0 / 3.0)
		return false;
	++NumTests;
	return true; 
}