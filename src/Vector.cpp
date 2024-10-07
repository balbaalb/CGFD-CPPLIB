#include "Vector.h"
#include "MathUtils.h"
void Vector::CopyBody(const bMatrix& A)
{
	if (A.GetDim2() != 1)
		throw "Vector::copyBody(): Input has more than one col";
	this->bMatrix::CopyBody(A);
}
Vector::Vector(int M) : bMatrix(M,1)
{

}
Vector::Vector(const bMatrix& A)
{
	this->CopyBody(A);
}
Vector& Vector::operator=(const Vector& U)
{
	this->DeleteBody();
	this->CopyBody(U);
	return *this;
}
double Vector::DotProduct(const Vector& U) const
{
	if (U.GetDim1() != this->GetDim1())
		throw "Vector::DotProduct(): Input size is incorrect";
	double dotProduct = 0;
	for (int i = 0; i < this->GetDim1(); ++i)
	{
		dotProduct += (*this)(i)* U(i);
	}
	return dotProduct;
}
bMatrix Vector::DiadicProduct(const Vector& U) const
{
	return this->Multiply(U.Transpose());
}
double& Vector::operator()(int i)
{
	return this->bMatrix::Set(i, 0);
}
const double& Vector::operator()(int i) const
{
	return this->bMatrix::Get(i, 0);
}
Vector Vector::operator*(double a) const
{
	return this->Multiply(a);
}
bMatrix Vector::operator*(const Vector& B) const
{
	return this->Multiply(B.Transpose());
}
double Vector::operator||(const Vector& U) const
{
	return this->DotProduct(U);
}
double Vector::abs() const
{
	double normSquared = 0;
	for (int i = 0; i < this->GetDim1(); ++i)
	{
		normSquared += (*this)(i)* (*this)(i);
	}
	return sqrt(normSquared);
}
Vector Vector::tangent() const
{
	Vector A(*this);
	double a = this->abs();
	if (a > 0)
	{
		A = A / a;
	}
	return A;
}
bool Vector::IsParallel(const Vector& L) const
{
	double dotProduct = (*this) || L;
	double thisAbs = (*this).abs();
	double L_abs = L.abs();
	if (!L_abs || !thisAbs)
		throw "Vector::IsParallel(): Input or (*this) abs is 0";
	double cosTheta = dotProduct / thisAbs / L_abs;
	return (fabs(fabs(cosTheta) - 1) < 1.0e-10);
}
bool Vector::IsPerpendicular(const Vector& L) const
{
	if (!this->abs())
		throw "Vector::IsPerpendicular(): this->abs() = 0";
	if (!L.abs())
		throw "Vector::IsPerpendicular(): Input abs is 0";
	double dotProduct = (*this) || L;
	return (fabs(dotProduct) < 1.0e-10);
}
double Vector::distance(const Vector& L) const
{
	Vector dist = (*this) - L;
	return dist.abs();
}