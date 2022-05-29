#ifndef VectorH
#define VectorH
#include "bMatrix.h"
using namespace std;
class Vector : public bMatrix
{
protected:
	void CopyBody(const bMatrix& A);
	double DotProduct(const Vector& U) const;
	bMatrix DiadicProduct(const Vector& U) const;
public:
	Vector(int M = 1);
	Vector(const bMatrix& A);
	void operator=(const Vector& U);
	double& operator()(int i);
	const double& operator()(int i) const;
	Vector operator*(double a) const;
	bMatrix operator*(const Vector& B) const;//diadic product
	double operator||(const Vector& U) const;//dot product
	double abs() const;
	Vector tangent() const;
	bool IsParallel(const Vector& L) const;
	bool IsPerpendicular(const Vector& L) const;
	double distance(const Vector& L) const;
};
bool tester_Vector(int& NumTests);
#endif