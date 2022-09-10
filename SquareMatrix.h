#ifndef SquareMatrixH
#define SquareMatrixH
#include "bMatrix.h"
#include "Vector.h"
class SquareMatrix : public bMatrix
{
	double MakeUpperTriangular(Vector* b = 0);
	Vector Solve_GaussEliniation(const Vector& B) const;
	Vector Solve_GaussSeidel(const Vector& B) const;
	double GaussSeidelTolerance{ 1.0e-6 };
public:
	enum METHOD {GAUSS_ELIMINATION, GAUSS_SEIDEL};
	SquareMatrix(int N = 1);
	SquareMatrix(const SquareMatrix& rhs);
	SquareMatrix& operator=(const SquareMatrix& rhs);
	int GetDim() const;
	void SetGaussSeidelTolerance(double value);
	bool CanUse_GaussSeidel(bool fixRoundOffError = false);
	Vector Solve(const Vector& B, METHOD method = METHOD::GAUSS_ELIMINATION) const;
};
bool tester_SquareMatrix(int& NumTests);
bool tester_SquareMatrix_1(int& NumTests);
bool tester_SquareMatrix_2(int& NumTests);
bool tester_SquareMatrix_3(int& NumTests);
bool tester_SquareMatrix_4(int& NumTests);
bool tester_SquareMatrix_5(int& NumTests);
#endif