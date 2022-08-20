#ifndef SquareMatrixH
#define SquareMatrixH
#include "bMatrix.h"
#include "Vector.h"
class SquareMatrix : public bMatrix
{
public:
	SquareMatrix(int N = 1);
	SquareMatrix(const SquareMatrix& rhs);
	SquareMatrix& operator=(const SquareMatrix& rhs);
	int GetDim() const;
	double MakeUpperTriangular(Vector* b = 0);
	Vector Solve(const Vector& B) const;
};
bool tester_SquareMatrix(int& NumTests);
bool tester_SquareMatrix_1(int& NumTests);
bool tester_SquareMatrix_2(int& NumTests);
bool tester_SquareMatrix_3(int& NumTests);
bool tester_SquareMatrix_4(int& NumTests);
#endif