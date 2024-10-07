#ifndef Test_SquareMatrixH
#define Test_SquareMatrixH
#include "../src/SquareMatrix.h"
bool tester_SquareMatrix_1(int& NumTests)
{
	SquareMatrix A(5);
	A(1, 1) = 4; 	A(1, 2) = 5; 	A(1, 3) = 12; 	A(1, 4) = 13.5;
	A(1, 1) = 0.2; 	A(1, 2) = 4;	A(2, 3) = 0.5; 	A(2, 4) = -17;
	A(3, 1) = 6;  	A(3, 2) = 14;	A(3, 3) = 13; 	A(3, 4) = -18;
	A(4, 1) = 5; 	A(4, 2) = 7; 	A(4, 3) = -9;	A(4, 4) = -9;
	Vector b(5);
	b(1) = 11;
	b(2) = -3;
	b(3) = 0;
	b(4) = 5;
	Vector X = A.Solve(b);
	Vector bb = A * X;
	if (bb != b)
		return false;
	for (int i = 0; i < 5; ++i){
		if (fabs(bb(i) - b(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_2(int& NumTests)
{
	//solution of A*xx=yy : Needs correct pivoting
	SquareMatrix A(4);
	A(1, 1) = 1.0; A(1, 2) = 0.0,	A(1, 3) = 0.0;
	A(2, 1) = 1.0; A(2, 2) = 0.0;	A(2, 3) = 1.0;
	A(3, 1) = 1.0; A(3, 2) = 1.0;	A(3, 3) = 0.0;
	Vector x(4); 
	x(1) = 1.0;
	x(2) = 0.5;
	x(3) = -2.0;
	Vector yy = A * x;
	Vector xx = A.Solve(yy);
	if (xx != x)
		return false;
	for (int i = 0; i < 4; ++i){
		if (fabs(xx(i) - x(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_3(int& NumTests)
{
	
	SquareMatrix A6(7);
	A6(1, 1) = 1; A6(1, 2) = 0; A6(1, 3) = 0; A6(1, 4) = 0; A6(1, 5) = 0; A6(1, 6) = 0;
	A6(2, 1) = 1; A6(2, 2) = 1; A6(2, 3) = 1; A6(2, 4) = 0; A6(2, 5) = 0; A6(2, 6) = 0;
	A6(3, 1) = 1; A6(3, 2) = 0; A6(3, 3) = 0; A6(3, 4) = 0; A6(3, 5) = 1; A6(3, 6) = 0;
	A6(4, 1) = 1; A6(4, 2) = 1; A6(4, 3) = 1; A6(4, 4) = 1; A6(4, 5) = 1; A6(4, 6) = 1;
	A6(5, 1) = 1; A6(5, 2) = 1; A6(5, 3) = 1; A6(5, 4) = 2; A6(5, 5) = 4; A6(5, 6) = 4;
	A6(6, 1) = 1; A6(6, 2) = 2; A6(6, 3) = 8; A6(6, 4) = 4; A6(6, 5) = 1; A6(6, 6) = 2;
	Vector x(7);
	x(1) = 0.5;
	x(2) = 1.6;
	x(3) = -1;
	x(4) = 2;
	x(5) = 3;
	x(6) = 1;
	Vector yy = A6 * x;
	Vector xx = A6.Solve(yy);
	if (xx != x)
		return false;
	for (int i = 0; i < 7; ++i){
		if (fabs(xx(i) - x(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_4(int& NumTests)
{
	
	SquareMatrix A(4);
	A(1, 1) = 5.0; A(1, 2) = 0.0, A(1, 3) = 0.0;
	A(2, 1) = -1.0; A(2, 2) = 0.0;	A(2, 3) = 16.0;
	A(3, 1) = 4.0; A(3, 2) = -6.0;	A(3, 3) = 0.0;
	Vector x(4);
	x(1) = 1.0;
	x(2) = 0.5;
	x(3) = -2.0;
	Vector yy = A * x;
	Vector xx = A.Solve(yy);
	if (xx != x)
		return false;
	for (int i = 0; i < 4; ++i){
		if (fabs(xx(i) - x(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_5(int& NumTests)
{

	SquareMatrix A(4);
	Vector B(4), X(4);
	A(0, 0) = 10;	A(0, 1) = 2;	A(0, 2) = -1;	A(0, 3) = 5;	X(0) = 5;	B(0) = 64.5;
	A(1, 0) = 5;	A(1, 1) = 12;	A(1, 2) = 2;	A(1, 3) = 4;	X(1) = -11; B(1) = -61;
	A(2, 0) = -3;	A(2, 1) = 5;	A(2, 2) = 18; 	A(2, 3) = 8;	X(2) = 6;	B(2) = 106;
	A(3, 0) = -2;	A(3, 1) = 4;	A(3, 2) = 3;	A(3, 3) = 20;	X(3) = 8.5;	B(3) = 134;
	if (!A.CanUse_GaussSeidel())
		return false;
	Vector xx = A.Solve(B, SquareMatrix::METHOD::GAUSS_SEIDEL);
	for (int i = 0; i < 4; ++i) {
		if (fabs(xx(i) - X(i)) > 1.0e-5)
			return false;
	}
	A(2, 2) = 15.99;
	if (A.CanUse_GaussSeidel())
		return false;
	++NumTests;
	return true;
}
bool tester_SquareMatrix(int& NumTests)
{
	SquareMatrix A(2);
	A(0, 0) = 11.5;	A(0, 1) = -2.35;
	A(1, 0) = 17;
	SquareMatrix A2(A), A3;
	A3 = A2;
	if (A3.GetDim1() != 2 || A3.GetDim1() != 2)
		return false;
	if (A3 != A)
		return false;
	if (!tester_SquareMatrix_1(NumTests))
		return false;
	if (!tester_SquareMatrix_2(NumTests))
		return false;
	if (!tester_SquareMatrix_3(NumTests))
		return false;
	if (!tester_SquareMatrix_4(NumTests))
		return false;
	if (!tester_SquareMatrix_5(NumTests))
		return false;
	++NumTests;
	return true;
}
#endif