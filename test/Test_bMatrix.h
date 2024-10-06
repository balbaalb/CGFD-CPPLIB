#ifndef Test_bMatrixH
#define Test_bMatrixH
#include "../src/bMatrix.h"
bool tester_bMatrix(int& NumTests)
{
	bMatrix A(3, 5);
	A(0, 0) = 12; A(0, 1) = 13; A(0, 2) = -1; A(0, 3) = -5.5; A(0, 4) = -13;
	A(1, 0) = 9; A(1, 1) = 7; A(1, 2) = 6.5; A(1, 3) = -7.5; A(1, 4) = 5.3;
	A(2, 0) = 0; A(2, 1) = 15; A(2, 2) = -18; A(2, 3) = 19; A(2, 4) = -11.5;
	bMatrix AA(A);
	bMatrix AAA(10,10);
	AAA = AA;
	if (AAA.GetDim1() != 3)
		return false;
	if (AAA.GetDim2() != 5)
		return false;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			if (fabs(AAA(i, j) - A(i, j)) > 1.0e-15)
				return false;
		}
	}
	if (fabs(AAA(2, 3) - 19) > 1.0e-10)
		return false;
	bool isEqual = (AAA == A);
	if (!isEqual)
		return false;
	bool notEqual = (AAA != A);
	if (notEqual)
		return false;
	bMatrix B = A * 2.56;
	if (fabs(B(2, 3) - 19 * 2.56) > 1.0e-10)
		return false;
	bMatrix C = B + A;
	if (fabs(C(2, 3) - 19 * 3.56) > 1.0e-10)
		return false;
	bMatrix D = C - B;
	if (fabs(D(2, 3) - 19) > 1.0e-10)
		return false;
	B = B / 6.16;
	if (fabs(B(2, 3) - 19 * 2.56 / 6.16) > 1.0e-10)
		return false;
	bMatrix E(5, 2);
	E(0, 0) = 1.5; E(0, 1) = 14.5;
	E(1, 0) = -6; E(1, 1) = 8.5;
	E(2, 0) = 5.3; E(2, 1) = 1.8;
	E(3, 0) = 8.5; E(3, 1) = -9.7;
	E(4, 0) = 4.3; E(4, 1) = 11.5;
	bMatrix F = A * E;
	if (F.GetDim1() != 3)
		return false;
	if (F.GetDim2() != 2)
		return false;
	if (fabs(F(0, 0) + 167.95) > 1.0e-10)
		return false;
	if (fabs(F(0, 1) - 186.55) > 1.0e-10)
		return false;
	if (fabs(F(1, 0) + 35.01) > 1.0e-10)
		return false;
	if (fabs(F(1, 1) - 335.4) > 1.0e-10)
		return false;
	if (fabs(F(2, 0) + 73.35) > 1.0e-10)
		return false;
	if (fabs(F(2, 1) + 221.45) > 1.0e-10)
		return false;
	AAA = AAA.ChangeDim1(4);
	AAA(3, 0) = 10; AAA(3, 1) = 1.5; AAA(3, 2) = 9; AAA(3, 3) = 36; AAA(3, 4) = 12.6;
	AAA = AAA.ChangeDim2(6);
	if (AAA.GetDim1() != 4)
		return false;
	if (AAA.GetDim2() != 6)
		return false;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			if (fabs(AAA(i, j) - A(i, j)) > 1.0e-15)
				return false;
		}
	}
	if (fabs(AAA(3, 0) - 10) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 1) - 1.5) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 2) - 9) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 3) - 36) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 4) - 12.6) > 1.0e-15)
		return false;
	for (int i = 0; i < 4; ++i)
	{
		if (fabs(AAA(i, 5)) > 1.0e-15)
			return false;
	}
	AAA(0, 5) = -1;
	AAA(1, 5) = -2;
	AAA(2, 5) = -3;
	AAA(3, 5) = -4;
	AAA = AAA.insertSubMatrix(1, 4, F);
	if (AAA.GetDim1() != 7)
		return false;
	if (AAA.GetDim2() != 8)
		return false;
	for (int j = 0; j < 4; ++j)
	{
		if (fabs(AAA(0, j) - A(0, j)) > 1.0e-15)
			return false;
	}
	if (fabs(AAA(0, 6) - A(0, 4)) > 1.0e-15)
		return false;
	if (fabs(AAA(0, 7) + 1) > 1.0e-15)
		return false;
	for (int i = 1; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (fabs(AAA(i, j)) > 1.0e-15)
				return false;
		}
	}
	for (int i = 1; i < 4; ++i)
	{
		for (int j = 4; j < 6; ++j)
		{
			if (fabs(AAA(i, j) - F(i - 1, j - 4) )> 1.0e-15)
				return false;
		}
	}
	for (int i = 1; i < 4; ++i)
	{
		for (int j = 6; j < 8; ++j)
		{
			if (fabs(AAA(i, j))> 1.0e-15)
				return false;
		}
	}
	for (int i = 4; i < 6; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (fabs(AAA(i, j) - A(i - 3,j))> 1.0e-15)
				return false;
		}
	}
	if (fabs(AAA(6, 0) - 10)> 1.0e-15)
			return false;
	if (fabs(AAA(6, 1) - 1.5)> 1.0e-15)
		return false;
	if (fabs(AAA(6, 2) - 9)> 1.0e-15)
		return false;
	if (fabs(AAA(6, 3) - 36)> 1.0e-15)
		return false;
	for (int i = 4; i < 7; ++i)
	{
		for (int j = 4; j < 6; ++j)
		{
			if (fabs(AAA(i, j))> 1.0e-15)
				return false;
		}
	}
	for (int i = 4; i < 6; ++i)
	{
		if (fabs(AAA(i, 6) - A(i - 3, 4))> 1.0e-15)
			return false;
	}
	if (fabs(AAA(6, 6) - 12.6)> 1.0e-15)
		return false;
	if (fabs(AAA(4, 7) + 2)> 1.0e-15)
		return false;
	if (fabs(AAA(5, 7) + 3)> 1.0e-15)
		return false;
	if (fabs(AAA(6, 7) + 4)> 1.0e-15)
		return false;
	bMatrix AT = A.Transpose();
	if (AT.GetDim1() != 5)
		return false;
	if (AT.GetDim2() != 3)
		return false;
	if (fabs(AT(2, 1) - 6.5)> 1.0e-15)
		return false;
	if (fabs(AT(4,2) - A(2,4))> 1.0e-15)
		return false;
	if (fabs(AT(1,1) -7)> 1.0e-15)
		return false;
	if (fabs(AT(3,2) - A(2,3))> 1.0e-15)
		return false;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			if (fabs(AT(j, i) - A(i, j))> 1.0e-15)
				return false;
		}
	}
	bMatrix A2 = A.ChangeDim1(2);
	A2 = A2.ChangeDim2(4);
	if (A2.GetDim1() != 2)
		return false;
	if (A2.GetDim2() != 4)
		return false;
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (fabs(A2(i, j) - A(i, j)) > 1.0e-15)
				return false;
		}
	}
	A(0, 0) = 12; A(0, 1) = 13; A(0, 2) = -1; A(0, 3) = -5.5; A(0, 4) = -13;
	A(1, 0) = 9; A(1, 1) = 7; A(1, 2) = 6.5; A(1, 3) = -7.5; A(1, 4) = 5.3;
	A(2, 0) = 0; A(2, 1) = 15; A(2, 2) = -18; A(2, 3) = 19; A(2, 4) = -11.5;
	bMatrix A_orig = A;
	bMatrix A_swaped(3, 5);
	A_swaped(2, 0) = 12; A_swaped(2, 1) = 13; A_swaped(2, 2) = -1; A_swaped(2, 3) = -5.5; A_swaped(2, 4) = -13;
	A_swaped(1, 0) = 9; A_swaped(1, 1) = 7; A_swaped(1, 2) = 6.5; A_swaped(1, 3) = -7.5; A_swaped(1, 4) = 5.3;
	A_swaped(0, 0) = 0; A_swaped(0, 1) = 15; A_swaped(0, 2) = -18; A_swaped(0, 3) = 19; A_swaped(0, 4) = -11.5;
	A.SwitchRows(0, 2);
	if (A == A_orig)
		return false;
	if (A != A_swaped)
		return false;
	A.SwitchRows(0, 2);
	if (A != A_orig)
		return false;
	if (A == A_swaped)
		return false;
	NumTests += 1;
	return true;
}
#endif