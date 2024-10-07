#ifndef Test_VectorH
#define Test_VectorH
#include "../src/Vector.h"
#include "../src/MathUtils.h"
bool tester_Vector(int& NumTests)
{
	Vector V(3);
	V(0) = -1; V(1) = 2.0; V(2) = 5.5;
	Vector VV(V);
	Vector VVV;
	VVV = VV;
	if (VVV(0) != -1 || VVV(1) != 2 || VVV(2) != 5.5)
		return false;
	if (VVV.GetDim1() != 3)
		return false;
	if (VVV.GetDim2() != 1)
		return false;
	bool isEqual = (VVV == V);
	if (!isEqual)
		return false;
	bool notEqual = (VVV != V);
	if (notEqual)
		return false;
	Vector W(3);
	W(0) = -18; W(1) = 1.35; W(2) = -8;
	double wv = W || V;
	if (fabs(wv - (18.0 * 1.0 + 1.35 * 2.0 - 8.0 * 5.5)) > 1.0e-15)
		return false;
	Vector U = V + W;
	if (U.GetDim1() != 3)
		return false;
	if (U.GetDim2() != 1)
		return false;
	if (fabs(U(0) + 19) > 1.0e-15)
		return false;
	if (fabs(U(1) - 3.35) > 1.0e-15)
		return false;
	if (fabs(U(2) + 2.5) > 1.0e-15)
		return false;
	U = V - W;
	if (U.GetDim1() != 3)
		return false;
	if (U.GetDim2() != 1)
		return false;
	if (fabs(U(0) -17) > 1.0e-15)
		return false;
	if (fabs(U(1) - 0.65) > 1.0e-15)
		return false;
	if (fabs(U(2) - 13.5) > 1.0e-15)
		return false;
	Vector UU = U * 5.0;
	UU = UU / 2.0;
	if (UU.GetDim1() != 3)
		return false;
	if (UU.GetDim2() != 1)
		return false;
	if (fabs(UU(0) - 17 * 2.5) > 1.0e-15)
		return false;
	if (fabs(UU(1) - 0.65 * 2.5) > 1.0e-15)
		return false;
	if (fabs(UU(2) - 13.5 * 2.5) > 1.0e-15)
		return false;
	Vector A(5);
	A(0) = 9.5;  A(1) = 11.0; A(2) = -3.0; A(3) = 4.5; A(4) = 3.14;
	bMatrix B = V * A;
	if (B.GetDim1() != 3)
		return false;
	if (B.GetDim2() != 5)
		return false;
	for (int i = 0; i < V.GetDim1(); ++i)
	{
		for (int j = 0; j < A.GetDim1(); ++j)
		{
			if (IsNotEqual(B(i, j),V(i)*A(j)))
				return false;
		}
	}
	V(0) = -1; V(1) = 2.0; V(2) = 5.5;
	double Vabs = sqrt(1.0 + 4.0 + 5.5*5.5);
	if (fabs(V.abs() - Vabs) > 1.0e-10)
		return false;
	Vector Vhat = V.tangent();
	if (Vhat.GetDim1() != 3)
		return false;
	if (Vhat.GetDim2() != 1)
		return false;
	if (fabs(Vhat(0) + 1.0 / Vabs) > 1.0e-10)
		return false;
	if (fabs(Vhat(1) - 2.0 / Vabs) > 1.0e-10)
		return false;
	if (fabs(Vhat(2) - 5.5 / Vabs) > 1.0e-10)
		return false;
	bMatrix C(3, 3);
	C(0, 0) = 11.0; C(0, 1) = -5.5; C(0, 2) = -1.56;
	C(1, 0) = 10.34; C(1, 1) = 108.85; C(1, 2) = 14.5;
	C(2, 0) = -11.75; C(2, 1) = 6.85; C(2, 2) = 5.98;
	Vector CV = C * V;
	if (CV.GetDim1() != 3)
		return false;
	if (CV.GetDim2() != 1)
		return false;
	if (fabs(CV(0) - (11.0 * (-1) - 5.5 * 2.0 - 1.56 * 5.5)) > 1.0e-10)
		return false;
	if (fabs(CV(1) - (10.34 * (-1) + 108.85 * 2.0 + 14.5 * 5.5)) > 1.0e-10)
		return false;
	if (fabs(CV(2) - (-11.75 * (-1) + 6.85 * 2.0 + 5.98 * 5.5)) > 1.0e-10)
		return false;
	Vector V3 = V * 3.0;
	bool isParallel = V3.IsParallel(V);
	if (!isParallel)
		return false;
	isParallel = V3.IsParallel(V * (-0.1569) );
	if (!isParallel)
		return false;
	V3(2) = V3(2) + 0.5;
	isParallel = V3.IsParallel(V);
	if (isParallel)
		return false;
	V(0) = -1; V(1) = 2.0; V(2) = 5.5;
	Vector VN(3);
	VN(0) = 1.0; VN(1) = 0.5;
	bool isPerpendicular = VN.IsPerpendicular(V);
	if (!isPerpendicular)
		return false;
	Vector V4 = V * 4.0;
	isPerpendicular = V4.IsPerpendicular(V);
	if (isPerpendicular)
		return false;
	W(0) = -18; W(1) = 1.35; W(2) = -8;
	if (fabs(V.distance(W) - sqrt(17 * 17 + 0.65 * 0.65 + 13.5 * 13.5)) > 1.0e-10)
		return false;
	NumTests += 1;
	return true;
}
#endif