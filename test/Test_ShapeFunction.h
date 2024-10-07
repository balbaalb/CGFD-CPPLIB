#ifndef Test_ShapeFunctionH
#define Test_ShapeFunctionH
#include "../src/ShapeFunction.h"
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
bool tester_ShapeFunction(int& NumTests)
{
	if (!tester_ShapeFunction_1(NumTests))
		return false;
	if (!tester_ShapeFunction_2(NumTests))
		return false;
	++NumTests;
	return true;
}
#endif