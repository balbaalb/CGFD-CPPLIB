#ifndef Test_CircleH
#define Test_CircleH
#include "../src/Circle.h"
#include "../src/MathUtils.h"
bool tester_Circle(int& NumTests)
{
	if (true)
	{
		Circle C;
		Vector3D A(2.5, 4.0);
		C.SetCenter(A);
		double R = 5.0;
		C.SetRadius(R);
		Circle C1(C), C2;
		C2 = C1;
		if (fabs(C2.GetRadius() - R) > 1.0e-15 || C2.GetCenter() != A)
			return false;
		if (C2 != C)
			return false;
		Vector3D t1(-1.5, 3.25);
		Vector3D t2(11.5, -14.75);
		Vector3D t3(2.5, 3.85);
		t1 = t1 / t1.abs();
		t2 = t2 / t2.abs();
		t3 = t3 / t3.abs();
		Vector3D P1 = A + t1 * R;
		Vector3D P2 = A + t2 * R;
		Vector3D P3 = A + t3 * R;
		Circle C3;
		C3.Set3Points(P1, P2, P3);
		if (C3 != C)
			return false;
		if (fabs(C3.Area() - pi * R * R) > 1.0e-10)
			return false;
		if (fabs(C3.Perimeter() - 2.0 * pi * R) > 1.0e-10)
			return false;
		Circle C4(A, 12.3658);
		if (C4 == C3)
			return false;
	}
	if (true)
	{
		Vector3D A(0,4.0);
		Vector3D B, C(1.0,0);
		Line L;
		L.SetTwoPoints(B, C);
		Circle C1;
		C1.SetCenterAndTangent(A, L);
		double R = 4;
		if (fabs(C1.GetRadius() - R) > 1e-10)
			return false;
	}
	if (true)
	{
		Vector3D A(4.0, 4.0);
		Vector3D B, C(0,-1.0);
		Line L;
		L.SetTwoPoints(B, C);
		Circle C1;
		C1.SetCenterAndTangent(A, L);
		double R = 4;
		if (fabs(C1.GetRadius() - R) > 1e-10)
			return false;
	}
	if (true)
	{
		Vector3D A(4.0, 8);
		Vector3D B(0,0), C(3,1);
		Line L;
		L.SetTwoPoints(B, C);
		Circle C1;
		C1.SetCenterAndTangent(A, L);
		double R = sqrt(40);
		if (fabs(C1.GetRadius() - R) > 1e-10)
			return false;
	}
	if (true)
	{
		Vector3D A(4.0, 9.5);
		Vector3D B(0, 1.5), C(3, 2.5);
		Line L;
		L.SetTwoPoints(B, C);
		Circle C1;
		C1.SetCenterAndTangent(A, L);
		double R = sqrt(40);
		if (fabs(C1.GetRadius() - R) > 1e-10)
			return false;
		Line t1, t2;
		Vector3D D(-3, 0.5);
		int n = C1.tangent(D, &t1, &t2);
		if (n != 2)
			return false;
		if (t1 != L && t2 != L)
			return false;
	}
	if (true)
	{
		Vector3D O1(0, 1.5);
		double R1 = 4.5;
		Circle C1(O1, R1);
		Vector3D O2(9, 13.5);
		double R2 = 10;
		Circle C2(O2, R2);
		int n = C1.NumIntersectionPoints(C2);
		if (n)
			return false;
		C1.SetRadius(5);
		n = C1.NumIntersectionPoints(C2);
		if (n != 1)
			return false;
		Line L2;
		C1.SetRadius(5.5);
		n = C1.NumIntersectionPoints(C2);
		if (n != 2)
			return false;
		Vector3D Q1(11.5, -5), Q2(21.5, -5);
		R1 = 5;
		R2 = sqrt(65);
		C1.SetCenter(Q1);
		C2.SetCenter(Q2);
		C1.SetRadius(R1);
		C2.SetRadius(R2);
		Vector3D P1, P2;
		n = C1.NumIntersectionPoints(C2, &P1, &P2);
		if (n != 2)
			return false;
		if (fabs(P1(0) - 14.5) > 1.0e-10 || fabs(P1(1) + 1) > 1.0e-10)
			return false;
		if (fabs(P2(0) - 14.5) > 1.0e-10 || fabs(P2(1) + 9) > 1.0e-10)
			return false;
	}
	if (true)
	{
		Vector3D O(9, 13.5);
		double R = 10;
		Circle C(O, R);
		Vector3D P(3, 5.5);
		Vector3D t(-4, 3);
		Line L;
		L.SetPointAndVector(P, t);
		Line L2;
		int n = C.tangent(P, &L2);
		if (n != 1)
			return false;
		if (L2 != L)
			return false;
	}
	NumTests += 1;
	return true;
}
#endif