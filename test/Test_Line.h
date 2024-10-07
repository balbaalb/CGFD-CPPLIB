#ifndef Test_LineH
#define Test_LineH
#include "../src/Line.h"
#include "../src/MathUtils.h"
bool tester_Line(int& NumTests)
{
	Vector3D P(2.0, 3.0, 5.0);
	Vector3D Q(1.0, -8.0, 10.0);
	Line L(P, Q);
	Line LL(L), LLL;
	LLL = LL;
	Vector3D PPP = LLL.GetStartPoint();
	if (PPP(0) != 2.0 || PPP(1) != 3.0 || PPP(2) != 5.0)
		return false;
	if (PPP != P)
		return false;
	if (LLL.tangent() != L.tangent())
		return false;
	double p = sqrt(1.0 + 121 + 25);
	Vector3D e(-1.0 / p, -11.0 / p, 5.0 / p);
	if (fabs(LLL.tangent()(0) - e(0)) > 1.0e-10 || fabs(LLL.tangent()(1) - e(1)) > 1.0e-10 || fabs(LLL.tangent()(2) - e(2)) > 1.0e-10)
		return false;
	if (LLL != L)
		return false;
	Line LLLL;
	LLLL.SetTwoPoints(P, Q);
	if (LLLL != L || LLLL != LL || LLLL != LLL)
		return false;
	Line L5;
	L5.SetPointAndVector(L.GetStartPoint(), L.tangent());
	if (L5 != L || !(L5 == L))
		return false;
	Line L6;
	L6.SetPointAndVector(L.GetStartPoint(), L.tangent() * (-2.5));
	if (L6 != L  || !(L6 == L))
		return false;
	Vector3D P2(1.5, 5.0, -6), Q2(-9, 4.75, -6), R2(14.78, -13.95, 0);
	Line L2, L3;
	L2.SetTwoPoints(P2, Q2);
	L3.SetPointAndVector(P2, R2);
	if (!L2.abs())
		return false;
	if (!L3.abs())
		return false;
	if (!L2.IsOnSamePlane(L3))
		return false;
	P2.SetCoords(-9.78, 12.35, -6);
	L3.SetPointAndVector(P2, R2);
	if (!L2.IsOnSamePlane(L3))
		return false;
	P2(2) = -8;
	L3.SetPointAndVector(P2, R2);
	if (L2.IsOnSamePlane(L3))
		return false;
	Vector3D P4(0, 0, -1.5), Q4(-3.5, 3.5);
	Line L4(P4, P4 + Q4);
	Vector3D A(4, 3, 1);
	Vector3D AP = L4.project(A);
	Vector3D AP_expected(0.5, -0.5, -1.5);
	if (AP != AP_expected)
		return false;
	if (fabs(AP(0) - 0.5) > 1.0e-10 || fabs(AP(1) + 0.5) > 1.0e-10 || fabs(AP(2) + 1.5) > 1.0e-10)
		return false;
	double dist4 = L4.distance(A);
	if (fabs(dist4 - sqrt(2 * 3.5 * 3.5 + 2.5 * 2.5)) > 1.0e-10)
		return false;
	Vector3D APP = AP - A;
	double dist44 = APP.abs();
	if (fabs(dist4 - dist44) > 1.0e-10)
		return false;
	if (true)
	{
		Vector3D S1(4, 3), e1(0, 0, 1);
		Line L1;
		L1.SetPointAndVector(S1, e1);
		Vector3D S2, e2(-1, 1);
		Line L2;
		L2.SetPointAndVector(S2, e2);
		Vector3D P1 = L1.NearestPointTo(L2);
		if (P1 != S1)
			return false;
		Vector3D P2 = L2.NearestPointTo(L1);
		if (fabs(P2(0) - 0.5) > 1.0e-10 || fabs(P2(1) + 0.5) > 1.0e-10 || fabs(P2(2)) > 1.0e-10)
			return false;
		double dist = L1.distance(L2);
		if (fabs(dist - sqrt(2 * 3.5 * 3.5)) > 1.0e-10)
			return false;
		if (L1.DoesIntersect(L2))
			return false;
		if (L1.intersection(L2))
			return false;
	}
	if (true)
	{
		Vector3D S1(3.5, 2.75, -10), e1(2, 1);
		Vector3D S2(3.5, 2.75, 7), e2(3, -1.6);
		Vector3D S3(3.5, 2.75, -10), e3(3, -1.6);
		Line L1, L2, L3;
		Vector3D S1_copy = S1;
		S1 = S1 + e1 * 110.21;
		L1.SetPointAndVector(S1, e1);
		S2 = S2 - e2 * 1.85;
		L2.SetPointAndVector(S2, e2);
		S3 = S3 - e3 * 1.85;
		L3.SetPointAndVector(S3, e3);
		double dist = L1.distance(L2);
		if (fabs(dist - 17) > 1.0e-10)
			return false;
		double dist23 = L1.distance(L3);
		if (fabs(dist23) > 1.0e-10)
			return false;
		Vector3D* xPtr = L1.intersection(L3);
		if (!xPtr)
			return false;
		if ((*xPtr) != S1_copy)
			return false;
		if (fabs((*xPtr)(0) - 3.5) > 1.0e-10 || fabs((*xPtr)(1) - 2.75) > 1.0e-10 || fabs((*xPtr)(2) + 10) > 1.0e-10)
			return false;
		Vector3D P1 = L1.NearestPointTo(L3);
		Vector3D P3 = L3.NearestPointTo(L1);
		if (P1 != (*xPtr) || P3 != (*xPtr))
			return false;
		double ty = pi / 12;//15 Deg
		double tz = 10.5 / 180.0 * pi;
		double tx = 47.0 / 180 * pi;
		S1 = S1.rotate_x(tx).rotate_y(ty).rotate_z(tz);
		e1 = e1.rotate_x(tx).rotate_y(ty).rotate_z(tz);
		S2 = S2.rotate_x(tx).rotate_y(ty).rotate_z(tz);
		e2 = e2.rotate_x(tx).rotate_y(ty).rotate_z(tz);
		L1.SetPointAndVector(S1, e1);
		L2.SetPointAndVector(S2, e2);
		dist = L1.distance(L2);
		if (fabs(dist - 17) > 1.0e-10)
			return false;
	}
	if (true)
	{
		Vector3D A(0, 1.5, 5.75), B(3, 2.5, 5.75);
		Line L(A, B);
		Vector3D P(4, 9.5, 5.75);
		Vector3D PP = L.project(P);
		if (fabs(PP(0) - 6) > 1.0e-10 || fabs(PP(1) - 3.5) > 1.0e-10 || fabs(PP(2) - 5.75) > 1.0e-10)
			return false;
	}
	NumTests += 1;
	return true;
}
#endif