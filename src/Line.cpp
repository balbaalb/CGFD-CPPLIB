#include "Line.h"
#include "MathUtils.h"
void Line::CopyBody(const Line& L)
{
	this->Vector3D::CopyBody(L);
	this->Start = L.Start;
}
Line::Line() : Vector3D(1, 0, 0) //x-axis is the default line
{

}
Line::Line(const Vector3D& x0, const Vector3D& x1) : Vector3D(x1 - x0)
{
	this->Start = x0;
}
Line::Line(const Line& L)
{
	this->CopyBody(L);
}
Line& Line::operator=(const Line& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
	return *this;
}
void Line::SetTwoPoints(const Vector3D& start, const Vector3D& end)
{
	this->Start = start;
	this->Vector3D::CopyBody(end - start);
}
void Line::SetPointAndVector(const Vector3D& start, const Vector3D& dx)
{
	this->Start = start;
	this->Vector3D::CopyBody(dx);
}
Vector3D Line::GetStartPoint() const
{
	return this->Start;
}
bool Line::IsOnSamePlane(const Line& L) const
{
	//BUG: When start points are the same it crashes
	Vector3D A = this->GetStartPoint();
	Vector3D B = A + (*this);
	Vector3D C = L.GetStartPoint();
	Vector3D D = C + L;
	if (A == C || A == D || B == C || B == D)
		return true;
	Vector3D AB = B - A;
	Vector3D AC = C - A;
	Vector3D AD = D - A;
	Vector3D normal = AB && AC;
	if (normal.abs() < 1.0e-10)
		return true;
	bool isPerpendicular = normal.IsPerpendicular(AD);
	return isPerpendicular;
}
Vector3D Line::project(const Vector3D& P) const
{
	Vector3D e = this->tangent();
	Vector3D P2 = P - this->Start;
	double delta = P2 || e;
	Vector3D PP = this->Start + e * delta;
	return PP;
}
double Line::distance(const Vector3D& P) const
{
	Vector3D PP = this->project(P);
	Vector3D PPP = PP - P;
	return PPP.abs();
}
Vector3D Line::NearestPointTo(const Line& L2) const
{
	//Consider a line L2, with start point S2 and tangent e2
	//Consider a point P not on L2
	//Q = P - S2
	//Delta = e2 * (Q || e2)
	//A (projection of P on L2) = S2 + Delta
	//D (vector from A to P) = Q - Delta
	//it can be shown that d2 := |D|^2 = Q^2 - (Q || e2)^2
	//Conside that P is on line L1 with start point S1 and tangent e1
	//define parameter t1 such that P = S1 + e1 * t1
	//Then d2 is a function of t1 (not shown here) and d/dt(d1) = 0 finds the parameter t1 that
	//minimizes d2, that minimized d2 is the distance between lines L1 and L2
	//And point P(t1) is the closest point to L2 on L1 ( = *this)
	Vector3D S12 = L2.GetStartPoint() - this->GetStartPoint();
	Vector3D e1 = this->tangent();
	Vector3D e2 = L2.tangent();
	double e1_dot_e2 = e1 || e2;
	double denominator = 1 - e1_dot_e2 * e1_dot_e2;
	if (!denominator)
		return this->Start;
	double numerator = S12 || (e1 - e2 * e1_dot_e2);
	double t1 = numerator / denominator;
	Vector3D P = this->Start + e1 * t1;
	return P;
}
double Line::distance(const Line& L2) const
{
	Vector3D P1 = this->NearestPointTo(L2);
	Vector3D P2 = L2.NearestPointTo(*this);
	return P1.distance(P2);
}
bool Line::DoesIntersect(const Line& L2) const
{
	if (this->IsParallel(L2))
		return false;
	if (!this->IsOnSamePlane(L2))
		return false;
	return true;
}
Vector3D* Line::intersection(const Line& L2) const
{
	Vector3D* intersectionPtr = 0;
	if (!this->DoesIntersect(L2))
		return intersectionPtr;
	Vector3D S12 = L2.GetStartPoint() - this->GetStartPoint();
	if (!S12.abs())
	{
		intersectionPtr = new Vector3D(this->GetStartPoint());
		return intersectionPtr;
	}
	Vector3D e1 = this->tangent();
	Vector3D e2 = L2.tangent();
	int i, j;
	if ((e1(0) || e2(0)) && (e1(1) || e2(1)) && fabs(S12(0)) + fabs(S12(1)))
	{
		i = 0;
		j = 1;
	}
	else if ((e1(0) || e2(0)) && (e1(2) || e2(2)) && fabs(S12(0)) + fabs(S12(2)))
	{
		i = 0;
		j = 2;
	}
	else
	{
		i = 1;
		j = 2;
	}
	double e1x = e1(i);
	double e2x = e2(i);
	double e1y = e1(j);
	double e2y = e2(j);
	double S12x = S12(i);
	double S12y = S12(j);
	double t1 = (S12x * e2y - S12y * e2x) / (e1x * e2y - e1y * e2x);
	Vector3D intersection = this->GetStartPoint() + e1 * t1;
	intersectionPtr = new Vector3D(intersection);
	return intersectionPtr;
}
bool Line::operator==(const Line& L2) const
{
	Vector3D e1 = this->tangent();
	Vector3D e2 = L2.tangent();
	if (!e1.IsParallel(e2))
		return false;
	if (L2.Start == this->Start)
		return true;
	Vector3D S1S2 = L2.Start - this->Start;
	if (!e1.IsParallel(S1S2))
		return false;
	return true;
}
bool Line::operator!=(const Line& L2) const
{
	return !((*this) == L2);
}
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