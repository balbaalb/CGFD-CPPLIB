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