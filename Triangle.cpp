#include "LineSegment2D.h"
#include "Triangle.h"
#include "MathUtils.h"
void Triangle::CopyBody(const Triangle& rhs)
{
	for (int i = 0; i < 3; ++i)
	{
		this->points[i] = rhs.points[i];
	}
}
Triangle::Triangle()
{

}
Triangle::Triangle(const Vector3D& A, const Vector3D& B, const Vector3D& C)
{
	this->points[0] = A;
	this->points[1] = B;
	this->points[2] = C;
}
Triangle::Triangle(const Triangle& rhs)
{
	this->CopyBody(rhs);
}
void Triangle::operator=(const Triangle& rhs)
{
	this->CopyBody(rhs);
}
double Triangle::GetArea() const
{
	Vector3D A = this->points[1] - this->points[0];
	Vector3D B = this->points[2] - this->points[1];
	return fabs(A.CrossProduct2D(B) / 2.0);
}
double Triangle::GetPerimeter() const
{
	double p = this->points[1].distance(this->points[0]);
	p += this->points[2].distance(this->points[1]);
	p += this->points[0].distance(this->points[2]);
	return p;
}
Vector3D Triangle::GetPoint(int i) const
{
	return this->points[i];
}
void Triangle::SetPoint(int i, const Vector3D& p)
{
	this->points[i] = p;
}
void Triangle::GetLineSegments(LineSegment2D& e0, LineSegment2D& e1, LineSegment2D& e2) const
{
	e0.SetPoint(0, this->points[0]);
	e1.SetPoint(0, this->points[1]);
	e2.SetPoint(0, this->points[2]);
	e0.SetPoint(1, this->points[1]);
	e1.SetPoint(1, this->points[2]);
	e2.SetPoint(1, this->points[0]);
}
bool Triangle::isInside(const Vector3D& p) const
{
	LineSegment2D AB(this->points[0], this->points[1]);
	LineSegment2D BC(this->points[1], this->points[2]);
	LineSegment2D CA(this->points[2], this->points[0]);
	if (AB.GetRelationTo(p) == BC.GetRelationTo(p) && BC.GetRelationTo(p) == CA.GetRelationTo(p))
		return true;
	return false;
}
bool Triangle::isOnOrInside(const Vector3D& p) const
{
	if (this->isInside(p))
		return true;
	LineSegment2D AB(this->points[0], this->points[1]);
	if (AB.GetRelationTo(p) == ON_LINE_SEGMENT)
		return true;
	LineSegment2D BC(this->points[1], this->points[2]);
	if (BC.GetRelationTo(p) == ON_LINE_SEGMENT)
		return true;
	LineSegment2D CA(this->points[2], this->points[0]);
	if (CA.GetRelationTo(p) == ON_LINE_SEGMENT)
		return true;
	return false;
}
bool Triangle::isIntersecting(const Triangle& T) const
{
	LineSegment2D e_this[3];
	LineSegment2D e_T[3];
	this->GetLineSegments(e_this[0], e_this[1], e_this[2]);
	T.GetLineSegments(e_T[0], e_T[1], e_T[2]);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			if (e_this[i].isIntersecting(e_T[j]) && e_this[i] != e_T[j])
				return true;
	Vector3D center_this = (this->points[0] + this->points[1] + this->points[2]) / 3.0;
	Vector3D center_T = (T.points[0] + T.points[1] + T.points[2]) / 3.0;
	if (this->isInside(T.points[0]) || this->isInside(T.points[1]) || this->isInside(T.points[2]) 
		|| this->isInside(center_T))
		return true;
	if (T.isInside(this->points[0]) || T.isInside(this->points[0]) || T.isInside(this->points[0])
		|| T.isInside(center_this))
		return true;
	return false;
}
bool Triangle::operator==(const Triangle& rhs) const
{
	if (this->points[0] == rhs.points[0])
	{
		if (this->points[1] == rhs.points[1] && this->points[2] == rhs.points[2])
			return true;
		if (this->points[1] == rhs.points[2] && this->points[2] == rhs.points[1])
			return true;
	}
	if (this->points[0] == rhs.points[1])
	{
		if (this->points[1] == rhs.points[0] && this->points[2] == rhs.points[2])
			return true;
		if (this->points[1] == rhs.points[2] && this->points[2] == rhs.points[0])
			return true;
	}
	if (this->points[0] == rhs.points[2])
	{
		if (this->points[1] == rhs.points[0] && this->points[2] == rhs.points[1])
			return true;
		if (this->points[1] == rhs.points[1] && this->points[2] == rhs.points[0])
			return true;
	}
	return false;
}
bool Triangle::operator!=(const Triangle& rhs) const
{
	return !(this->operator==(rhs));
}
bool Triangle::isCongruent(const Triangle& rhs) const
{
	double L0 = this->points[0].distance(this->points[1]);
	double L1 = this->points[1].distance(this->points[2]);
	double L2 = this->points[2].distance(this->points[0]);
	double S0 = rhs.points[0].distance(rhs.points[1]);
	double S1 = rhs.points[1].distance(rhs.points[2]);
	double S2 = rhs.points[2].distance(rhs.points[0]);
	if (IsEqual(L0, S0))
	{
		if (IsEqual(L1, S1) && IsEqual(L2, S2))
			return true;
		if (IsEqual(L1, S2) && IsEqual(L2, S1))
			return true;
	}
	if (IsEqual(L0, S1))
	{
		if (IsEqual(L1, S0) && IsEqual(L2, S2))
			return true;
		if (IsEqual(L1, S2) && IsEqual(L2, S0))
			return true;
	}
	if (IsEqual(L0, S2))
	{
		if (IsEqual(L1, S1) && IsEqual(L2, S0))
			return true;
		if (IsEqual(L1, S0) && IsEqual(L2, S1))
			return true;
	}
	return false;
}
double Triangle::GetShapeFactor() const
{
	//Based on Lo 1991 papaer
	double Area = this->GetArea();
	if (Area < 1.0e-10)
		return 0;
	LineSegment2D e0, e1, e2;
	this->GetLineSegments(e0, e1, e2);
	double L0 = e0.GetPerimeter();
	double L1 = e1.GetPerimeter();
	double L2 = e2.GetPerimeter();
	double denominator = L0 * L0 + L1 * L1 + L2 * L2;
	if (denominator < 1.0e-10)
		return 0;
	double shapeFactor = 4.0 * sqrt(3.0) * Area / denominator;
	return shapeFactor;
}
double Triangle::GetMaxAngle() const
{
	double angle1 = Vector3D::GetAngle(this->points[0], this->points[1], this->points[2]);
	double angle2 = Vector3D::GetAngle(this->points[1], this->points[2], this->points[0]);
	double angle3 = Vector3D::GetAngle(this->points[2], this->points[0], this->points[1]);
	return fmax(angle1,fmax(angle2,angle3));
}
double Triangle::GetMinAngle() const
{
	double angle1 = Vector3D::GetAngle(this->points[0], this->points[1], this->points[2]);
	double angle2 = Vector3D::GetAngle(this->points[1], this->points[2], this->points[0]);
	double angle3 = Vector3D::GetAngle(this->points[2], this->points[0], this->points[1]);
	return fmin(angle1, fmin(angle2, angle3));
}
double Triangle::GetPiAngles() const
{
	double angle1 = Vector3D::GetAngle(this->points[0], this->points[1], this->points[2]);
	double angle2 = Vector3D::GetAngle(this->points[1], this->points[2], this->points[0]);
	double angle3 = Vector3D::GetAngle(this->points[2], this->points[0], this->points[1]);
	return angle1 * angle2 * angle3;
}
double Triangle::GetIncircleRadius() const
{
	double a = this->points[0].distance(this->points[1]);
	double b = this->points[1].distance(this->points[2]);
	double c = this->points[2].distance(this->points[0]);
	double s = (a + b + c) / 2.0;
	double r = sqrt((s - a) * (s - b) * (s - c) / s);
	return r;
}
bool tester_Triangle(int& NumTests)
{
	Vector3D A, B(2, 1), C(0, 3);
	Triangle t0(A, B, B);
	t0.SetPoint(2, C);
	Triangle t1(t0), t2;
	t2 = t1;
	if (t0 != t1 || t1 != t2)
		return false;
	Vector3D D(1, 1), E(2, 0), F(1, 0.5), GG(3,1.5);//F and GG are on line AB;
	if (!t2.isInside(D))
		return false;
	if (t2.isInside(E))
		return false;
	if (t2.isInside(F))
		return false;
	if (t2.isInside(GG))
		return false;
	if (!t2.isOnOrInside(F))
		return false;
	if (t2.isOnOrInside(GG))
		return false;
	if (!t2.isOnOrInside(A) || !t2.isOnOrInside(B) || !t2.isOnOrInside(C))
		return false;
	if (!t2.isOnOrInside(A *0.5 + B * 0.5) || !t2.isOnOrInside(B * 0.1 + C * 0.9) || !t2.isOnOrInside(C * 0.6 + A * 0.4))
		return false;
	Vector3D G(0, -1), H(0,1), I(1,1), J(2,0);
	Triangle t3(A, G, B), t4(A,B,H), t5(A,B,I), t6(G,I,J);
	if (t2.isIntersecting(t3))
		return false;
	if (!t2.isIntersecting(t4))
		return false;
	if (!t2.isIntersecting(t5))
		return false;
	if (!t2.isIntersecting(t6))
		return false;
	Vector3D dispalcement(5.5, 11.5, 18);
	Vector3D A2 = A + dispalcement;
	Vector3D B2 = B + dispalcement;
	Vector3D C2 = C + dispalcement;
	double theta_z = 53 / 180 * pi;
	A2 = A2.rotate_z(theta_z);
	B2 = B2.rotate_z(theta_z);
	C2 = C2.rotate_z(theta_z);
	Triangle t7(A2, C2, B2);
	if (!t2.isCongruent(t7))
		return false;
	if (t2.isCongruent(t6))
		return false;
	double a = 1.56;
	Vector3D A3, B3(a), C3(0.5 * a, sqrt(3.0) / 2.0 * a);
	Triangle t8(A3, B3, C3);
	double sf = t8.GetShapeFactor();
	if (fabs(sf - 1) > 1.0e-10)
		return false;
	Vector3D A4, B4(a), C4(0,a);
	Triangle t9(A4, B4, C4);
	double sf2 = t9.GetShapeFactor();
	if (fabs(sf2 - sqrt(3.0)/2.0) > 1.0e-10)
		return false;
	NumTests += 1;
	Vector3D A5, B5(1), C5(0, 1);
	Triangle t10(A5, B5, C5);
	if (fabs(t10.GetMaxAngle() - pi / 2.0) > 1.0e-10)
		return false;
	Vector3D A6, B6(1), C6(-cos(pi/6.0),sin(pi/6.0));
	Triangle t11(A6, B6, C6);
	if (fabs(t11.GetMaxAngle() - 5 * pi / 6.0) > 1.0e-10)
		return false;
	return true;
}