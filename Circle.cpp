#include "Circle.h"
#include "Line.h"
#include "MathUtils.h"
void Circle::CopyBody(const Circle& c)
{
	this->center = c.center;
	this->radius = c.radius;
}
Circle::Circle()
{
	this->radius = 0;
}
Circle::Circle(Vector3D Center, double rad)
{
	this->center = Center;
	this->radius = rad;
}
Circle::Circle(const Circle& c)
{
	this->CopyBody(c);
}
Circle& Circle::operator=(const Circle& c)
{
	this->CopyBody(c);
	return *this;
}
double Circle::Area() const
{
	return pi * this->radius * this->radius;
}
double Circle::Perimeter() const
{
	return 2.0 * pi * this->radius;
}
void Circle::SetCenter(const Vector3D& Center)
{
	this->center = Center;
}
Vector3D Circle::GetCenter() const
{
	return this->center;
}
void Circle::SetRadius(double R)
{
	this->radius = R;
}
double Circle::GetRadius() const
{
	return this->radius;
}
void Circle::Set3Points(const Vector3D& P1, const Vector3D& P2, const Vector3D& P3)
{
	Vector3D A = (P1 + P2) / 2.0;
	Vector3D B = (P2 + P3) / 2.0;
	Line L12(P1, P2);
	Line L23(P2, P3);
	Vector3D e3(0, 0, 1);
	Vector3D L12n = e3 && L12.tangent();
	Vector3D L23n = e3 && L23.tangent();
	Line LA, LB;
	LA.SetPointAndVector(A, L12n);
	LB.SetPointAndVector(B, L23n);
	Vector3D* Center = LA.intersection(LB);
	if (!Center)
		return;
	this->center = *Center;
	delete Center;
	this->radius = this->center.distance(P1);
}
bool Circle::operator==(const Circle& c) const
{
	return (fabs(this->radius - c.radius) < 1.0e-15 && this->center == c.center);
}
bool Circle::operator!=(const Circle& c) const
{
	return !((*this) == c);
}
int Circle::NumIntersectionPoints(const Circle& c, Vector3D* P1, Vector3D* P2) const
{
	double d = this->center.distance(c.center);
	if (d > this->radius + c.radius)
		return 0;
	if (d == this->radius + c.radius)
	{
		if (P1)
		{
			double alpha = this->radius / (this->radius + c.radius);
			Vector3D A = this->center * alpha + c.center * (1 - alpha);
			*P1 = A;
		}
		return 1;
	}
	if (P1 || P2)
	{
		double r = this->radius;
		double R = c.radius;
		double cosAlpha = (r * r + d * d - R * R) / (2.0 * r * d);//alpha = angle between r and d
		double dr = r * cosAlpha;
		double h = sqrt(r * r - dr * dr);
		Line centerLine;
		centerLine.SetTwoPoints(this->center, c.center);
		Vector3D e3(0, 0, 1);
		Vector3D normal = e3 && centerLine.tangent();
		Vector3D A = this->center + centerLine.tangent() * dr + normal * h;
		Vector3D B = this->center + centerLine.tangent() * dr - normal * h;
		if(P1)
			*P1 = A;
		if (P2)
			*P2 = B;
	}
	return 2;
}
void Circle::SetCenterAndTangent(const Vector3D& Center, const Line& tangent)
{
	this->center = Center;
	this->radius = tangent.distance(Center);
}
int Circle::tangent(const Vector3D& point, Line* t1, Line* t2) const
{
	double d = this->center.distance(point);
	if (d < this->radius)
		return 0;
	if (d == this->radius)
	{
		if (t1)
		{
			Line L;
			L.SetTwoPoints(point, this->center);
			Vector3D e3(0, 0, 1);
			Vector3D n = e3 && L.tangent();
			t1->SetPointAndVector(point, n);
		}
		return 1;
	}
	if (t1 || t2)
	{
		double S = sqrt(d * d - this->radius* this->radius);
		Circle C2(point, S);
		Vector3D A, B;
		this->NumIntersectionPoints(C2, &A, &B);
		if (t1)
			t1->SetTwoPoints(point, A);
		if (t2)
			t2->SetTwoPoints(point, B);
	}
	return 2;
}
bool Circle::isInside(const Vector3D& p, double tolerance) const
{//if on the circle will return false
	double dist = this->center.distance(p);
	return (dist < this->radius - tolerance);
}
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