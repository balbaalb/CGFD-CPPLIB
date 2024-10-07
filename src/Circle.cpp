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