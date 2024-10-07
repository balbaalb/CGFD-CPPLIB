#include "Plane.h"
#include "MathUtils.h"
void Plane::CopyBody(const Plane& P)
{
	this->Vector3D::CopyBody(P);
	this->Ref = P.Ref;
}
Plane::Plane()
{
	this->SetCoords(0, 0, 1);
}
Plane::Plane(const Plane& P)
{
	this->CopyBody(P);
}
Plane::Plane(const Vector3D& n, double D)
{
	this->SetNormalAndConst(n, D);
}
Plane::Plane(const Vector3D& n, const Vector3D& ref)
{
	this->SetNormalAndRef(n, ref);
}
Plane::Plane(const Vector3D& A, const Vector3D& B, const Vector3D& C)
{
	this->SetThreePoints(A, B, C);
}
Plane& Plane::operator=(const Plane& P)
{
	this->DeleteBody();
	this->CopyBody(P);
	return *this;
}
void Plane::SetNormalAndConst(const Vector3D& n, double D)
{
	if (!n.abs())
		throw "Plane::SetNormalAndConst(const Vector3D& n, double D) : n is zero";
	Vector3D n1 = n / n.abs();
	this->SetCoords(n1(0), n1(1), n1(2));
	double x0 = 0;
	double y0 = 0;
	double z0 = 0;
	if (n(2))
	{
		z0 = -D / n(2);
	}
	else if (n(1))
	{
		y0 = -D / n(1);
	}
	else
	{
		x0 = D / n(0);
	}
	this->Ref.SetCoords(x0, y0, z0);
}
void Plane::SetNormalAndRef(const Vector3D& n, const Vector3D& ref)
{
	if (!n.abs())
		throw "Plane::SetNormalAndRef(const Vector3D& n, const Vector3D& ref) : n is zero.";
	Vector3D n1 = n / n.abs();
	this->SetCoords(n1(0), n1(1), n1(2));
	this->Ref = ref;
}
void Plane::SetThreePoints(const Vector3D& A, const Vector3D& B, const Vector3D& C)
{
	Vector3D AB = B - A;
	Vector3D AC = C - A;
	Vector3D n = AB && AC;
	if (!n.abs())
		throw "Plane::SetThreePoints(A,B,C) : n is zero.";
	n = n / n.abs();
	this->SetCoords(n(0), n(1), n(2));
	this->Ref = A;
}
Vector3D Plane::GetRefPoint() const
{
	return this->Ref;
}
Vector3D Plane::normal() const
{
	return this->Vector3D::tangent();
}
Vector3D Plane::tangent() const//OverWriting
{
	Vector3D A, B(1,0,0);
	Vector3D AP = this->project(A);
	Vector3D BP = this->project(B);
	Vector3D AB = BP - AP;
	if (AB.abs() > 1.0e-10)
		return AB.tangent();
	Vector3D C(0, 1, 0);
	Vector3D CP = this->project(C);
	Vector3D AC = CP - AP;
	return AC.tangent();
}
double Plane::GetConst() const
{
	double D = -((*this) || this->Ref);
	return D;
}
bool Plane::IsParallel(const Line& L) const
{
	return this->Vector3D::IsPerpendicular(L);
}
bool Plane::IsPerpendicular(const Line& L) const
{
	return this->Vector3D::IsParallel(L);
}
bool Plane::IsParallel(const Plane& P2) const
{
	return this->Vector3D::IsParallel(P2);
}
bool Plane::IsPerpendicular(const Plane& P2) const
{
	return this->Vector3D::IsPerpendicular(P2);
}
bool Plane::operator==(const Plane& P) const
{
	if (!this->IsParallel(P))
		return false;
	if (!this->IsOnPlane(P.GetRefPoint()))
		return false;
	return true;
}
bool Plane::operator!=(const Plane& P) const
{
	return !this->operator==(P);
}
bool Plane::IsOnPlane(const Vector3D& Point) const
{
	Vector3D n = this->normal();
	double D = this->GetConst();
	double c = (n || Point) + D;
	if (fabs(c) > 1.0e-10)
		return false;
	else
		return true;
}
Line* Plane::intersection(const Plane& P2) const
{
	if (this->IsParallel(P2))
		return 0;
	Vector3D n1 = this->normal();
	Vector3D n2 = P2.normal();
	double A1 = n1(0);
	double B1 = n1(1);
	double C1 = n1(2);
	double A2 = n2(0);
	double B2 = n2(1);
	double C2 = n2(2);
	double D1 = -(n1 || this->Ref);
	double D2 = -(n2 || P2.Ref);
	Vector3D Start, e;
	if (fabs(A1 * B2 - A2 * B1) > 1.0e-10)
	{
		double denum = A1 * B2 - A2 * B1;
		double x0 = (-D1 * B2 + D2 * B1) / denum;
		double y0 = (-A1 * D2 + A2 * D1) / denum;
		double alpha = (-C1 * B2 + C2 * B1) / denum;
		double beta = (-A1 * C2 + A2 * C1) / denum;
		Start.SetCoords(x0, y0, 0);
		e.SetCoords(alpha, beta, 1);
	}
	else if (fabs(A1 * C2 - A2 * C1) > 1.0e-10)
	{
		double denum = A1 * C2 - A2 * C1;
		double x0 = (-D1 * C2 + D2 * C1) / denum;
		double z0 = (-A1 * D2 + A2 * D1) / denum;
		double alpha = (-B1 * C2 + B2 * C1) / denum;
		double gamma = (-A1 * B2 + A2 * B1) / denum;
		Start.SetCoords(x0, 0, z0);
		e.SetCoords(alpha, 1, gamma);
	}
	else
	{
		double denum = B1 * C2 - B2 * C1;
		double y0 = (-D1 * C2 + D2 * C1) / denum;
		double z0 = (-B1 * D2 + B2 * D1) / denum;
		double beta = (-A1 * C2 + A2 * C1) / denum;
		double gamma = (-B1 * A2 + B2 * A1) / denum;
		Start.SetCoords(0, y0, z0);
		e.SetCoords(1, beta, gamma);
	}
	Line* Lptr = new Line(Start, Start + e);
	return Lptr;
}
double Plane::distance(const Vector3D& Point) const
{
	Vector3D Q = Point - this->Ref;
	double dist = Q || this->normal();
	return dist;
}
Vector3D Plane::project(const Vector3D& Point) const
{
	Vector3D Q = Point - this->Ref;
	double dist = Q || this->normal();
	Vector3D PP = this->Ref + Q - this->normal() * dist;
	return PP;
}
Vector3D* Plane::intersection(const Line& L) const
{
	//It is obtained by minimizing the distance between a point on line L and its projection on (*this) plane.
	if (this->IsParallel(L))
		return 0;
	Vector3D n = this->normal();
	Vector3D e = L.tangent();
	Vector3D SR = this->Ref - L.GetStartPoint();
	double t = (SR || n) / (e || n);
	Vector3D A = L.GetStartPoint() + e * t;
	Vector3D* Aptr = new Vector3D(A);
	return Aptr;
}
Line* Plane::project(const Line& L) const
{
	if (this->IsPerpendicular(L))
		return 0;
	Vector3D P1 = this->project(L.GetStartPoint());
	Vector3D P2 = this->project(L.GetStartPoint() + L.tangent());
	Line* Lptr = new Line(P1, P2);
	return Lptr;
}
double Plane::distance(const Line& L) const
{
	if (!this->IsParallel(L))
		return 0;
	Vector3D S = L.GetStartPoint();
	Vector3D SP = this->project(S);
	return S.distance(SP);
}
double Plane::distance(const Plane& P) const
{
	if (!this->IsParallel(P))
		return 0;
	Vector3D R = P.Ref;
	Vector3D RP = this->project(R);
	return R.distance(RP);
}

