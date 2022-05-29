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
void Plane::operator=(const Plane& P)
{
	this->DeleteBody();
	this->CopyBody(P);
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
bool tester_Plane(int& NumTests)
{
	if (true)
	{
		Vector3D n(2, 5, 8);
		Vector3D R(0, 0, 8);
		Plane P(n, R);
		Plane PP(P);
		Plane PPP;
		PPP = PP;
		Vector3D nnn = PPP.normal();
		Vector3D RRR = PPP.GetRefPoint();
		if (!nnn.IsParallel(n) || RRR != R)
			return false;
		if (!PPP.IsOnPlane(R))
			return false;
	}
	if (true)
	{
		Vector3D A(2, 0, 0);
		Vector3D B(0, 5, 0);
		Vector3D C(0, 0, 8);
		Plane P(A, B, C);
		if (!P.IsOnPlane(A) || !P.IsOnPlane(B) || !P.IsOnPlane(C))
			return false;
		Vector3D AA = A, O; AA(2) = 1;
		if (P.IsOnPlane(AA) || P.IsOnPlane(O))
			return false;
		Line e1(O, A);
		Line e2(O, B);
		Line e3(O, C);
		Vector3D* AAA = P.intersection(e1);
		Vector3D* BBB = P.intersection(e2);
		Vector3D* CCC = P.intersection(e3);
		if (!AAA || !BBB || !CCC)
			return false;
		if ((*AAA) != A || (*BBB) != B || (*CCC) != C)
			return false;
		if (fabs((*AAA)(0) - 2) > 1.0e-10 || fabs((*AAA)(1)) > 1.0e-10 || fabs((*AAA)(2)) > 1.0e-10)
			return false;
		if (fabs((*BBB)(0)) > 1.0e-10 || fabs((*BBB)(1) - 5) > 1.0e-10 || fabs((*BBB)(2)) > 1.0e-10)
			return false;
		if (fabs((*CCC)(0)) > 1.0e-10 || fabs((*CCC)(1)) > 1.0e-10 || fabs((*CCC)(2) - 8) > 1.0e-10)
			return false;
		Line AB(A, B);
		Line AC(A, C);
		Line BC(B, C);
		Vector3D A_B = B - A;
		Vector3D A_C = C - A;
		Vector3D B_C = C - B;
		Plane xy(O, A, B);
		Plane xz(e2, O);
		Plane yz(e1, 0.0);
		Line* L1 = P.intersection(yz);
		Line* L2 = P.intersection(xz);
		Line* L3 = P.intersection(xy);
		Vector3D n = P.normal();
		if (fabs(n.abs() - 1) > 1.0e-10)
			return false;
		if (fabs(AB.abs() - sqrt(4+25)) > 1.0e-10 || fabs(BC.abs() - sqrt(25+64)) > 1.0e-10 || fabs(AC.abs() - sqrt(4+64)) > 1.0e-10)
			return false;
		if (!L1 || !L2 || !L3)
			return false;
		if ((*L1) != BC || (*L2) != AC || (*L3) != AB)
			return false;
		Vector3D* OA = xy.intersection(AC);
		if (!OA || (*OA) != A)
			return false;
		Vector3D* OB = xy.intersection(BC);
		if (!OB || (*OB) != B)
			return false;
		Vector3D* OC = yz.intersection(AC);
		if (!OC || (*OC) != C)
			return false;
		OB = yz.intersection(AB);
		if (!OB || (*OB) != B)
			return false;
		OA = xz.intersection(AB);
		if (!OA || (*OA) != A)
			return false;
		OC = xz.intersection(BC);
		if (!OC || (*OC) != C)
			return false;
	}
	if (true)
	{
		double tx = pi / 25, ty = 45 / 180 * pi, tz = pi / 13.5;
		Vector3D O1, A1(2, 0, 0), B1(0, 5, 0), C1(11.5, 12, 8), CC1(11.5,12);
		Vector3D F1(0, 0, 7), G1(2, 1.0 / 3.0, 3.5), H1(4.0, 2.0 / 3.0), I1(-1.05, 12, 7), J1(14.5,0,7), K1(2,0,7);
		Vector3D FP1, GP1(2, 1.0 / 3.0);
		Vector3D Delta(-10.5, 12, 15);
		Vector3D A = A1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D B = B1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D C = C1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		CC1 = CC1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D F = F1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D G = G1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D H = H1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D FP = FP1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D GP = GP1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D I = I1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D J = J1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D K = K1.rotate_x(tx).rotate_y(ty).rotate_z(tz) + Delta;
		Vector3D O = O1 + Delta;
		Plane P(A, B, O);
		double CP = P.distance(C);
		if (fabs(CP - 8) > 1.0e-10)
			return false;
		Vector3D CC = P.project(C);
		if (CC != CC1)
			return false;
		Vector Test1 = CC - CC1;
		if (Test1.abs() > 1.0e-10)
			return false;
		Line L(F, G);
		Vector3D* HH = P.intersection(L);
		if (!HH)
			return false;
		Vector3D Test2 = (*HH) - H;
		if (Test2.abs() > 1.0e-10)
			return false;
		Line LPP(FP, GP);
		Line* LP = P.project(L);
		if (!LP || LPP != (*LP))
			return false;
		if (P.distance(L) > 1.0e-10)
			return false;
		if (fabs(P.distance(I) - 7) > 1.0e-10)
			return false;
		Line FI(F, I);
		if (fabs(P.distance(FI) - 7) > 1.0e-10)
			return false;
		Plane P2(F, I, J);
		if (fabs(P.distance(P2) - 7) > 1.0e-10)
			return false;
		if (!P.IsOnPlane(CC1) || P.IsOnPlane(F))
			return false;
		if (!P.IsParallel(P2) || P.IsPerpendicular(P2))
			return false;
		Plane xz(O1, A1, F1);
		xz.SetThreePoints(O, A, F);
		if (P.IsParallel(xz) || !P.IsPerpendicular(xz))
			return false;
		Line L2(A, K);
		if (P.IsParallel(L2) || !P.IsPerpendicular(L2))
			return false;
		if (P.IsParallel(L) || P.IsPerpendicular(L))
			return false;
		Line IK(I, K);
		if (!P.IsParallel(IK) || P.IsPerpendicular(IK))
			return false;
		double aa = P.normal() || P.tangent();
		if (fabs(aa) > 1.0e-10)
			return false;
		Vector3D Ref = P.GetRefPoint();
		if (!P.IsOnPlane(Ref))
			return false;
	}
	NumTests += 1;
	return true;
}

