#include "MathUtils.h"
#include "LineSegment2D.h"
void LineSegment2D::CopyBody(const LineSegment2D& rhs)
{
	for (int i = 0; i < 2; ++i)
	{
		this->points[i] = rhs.points[i];
	}
}
LineSegment2D::LineSegment2D()
{

}
LineSegment2D::LineSegment2D(const Vector3D& A, const Vector3D& B)
{
	this->points[0] = A;
	this->points[1] = B;
}
LineSegment2D::LineSegment2D(const LineSegment2D& rhs)
{
	this->CopyBody(rhs);
}
void LineSegment2D::operator=(const LineSegment2D& rhs)
{
	this->CopyBody(rhs);
}
double LineSegment2D::GetArea() const
{
	return 0;
}
double LineSegment2D::GetPerimeter() const
{
	return this->points[0].distance(this->points[1]);
}
Vector3D LineSegment2D::GetPoint(int i) const
{
	return this->points[i];
}
void LineSegment2D::SetPoint(int i, const Vector3D& p)
{
	this->points[i] = p;
}
bool LineSegment2D::isOn(const Vector3D& P) const
{
	if (this->GetRelationTo(P) == ON_LINE_SEGMENT)
		return true;
	return false;
}
bool LineSegment2D::isIntersecting(const LineSegment2D& L) const
{
	if (this->operator==(L))
		return true;
	/*Vector3D V = this->points[1] - this->points[0];
	Vector3D L0 = L.points[0] - this->points[0];
	Vector3D L1 = L.points[1] - this->points[0];
	double c0 = V.CrossProduct2D(L0);
	double c1 = V.CrossProduct2D(L1);
	Vector3D W = L.points[1] - L.points[0];
	Vector3D S0 = this->points[0] - L.points[0];
	Vector3D S1 = this->points[1] - L.points[0];
	double c2 = W.CrossProduct2D(S0);
	double c3 = W.CrossProduct2D(S1);
	if (fabs(c0) > 1.0e-10 && fabs(c1) > 1.0e-10 && c0 * c1 < 0 && c2 * c3 < 0)
		return true;
	//Needs expansion*/
	Vector3D A = this->GetPoint(0);
	Vector3D B = this->GetPoint(1);
	Vector3D C = L.GetPoint(0);
	Vector3D D = L.GetPoint(1);
	if (this->GetRelationTo(C) == ON_LINE_SEGMENT && this->GetRelationTo(D) == ON_LINE_SEGMENT)
		return true;
	if (L.GetRelationTo(A) == ON_LINE_SEGMENT && L.GetRelationTo(B) == ON_LINE_SEGMENT)
		return true;
	if (this->GetRelationTo(C) == ON_LINE_SEGMENT && C != A && C != B)
		return true;
	if (this->GetRelationTo(D) == ON_LINE_SEGMENT && D != A && D != B)
		return true;
	if (L.GetRelationTo(A) == ON_LINE_SEGMENT && A != C && A != D)
		return true;
	if (L.GetRelationTo(B) == ON_LINE_SEGMENT && B != C && B != D)
		return true;
	if (this->GetRelationTo(C) != ON_LINE_SEGMENT && this->GetRelationTo(D) != ON_LINE_SEGMENT
		&& L.GetRelationTo(A) != ON_LINE_SEGMENT && L.GetRelationTo(B) != ON_LINE_SEGMENT
		&& this->GetRelationTo(C) != this->GetRelationTo(D) 
		&& L.GetRelationTo(A) != L.GetRelationTo(B))
		return true;
	//Error: 1 on line segment, another on extension but covers the line Example: A(0,0), B(1,0), C = A & D = (2,0)
	return false;
}
bool LineSegment2D::isIntersectingHorizontalRayFrom(const Vector3D& P)
{
	return this->isIntersectingRayFrom(P,0);
}
bool LineSegment2D::isIntersectingRayFrom(const Vector3D& P, double theta)
{//Does intersect the line from P(x,y)
	if (theta > 2.0 * pi)
		theta -= floor(theta / 2.0 / pi) * 2.0 * pi;
	if (theta == 2.0 * pi)
		theta = 0;
	Vector3D PA = this->points[0] - P;
	Vector3D PB = this->points[1] - P;
	double maxDist = std::fmax(PA.abs(), PB.abs());
	double farDistance = maxDist * 2.0;
	Vector3D Q;
	if (theta == pi / 2.0)
	{
		Q(0) = P(0);
		Q(1) = P(1) + farDistance;
	}
	else if (theta == 30 * pi / 2.0)
	{
		Q(0) = P(0);
		Q(1) = P(1) - farDistance;
	}
	else
	{
		double m = tan(theta);
		double b = P(1) - m * P(0);
		Q(0) = P(0) + maxDist * cos(theta);
		Q(1) = P(1) + maxDist * sin(theta);
	}
	LineSegment2D L2(P, Q);
	return this->isIntersecting(L2);
}
bool LineSegment2D::operator==(const LineSegment2D& rhs) const
{
	if (this->points[0] == rhs.points[0] && this->points[1] == rhs.points[1])
		return true;
	if (this->points[0] == rhs.points[1] && this->points[1] == rhs.points[0])
		return true;
	return false;
}
bool LineSegment2D::operator!=(const LineSegment2D& rhs) const
{
	return !(this->operator==(rhs));
}
POINT_LINESEGMENT2D_RELATION LineSegment2D::GetRelationTo(const Vector3D& P) const
{
	Vector3D L = this->points[1] - this->points[0];
	Vector3D LP = P - this->points[0];
	double LxP = L.CrossProduct2D(LP);
	if (LxP > 1.0e-10)
		return ON_LEFT;
	if (LxP < -1.0e-10)
		return ON_RIGHT;
	if (P == this->points[0] || P == this->points[1])
		return ON_LINE_SEGMENT;
	if ((L || LP) > -1.0e-10 && LP.abs() <= L.abs())
		return ON_LINE_SEGMENT;
	return ON_EXTENSION;
}
bool tester_LineSegment2D(int& NumTests)
{
	Vector3D A(1, 1), B(3, 2), C(3, 1.5), D(3, 0), E(1, 2), F(2.5,2.5),G(2,2),H(1.5,1.5);
	LineSegment2D L0(A,B);
	LineSegment2D L1 = L0;
	if (L0 != L1)
		return false;
	L0.SetPoint(1, C);
	LineSegment2D L2 = L0;
	if (L0 == L1 || L0 != L2)
		return false;
	L0.SetPoint(0, D);
	LineSegment2D L3 = L0;
	if (L0 == L1 || L0 == L2 || L0 != L3)
		return false;
	L0.SetPoint(0, E);
	LineSegment2D L4 = L0;
	if (L0 == L1 || L0 == L2 || L0 == L3 || L0 != L4)
		return false;
	if (L1.GetPoint(0) != A || L1.GetPoint(1) != B)
		return false;
	if (L2.GetPoint(0) != A || L2.GetPoint(1) != C)
		return false;
	if (L3.GetPoint(0) != D || L3.GetPoint(1) != C)
		return false;
	if (L4.GetPoint(0) != E || L4.GetPoint(1) != C)
		return false;
	if (L1.isIntersecting(L2))
		return false;
	if (!L1.isIntersecting(L1))
		return false;
	if (L1.isIntersecting(L3))
		return false;
	if (!L1.isIntersecting(L4))
		return false;
	if (fabs(L3.GetPerimeter() - 1.5) > 1.0e-10)
		return false;
	if (L3.GetArea())
		return false;
	if (fabs(L4.GetPerimeter() - sqrt(4 + 0.25)) > 1.0e-10)
		return false;
	LineSegment2D AF(A, F), GH(G, H), HG(H, G);
	if (!AF.isIntersecting(GH) || !AF.isIntersecting(HG) || !GH.isIntersecting(AF) || !HG.isIntersecting(AF))
		return false;
	if (!tester_LineSegment2D_1(NumTests))
		return false;
	if (!tester_LineSegment2D_2(NumTests))
		return false;
	if (!tester_LineSegment2D_3(NumTests))
		return false;
	NumTests += 1;
	return true;
}
bool tester_LineSegment2D_1(int& NumTests)
{
	Vector3D A, B(1, 0), C(0.5, 1), D(0.5, 0.5), E(0, 0.5), F(2.5,2.5),G(2,2),H(1.5,1.5);
	Vector3D Delta(1.5, -5.8);
	A = A + Delta;
	B = B + Delta;
	C = C + Delta;
	D = D + Delta;
	E = E + Delta;
	F = F + Delta;
	G = G + Delta;
	H = H + Delta;
	double theta = 35.75;
	A = A.rotate_z(theta);
	B = B.rotate_z(theta);
	C = C.rotate_z(theta);
	D = D.rotate_z(theta);
	E = E.rotate_z(theta);
	F = F.rotate_z(theta);
	G = G.rotate_z(theta);
	H = H.rotate_z(theta);
	LineSegment2D DA(D, A), BC(B, C), BE(B,E), AF(A,F), GH(G,H), HG(H,G);
	if (DA.isIntersecting(BC) || BC.isIntersecting(DA))
		return false;
	if (!DA.isIntersecting(BE) || !BE.isIntersecting(DA))
		return false;
	if (DA.GetRelationTo(B) != ON_LEFT || DA.GetRelationTo(E) != ON_RIGHT || DA.GetRelationTo(F) != ON_EXTENSION)
		return false;
	if (AF.GetRelationTo(D) != ON_LINE_SEGMENT || AF.GetRelationTo(A) != ON_LINE_SEGMENT || AF.GetRelationTo(F) != ON_LINE_SEGMENT)
		return false;
	if (!AF.isIntersecting(GH) || !AF.isIntersecting(HG) || !GH.isIntersecting(AF) || !HG.isIntersecting(AF))
		return false;
	NumTests += 1;
	return true;
}
bool tester_LineSegment2D_2(int& NumTests)
{
	double h = sqrt(3.0) / 2.0;
	Vector3D A(6, 5), C(0, 1), E(6, 4), F(6 - h, 4.5);
	LineSegment2D AC(A, C), EF(E, F);
	if (!AC.isIntersecting(EF))
		return false;
	if (!EF.isIntersecting(AC))
		return false;
	NumTests += 1;
	return true;
}
bool tester_LineSegment2D_3(int& NumTests)
{
	Vector3D A, B(0, 2), C(1, 2);
	LineSegment2D AB(A, C), BC(B, C), AC(A,C);
	Vector3D* P = new Vector3D(-1, 1);
	if (!AB.isIntersectingHorizontalRayFrom(*P)) 
		return false;
	delete P; P = new Vector3D(+1, 1);
	if (AB.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(-1, 0);
	if (!AB.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(0, 2);
	if (!AB.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(0, 1);
	if (!AB.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(-1, 2);
	if (!BC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P;P = new Vector3D(0.5, 2);
	if (!BC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(0, 0);
	if (BC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(2, 2);
	if (BC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(-1, -1);
	if (AC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(-1, 3);
	if (AC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(-1, 1);
	if (!AC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(0.44445, 1);
	if (!AC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P; P = new Vector3D(0.50001, 1);
	if (AC.isIntersectingHorizontalRayFrom(*P))
		return false;
	delete P;
	NumTests += 1;
	return true;
}
