#ifndef Test_LineSegment2DH
#define Test_LineSegment2DH
#include "../src/LineSegment2D.h"
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
#endif