#ifndef Test_PlaneH
#define Test_PlaneH
#include "../src/Plane.h"
#include "../src/MathUtils.h"
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
#endif