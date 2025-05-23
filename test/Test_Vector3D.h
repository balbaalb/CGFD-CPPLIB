#ifndef Test_Vector3DH
#define Test_Vector3DH
#include "../src/Vector3D.h"
#include "../src/MathUtils.h"
bool tester_Vector3D(int& NumTests)
{
	Vector3D P(1);
	P(2) = 3;
	Vector3D PP(P), PPP;
	PPP = PP;
	if (PPP.GetDim1() != 3 || PPP.GetDim2() != 1)
		return false;
	if (PPP(0) != 1 || PPP(1) != 0 || PPP(2) != 3)
		return false;
	bool isEqual = (PPP == P);
	if (!isEqual)
		return false;
	bool notEqual = (PPP != P);
	if (notEqual)
		return false;
	Vector3D Q(1.5, 2, -5);
	Vector3D PQ = Q - P;
	if (PQ.GetDim1() != 3 || PQ.GetDim2() != 1)
		return false;
	if (fabs(PQ(0) - 0.5) > 1.e15)
		return false;
	if (fabs(PQ(1) - 2) > 1.e15)
		return false;
	if (fabs(PQ(2) + 8) > 1.e15)
		return false;
	Vector3D QQ = PQ + P;
	if (QQ.GetDim1() != 3 || QQ.GetDim2() != 1)
		return false;
	if (fabs(QQ(0) - 1.5) > 1.e15)
		return false;
	if (fabs(QQ(1) - 2) > 1.e15)
		return false;
	if (fabs(QQ(2) + 5) > 1.e15)
		return false;
	double a = (P * 2.5) || (Q / 1.5);
	if (fabs(a - (1.0 * 2.5 - 3.0 * 2.5 * 5 / 1.5)) > 1.0e-10)
		return false;
	Vector3D V(-1, 2, 5.5);
	double Vabs = sqrt(1.0 + 4.0 + 5.5*5.5);
	if (fabs(V.abs() - Vabs) > 1.0e-10)
		return false;
	Vector3D Vhat = V.tangent();
	if (Vhat.GetDim1() != 3)
		return false;
	if (Vhat.GetDim2() != 1)
		return false;
	if (fabs(Vhat(0) + 1.0 / Vabs) > 1.0e-10)
		return false;
	if (fabs(Vhat(1) - 2.0 / Vabs) > 1.0e-10)
		return false;
	if (fabs(Vhat(2) - 5.5 / Vabs) > 1.0e-10)
		return false;
	bMatrix C(3, 3);
	C(0, 0) = 11.0; C(0, 1) = -5.5; C(0, 2) = -1.56;
	C(1, 0) = 10.34; C(1, 1) = 108.85; C(1, 2) = 14.5;
	C(2, 0) = -11.75; C(2, 1) = 6.85; C(2, 2) = 5.98;
	Vector3D CV = C * V;
	if (CV.GetDim1() != 3)
		return false;
	if (CV.GetDim2() != 1)
		return false;
	if (fabs(CV(0) - (11.0 * (-1) - 5.5 * 2.0 - 1.56 * 5.5)) > 1.0e-10)
		return false;
	if (fabs(CV(1) - (10.34 * (-1) + 108.85 * 2.0 + 14.5 * 5.5)) > 1.0e-10)
		return false;
	if (fabs(CV(2) - (-11.75 * (-1) + 6.85 * 2.0 + 5.98 * 5.5)) > 1.0e-10)
		return false;
	Vector3D VV = V * 1.75;
	Vector3D A = V && VV;
	if (A.GetDim1() != 3)
		return false;
	if (A.GetDim2() != 1)
		return false;
	if (A.abs() > 1.0e-10)
		return false;
	//P(0) = 1;  P(1) = -4;  P(2) = 3;
	P.SetCoords(1.0, -4.0, 3.0);
	A = P && V;
	if (fabs(A(2) - P.CrossProduct2D(V)) > 1.0e-10)
		return false;
	if (A.GetDim1() != 3)
		return false;
	if (A.GetDim2() != 1)
		return false;
	if (fabs(A(0) + 28) > 1.0e-10)
		return false;
	if (fabs(A(1) + 8.5) > 1.0e-10)
		return false;
	if (fabs(A(2) + 2) > 1.0e-10)
		return false;
	if (true)
	{
		Vector3D A(2, 1, 3);
		Vector3D B = A.rotate_z(pi / 2.0);
		if (fabs(B(0) + 1) > 1.0e-10 || fabs(B(1) - 2) > 1.0e-10 || B(2) != 3)
			return false;
		Vector3D C = B.rotate_x(pi / 2.0);
		if (fabs(C(0) + 1) > 1.0e-10 || fabs(C(1) + 3) > 1.0e-10 || fabs(C(2) - 2) > 1.0e-10)
			return false;
		Vector3D D = C.rotate_y(pi / 2);
		if (fabs(D(0) + 2) > 1.0e-10 || fabs(D(1) + 3) > 1.0e-10 || fabs(D(2) + 1) > 1.0e-10)
			return false;
	}
	if (true)
	{
		Vector3D A(-1.85, 0);
		Vector3D B;
		double theta = 65.3 / 180 * pi;
		double r = 4.35;
		Vector3D C(r * cos(theta), r * sin(theta));
		double beta = Vector3D::RotationAngle_z(A, B, C);
		if (fabs(beta - theta) > 1.0e-10)
			return false;
		theta += pi / 2.0;
		C(0) = r * cos(theta); C(1) = r * sin(theta);
		beta = Vector3D::RotationAngle_z(A, B, C);
		if (fabs(beta - theta) > 1.0e-10)
			return false;
		theta = -65.3 / 180 * pi;
		C(0) = r * cos(theta); C(1) = r * sin(theta);
		beta = Vector3D::RotationAngle_z(A, B, C);
		if (fabs(beta - theta) > 1.0e-10)
			return false;
		theta -= pi / 2.0;
		C(0) = r * cos(theta); C(1) = r * sin(theta);
		beta = Vector3D::RotationAngle_z(A, B, C);
		if (fabs(beta - theta) > 1.0e-10)
			return false;
	}
	if (true)
	{
		Vector3D A(cos(pi / 2.5), sin(pi / 2.5)), B, C(2.75, 0);
		double theta1 = Vector3D::GetAngle(A, B, C);
		if (fabs(theta1 - pi / 2.5) > 1.e-10)
			return false;
		Vector3D AA(-cos(pi / 2.5), sin(pi / 2.5));
		double theta2 = Vector3D::GetAngle(AA, B, C);
		if (fabs(theta2 - 3.0 * pi / 5.0) > 1.e-10)
			return false;
	}
	if (true)
	{
		Vector3D A(5, 0), B(1.5 * cos(pi / 3.0), 1.5 * sin(pi / 3.0)), C(0, 6.5), D(-1.85 * cos(pi / 6.0), 1.85 * sin(pi / 6.0)), E(-6.1, 0), F(-1.85 * cos(pi / 4.0), -1.85 * sin(pi / 4.0));
		Vector3D G(0, -1.6), H(1.85 * cos(pi / 5.0), -1.85 * sin(pi / 5.0));
		if (fabs(A.GetAngleAboutXAxis()) > 1.0e-10 || fabs(A.GetOrientationAboutXAxis()) > 1.0e-10)
			return false;
		if (fabs(B.GetAngleAboutXAxis() - pi / 3.0) > 1.0e-10 || fabs(B.GetOrientationAboutXAxis() - pi / 3.0) > 1.0e-10)
			return false;
		if (fabs(C.GetAngleAboutXAxis() - pi / 2.0) > 1.0e-10 || fabs(C.GetOrientationAboutXAxis() - pi / 2.0) > 1.0e-10)
			return false;
		if (fabs(D.GetAngleAboutXAxis() - 5.0 * pi / 6.0) > 1.0e-10 || fabs(D.GetOrientationAboutXAxis() - 5.0 * pi / 6.0) > 1.0e-10)
			return false;
		if (fabs(E.GetAngleAboutXAxis() - pi) > 1.0e-10 || fabs(E.GetOrientationAboutXAxis()) > 1.0e-10)
			return false;//Here
		if (fabs(F.GetAngleAboutXAxis() - 5.0 * pi / 4.0) > 1.0e-10 || fabs(F.GetOrientationAboutXAxis() - pi / 4.0) > 1.0e-10)
			return false;
		if (fabs(G.GetAngleAboutXAxis() - 3.0 * pi / 2.0) > 1.0e-10 || fabs(G.GetOrientationAboutXAxis() - pi / 2.0) > 1.0e-10)
			return false;
		if (fabs(H.GetAngleAboutXAxis() - 9.0 * pi / 5.0) > 1.0e-10 || fabs(H.GetOrientationAboutXAxis() - 4.0 * pi / 5.0) > 1.0e-10)
			return false;
		Vector3D O(0.5, 0.5), P00, P10(1), P11(1, 1), P01(0, 1);
		Vector3D OP00 = P00 - O;
		Vector3D OP10 = P10 - O;
		Vector3D OP11 = P11 - O;
		Vector3D OP01 = P01 - O;
		if (fabs(OP00.GetAngleAboutXAxis() - 1.25 * pi) > 1.0e-10 || fabs(OP00.GetOrientationAboutXAxis() - pi / 4.0) > 1.0e-10)
			return false;
		if (fabs(OP10.GetAngleAboutXAxis() - 1.75 * pi) > 1.0e-10 || fabs(OP10.GetOrientationAboutXAxis() - 3.0 * pi / 4.0) > 1.0e-10)
			return false;
		if (fabs(OP11.GetAngleAboutXAxis() - pi / 4.0) > 1.0e-10 || fabs(OP11.GetOrientationAboutXAxis() - pi / 4.0) > 1.0e-10)
			return false;
		if (fabs(OP01.GetAngleAboutXAxis() - 3.0 * pi / 4.0) > 1.0e-10 || fabs(OP01.GetOrientationAboutXAxis() - 3.0 * pi / 4.0) > 1.0e-10)
			return false;
	}
	NumTests += 1;
	return true;
}
#endif