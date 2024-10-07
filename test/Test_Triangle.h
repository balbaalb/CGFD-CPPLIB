#ifndef Test_TriangleH
#define Test_TriangleH
#include "../src/Triangle.h"
#include "../src/MathUtils.h"
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
#endif