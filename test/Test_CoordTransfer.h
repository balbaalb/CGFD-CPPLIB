#ifndef Test_CoordTransferH
#define Test_CoordTransferH
#include "../src/CoordTransfer.h"
bool tester_CoordTransfer(int& NumTests)
{
	bMatrix M(3, 3);
	M(0, 0) = 11.2; M(0, 1) = 6.8;
	M(1, 0) = -14.5; M(1, 1) = -5.9;
	M(2, 2) = 1.0;
	Vector3D C(12.5, -9.8);
	CoordTransfer CTr;
	CTr.SetRotationMatrix(M);
	CTr.Translate(C);
	CoordTransfer CTr1(CTr), CTr2;
	CTr2 = CTr1;
	Vector3D p0(12.5, -101);
	Vector3D p1 = M * p0 + C;
	Vector3D q1 = CTr2.OriginalToNew(p0);
	if (p1 != q1)
		return false;
	Vector3D q0 = CTr2.NewToOriginal(q1);
	if (q0 != p0)
		return false;
	NumTests += 1;
	return true;
}
#endif