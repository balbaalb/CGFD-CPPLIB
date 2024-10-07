#ifndef Test_FaceH
#define Test_FaceH
#include "../src/Face.h"
#include "../src/Edge.h"
bool tester_Face(int& NumTests)
{
	EdgeBasic e1;
	FaceBasic f1;
	FaceBasic f2(f1), f3;
	f3 = f2;
	f3.SetEdge(&e1);
	if (!f3.GetEdge() || f3.GetEdge() != &e1)
		return false;
	NumTests += 1;
	return true;
}
#endif