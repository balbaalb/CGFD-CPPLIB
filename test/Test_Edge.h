#ifndef Test_EdgeH
#define Test_EdgeH
#include "../src/Edge.h"
#include "../src/Vertex.h"
#include "../src/Face.h"
bool tester_Edge(int& NumTests)
{
	//A C3 graph
	VertexBasic v1, v2;
	v1.SetX(11.5); v1.SetY(-21.8); v1.SetZ(14.36);
	v2.SetX(v1.GetX() + 3.0); v2.SetY(v1.GetY() + 4.0); v2.SetZ(v1.GetZ());
	EdgeBasic e0, e1, e2, e3, e4, e5, e6, e7, e8;
	FaceBasic f1, f2;
	EdgeBasic e00(e0), e000;
	e000 = e00;
	e000.SetOrig(&v1);
	e000.SetDest(&v2);
	if (fabs(e000.GetLength() - 5.0) > 1.0e-10)
		return false;
	e000.SetRight(&f1);
	e000.SetLeft(&f2);
	e000.SetNextOrig(&e1);
	e000.SetPrevOrig(&e2);
	e000.SetNextDest(&e3);
	e000.SetPrevDest(&e4);
	e000.SetNextRight(&e5);
	e000.SetPrevRight(&e6);
	e000.SetNextLeft(&e7);
	e000.SetPrevLeft(&e8);
	if (!e000.GetOrig() || e000.GetOrig() != &v1)
		return false;
	if (!e000.GetDest() || e000.GetDest() != &v2)
		return false;
	if (!e000.GetRight() || e000.GetRight() != &f1)
		return false;
	if (!e000.GetLeft() || e000.GetLeft() != &f2)
		return false;
	if (!e000.GetNextOrig() || e000.GetNextOrig() != &e1)
		return false;
	if (!e000.GetPrevOrig() || e000.GetPrevOrig() != &e2)
		return false;
	if (!e000.GetNextDest() || e000.GetNextDest() != &e3)
		return false;
	if (!e000.GetPrevDest() || e000.GetPrevDest() != &e4)
		return false;
	if (!e000.GetNextRight() || e000.GetNextRight() != &e5)
		return false;
	if (!e000.GetPrevRight() || e000.GetPrevRight() != &e6)
		return false;
	if (!e000.GetNextLeft() || e000.GetNextLeft() != &e7)
		return false;
	if (!e000.GetPrevLeft() || e000.GetPrevLeft() != &e8)
		return false;
	e000.Reverse();
	if (!e000.GetOrig() || e000.GetOrig() != &v2)
		return false;
	if (!e000.GetDest() || e000.GetDest() != &v1)
		return false;
	if (!e000.GetRight() || e000.GetRight() != &f2)
		return false;
	if (!e000.GetLeft() || e000.GetLeft() != &f1)
		return false;
	if (!e000.GetNextOrig() || e000.GetNextOrig() != &e3)
		return false;
	if (!e000.GetPrevOrig() || e000.GetPrevOrig() != &e4)
		return false;
	if (!e000.GetNextDest() || e000.GetNextDest() != &e1)
		return false;
	if (!e000.GetPrevDest() || e000.GetPrevDest() != &e2)
		return false;
	if (!e000.GetNextRight() || e000.GetNextRight() != &e7)
		return false;
	if (!e000.GetPrevRight() || e000.GetPrevRight() != &e8)
		return false;
	if (!e000.GetNextLeft() || e000.GetNextLeft() != &e5)
		return false;
	if (!e000.GetPrevLeft() || e000.GetPrevLeft() != &e6)
		return false;
	NumTests += 1;
	return true;
}
#endif