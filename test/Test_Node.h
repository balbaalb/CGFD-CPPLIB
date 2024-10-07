#ifndef Test_NodeH
#define Test_NodeH
#include "../src/Node.h"
#include "../src/Vertex.h"
bool tester_Node(int& NumTests)
{
	Vector3D P(1.5, 3.5);
	VertexBasic* v = new VertexBasic(P);
	Node n0(v);
	Node n1(n0), n2(0);
	n2 = n1;
	if (n2.GetPoint() != P)
		return false;
	NodeT n3(v);
	n3.T = 10.5;
	NodeT n4(n3), n5(0);
	n5 = n4;
	if (n5.T != 10.5)
		return false;
	//---------------
	NodeTV n6(v);
	n6.T = 10.5;
	n6.vx = -1.5;
	n6.vy = 8.13;
	NodeTV n7(n6), n8(0);
	n8 = n7;
	if (n8.T != 10.5)
		return false;
	if (n8.vx != -1.5)
		return false;
	if (n8.vy != 8.13)
		return false;
	++NumTests;
	return true;
}
#endif