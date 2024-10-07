#ifndef Test_VertexIteratorH
#define Test_VertexIteratorH
#include "../src/VertexIterator.h"
#include "../src/Vertex.h"
#include "../src/QuadEdge.h"
bool tester_VertexIterator(int& NumTests)
{
	Vector3D P0, P1(1), P2(0, 1), P3(0, 0, 1);
	QuadEdge T;
	T.SetAsTetrahedron(P0, P1, P2, P3);
	VertexIterator itv(&T);
	vector<Vertex*> v;
	Vertex* vp = itv.Next();
	while (vp)
	{
		v.push_back(vp);
		vp = itv.Next();
	}
	if (v.size() != 4)
		return false;
	if (itv.Next() || itv.Next())
		return false;
	if (v[0] == v[1] || v[0] == v[2] || v[0] == v[3])
		return false;
	if (v[1] == v[2] || v[1] == v[3] || v[2] == v[3])
		return false;
	//Testing GetNext and GetPrev
	VertexIterator itv2(&T);
	if (itv2.Next() != v[0] || itv2.GetNext() != v[1] || itv2.GetPrev() != v[3])
		return false;
	if (itv2.Next() != v[1] || itv2.GetNext() != v[2] || itv2.GetPrev() != v[0])
		return false;
	if (itv2.Next() != v[2] || itv2.GetNext() != v[3] || itv2.GetPrev() != v[1])
		return false;
	if (itv2.Next() != v[3] || itv2.GetNext() != v[0] || itv2.GetPrev() != v[2])
		return false;
	if (itv2.Next())
		return false;
	//The next tests rely on the order that tetrahedron vertices, faces and edges are ordered
	Face* f[4];
	Vertex* vv[4];
	for (int i = 0; i < 4; ++i)
	{
		f[i] = T.GetFace(i);
		vv[i] = T.GetVertex(i);
		if (v[i] != vv[i])
			return false;
	}
	VertexIterator itvf0(f[0]);
	if (itvf0.Next() != vv[1] || itvf0.Next() != vv[0] || itvf0.Next() != vv[2])
		return false;
	if (itvf0.Next() || itvf0.Next())
		return false;
	//Testing GetNext and GetPrev
	VertexIterator itvf00(f[0]);
	if (itvf00.Next() != vv[1] || itvf00.GetNext() != vv[0] || itvf00.GetPrev() != vv[2])
		return false;
	if (itvf00.Next() != vv[0] || itvf00.GetNext() != vv[2] || itvf00.GetPrev() != vv[1])
		return false;
	if (itvf00.Next() != vv[2] || itvf00.GetNext() != vv[1] || itvf00.GetPrev() != vv[0])
		return false;
	if (itvf00.Next())
		return false;
	VertexIterator itvf1(f[1]);
	if (itvf1.Next() != vv[0] || itvf1.Next() != vv[1] || itvf1.Next() != vv[3])
		return false;
	if (itvf1.Next() || itvf1.Next())
		return false;
	//Testing GetNext and GetPrev
	VertexIterator itvf11(f[1]);
	if (itvf11.Next() != vv[0] || itvf11.GetNext() != vv[1] || itvf11.GetPrev() != vv[3])
		return false;
	if (itvf11.Next() != vv[1] || itvf11.GetNext() != vv[3] || itvf11.GetPrev() != vv[0])
		return false;
	if (itvf11.Next() != vv[3] || itvf11.GetNext() != vv[0] || itvf11.GetPrev() != vv[1])
		return false;
	if (itvf11.Next())
		return false;
	VertexIterator itvf2(f[2]);
	if (itvf2.Next() != vv[1] || itvf2.Next() != vv[2] || itvf2.Next() != vv[3])
		return false;
	if (itvf2.Next() || itvf2.Next())
		return false;
	//Testing GetNext and GetPrev
	VertexIterator itvf22(f[2]);
	if (itvf22.Next() != vv[1] || itvf22.GetNext() != vv[2] || itvf22.GetPrev() != vv[3])
		return false;
	if (itvf22.Next() != vv[2] || itvf22.GetNext() != vv[3] || itvf22.GetPrev() != vv[1])
		return false;
	if (itvf22.Next() != vv[3] || itvf22.GetNext() != vv[1] || itvf22.GetPrev() != vv[2])
		return false;
	if (itvf22.Next())
		return false;
	VertexIterator itvf3(f[3]);
	if (itvf3.Next() != vv[2] || itvf3.Next() != vv[0] || itvf3.Next() != vv[3])
		return false;
	if (itvf3.Next() || itvf3.Next())
		return false;
	//Testing GetNext and GetPrev
	VertexIterator itvf33(f[3]);
	if (itvf33.Next() != vv[2] || itvf33.GetNext() != vv[0] || itvf33.GetPrev() != vv[3])
		return false;
	if (itvf33.Next() != vv[0] || itvf33.GetNext() != vv[3] || itvf33.GetPrev() != vv[2])
		return false;
	if (itvf33.Next() != vv[3] || itvf33.GetNext() != vv[2] || itvf33.GetPrev() != vv[0])
		return false;
	if (itvf33.Next())
		return false;
	NumTests += 1;
	return true;
}
#endif