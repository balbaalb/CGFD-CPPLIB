#ifndef Test_FaceIteratorH
#define Test_FaceIteratorH
#include "../src/FaceIterator.h"
#include "../src/QuadEdge.h"
bool tester_FaceIterator(int& NumTests)
{
	Vector3D P0, P1(1), P2(0, 1), P3(0, 0, 1);
	QuadEdge T;
	T.SetAsTetrahedron(P0, P1, P2, P3);
	FaceIterator itf(&T);
	vector<Face*> f;
	Face* fp = itf.Next();
	while (fp)
	{
		f.push_back(fp);
		fp = itf.Next();
	}
	if (f.size() != 4)
		return false;
	if (itf.Next() || itf.Next())
		return false;
	if (f[0] == f[1] || f[0] == f[2] || f[0] == f[3])
		return false;
	if (f[1] == f[2] || f[1] == f[3] || f[2] == f[3])
		return false;
	//The next tests rely on the order that tetrahedron vertices, faces and edges are ordered
	Vertex* v[4];
	Face* ff[4];
	for (int i = 0; i < 4; ++i)
	{
		v[i] = T.GetVertex(i);
		ff[i] = T.GetFace(i);
		if (f[i] != ff[i])
			return false;
	}
	FaceIterator itfv0(v[0]);
	if (itfv0.Next() != ff[1] || itfv0.Next() != ff[3] || itfv0.Next() != ff[0])
		return false;
	if (itfv0.Next() || itfv0.Next())
		return false;
	FaceIterator itfv1(v[1]);
	if (itfv1.Next() != ff[2] || itfv1.Next() != ff[1] || itfv1.Next() != ff[0])
		return false;
	if (itfv1.Next() || itfv1.Next())
		return false;
	FaceIterator itfv2(v[2]);
	if (itfv2.Next() != ff[3] || itfv2.Next() != ff[2] || itfv2.Next() != ff[0])
		return false;
	if (itfv2.Next() || itfv2.Next())
		return false;
	FaceIterator itfv3(v[3]);
	if (itfv3.Next() != ff[3] || itfv3.Next() != ff[1] || itfv3.Next() != ff[2])
		return false;
	if (itfv3.Next() || itfv3.Next())
		return false;
	NumTests += 1;
	return true;
}
#endif