#ifndef Test_QuadEdgeIndexH
#define Test_QuadEdgeIndexH
#include "../src/QuadEdgeIndex.h"
#include "../src/QuadEdge.h"
#include "../src/Triangulation.h"
#include "../src/VertexIterator.h"
#include "../src/EdgeIterator.h"
#include "../src/FaceIterator.h"
bool tester_QuadEdgeIndex(int& NumTests)
{
	Triangulation T = Triangulation::OffDiagonalGrid(3, 5, 3.0, 5.0);
	QuadEdge* qe = T.GetMesh2D();
	int Nv = qe->NumVertices();
	int Ne = qe->NumEdges();
	int Nf = qe->NumFaces();
	QuadEdgeIndex ind(qe);
	VertexIterator itv(qe);
	Vertex* v = itv.Next();
	while (v)
	{
		int i = ind.GetIndex(v);
		if (i >= Nv)
			return false;
		if (ind.GetVertex(i) != v)
			return false;
		v = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		int i = ind.GetIndex(e);
		if (i >= Ne)
			return false;
		if (ind.GetEdge(i) != e)
			return false;
		e = ite.Next();
	}
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		int i = ind.GetIndex(f);
		if (i >= Nf)
			return false;
		if (ind.GetFace(i) != f)
			return false;
		f = itf.Next();
	}
	for (int i = 0; i < Nv; ++i)
	{
		Vertex* v = ind.GetVertex(i);
		if (!v)
			return false;
		if (ind.GetIndex(v) != i)
			return false;
	}
	for (int i = 0; i < Ne; ++i)
	{
		Edge* e = ind.GetEdge(i);
		if (!e)
			return false;
		if (ind.GetIndex(e) != i)
			return false;
	}
	for (int i = 0; i < Nf; ++i)
	{
		Face* f = ind.GetFace(i);
		if (!f)
			return false;
		if (ind.GetIndex(f) != i)
			return false;
	}
	++NumTests;
	return true;
}

#endif