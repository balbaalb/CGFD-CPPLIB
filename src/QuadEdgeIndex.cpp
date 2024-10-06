#include "QuadEdgeIndex.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "EdgeIterator.h"
#include "QuadEdge.h"
#include "Triangulation.h"//for testing
QuadEdgeIndex::QuadEdgeIndex(QuadEdge* qe_) : qe(qe_)
{
	int Nv = this->qe->NumVertices();
	int Ne = this->qe->NumEdges();
	int Nf = this->qe->NumFaces();
	this->VertexMap.rehash(Nv);
	this->EdgeMap.rehash(Ne);
	this->FaceMap.rehash(Nf);
	this->VertexArr.resize(Nv, 0);
	this->EdgeArr.resize(Ne, 0);
	this->FaceArr.resize(Nf, 0);
	VertexIterator itv(this->qe);
	Vertex* v = itv.Next();
	int iv = 0;
	while (v)
	{
		this->VertexArr[iv] = v;
		this->VertexMap[v] = iv;
		++iv;
		v = itv.Next();
	}
	EdgeIterator ite(this->qe);
	Edge* e = ite.Next();
	int ie = 0;
	while (e)
	{
		this->EdgeArr[ie] = e;
		this->EdgeMap[e] = ie;
		++ie;
		e = ite.Next();
	}
	FaceIterator itf(this->qe);
	Face* f = itf.Next();
	int ifc = 0;
	while (f)
	{
		this->FaceArr[ifc] = f;
		this->FaceMap[f] = ifc;
		++ifc;
		f = itf.Next();
	}
}
Vertex* QuadEdgeIndex::GetVertex(int vi)
{
	return this->VertexArr[vi];
}
Face* QuadEdgeIndex::GetFace(int fi)
{
	return this->FaceArr[fi];
}
Edge* QuadEdgeIndex::GetEdge(int ei)
{
	return this->EdgeArr[ei];
}
int QuadEdgeIndex::GetIndex(Vertex* v)
{
	return this->VertexMap[v];
}
int QuadEdgeIndex::GetIndex(Face* f)
{
	return this->FaceMap[f];
}
int QuadEdgeIndex::GetIndex(Edge* e)
{
	return this->EdgeMap[e];
}
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
