#ifndef Test_TriangulationH
#define Test_TriangulationH
#include "../src/Triangulation.h"
#include "../src/QuadEdge.h"
#include "../src/VertexIterator.h"
#include "../src/EdgeIterator.h"
#include "../src/FaceIterator.h"
#include "../src/DelaunayLifting.h"
#include "../src/Vertex.h"
#include "../src/Edge.h"
#include "../src/Face.h"
bool tester_Triangulation_1_0(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "0"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false;
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 17 || T.GetMesh2D()->NumFaces() != 11)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 4)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 3)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 4)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[10] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 6)
		return false;
	if (f[0]->GetEdge() != e[2] || f[0]->GetDegree() != 3)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 4)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[6] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[7] || f[6]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[1] || f[7]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[9] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[3])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[13] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[14] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[13] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[0])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[3] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[3])
		return false;
	if (e[3]->GetOrig() != v[0] || e[3]->GetDest() != v[3] || e[3]->GetLeft() != f[1] || e[3]->GetRight() != f[0])
		return false;
	if (e[3]->GetNextOrig() != e[2] || e[3]->GetNextDest() != e[4] || e[3]->GetNextLeft() != e[5] || e[3]->GetNextRight() != e[2])
		return false;
	if (e[3]->GetPrevOrig() != e[6] || e[3]->GetPrevDest() != e[15] || e[3]->GetPrevLeft() != e[0] || e[3]->GetPrevRight() != e[4])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[0])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[3])
		return false;
	if (e[4]->GetPrevOrig() != e[3] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[14] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[3])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[4] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[3] || e[6]->GetNextDest() != e[12] || e[6]->GetNextLeft() != e[10] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[9] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[6] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[13] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[10] || e[7]->GetPrevDest() != e[13] || e[7]->GetPrevLeft() != e[12] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[0] || e[9]->GetDest() != v[6] || e[9]->GetLeft() != f[3] || e[9]->GetRight() != f[4])
		return false;
	if (e[9]->GetNextOrig() != e[6] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[8] || e[9]->GetNextRight() != e[6])
		return false;
	if (e[9]->GetPrevOrig() != e[0] || e[9]->GetPrevDest() != e[8] || e[9]->GetPrevLeft() != e[2] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[6] || e[10]->GetDest() != v[5] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[7] || e[10]->GetNextLeft() != e[12] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[12] || e[10]->GetPrevLeft() != e[11] || e[10]->GetPrevRight() != e[6])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[10] || e[11]->GetRight() != f[5])
		return false;
	if (e[11]->GetNextOrig() != e[8] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[16] || e[11]->GetNextRight() != e[10])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[16] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[12])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[10] || e[12]->GetNextDest() != e[13] || e[12]->GetNextLeft() != e[11] || e[12]->GetNextRight() != e[7])
		return false;
	if (e[12]->GetPrevOrig() != e[6] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[10] || e[12]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[1] || e[13]->GetDest() != v[7] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[7] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[12] || e[13]->GetNextRight() != e[1])
		return false;
	if (e[13]->GetPrevOrig() != e[1] || e[13]->GetPrevDest() != e[12] || e[13]->GetPrevLeft() != e[7] || e[13]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[7] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[13] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[5] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[1] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[3] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[11] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	if (T.GetOppositeFace(f[5], e[12]) != f[6])
		return false;
	if (T.GetOppositeFace(f[6], e[12]) != f[5])
		return false;
	if (T.GetOppositeFace(f[1], e[5]) != f[8])
		return false;
	if (T.GetOppositeFace(f[8], e[5]) != f[1])
		return false;
	if (T.GetOppositeFace(f[0], e[15]))
		return false;
	if (T.GetOppositeEdge(f[8], v[7]) != e[5])
		return false;
	if (T.GetOppositeEdge(f[8], v[2]) != e[15])
		return false;
	if (T.GetOppositeEdge(f[8], v[3]) != e[14])
		return false;
	if (T.GetOppositeEdge(f[8], v[4]))
		return false;
	if (T.GetRightFaceOpposite(e[10]) != v[0])
		return false;
	if (T.GetLeftFaceOpposite(e[10]) != v[7])
		return false;
	if (T.GetRightFaceOpposite(e[2]) != v[3])
		return false;
	if (T.GetLeftFaceOpposite(e[2]) != v[6])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_1(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "1_removeE3"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 4)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 4)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 3)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[10] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 6)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[6] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[7] || f[6]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[1] || f[7]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[9] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[13] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[14] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[13] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[15] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[14] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[4] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[12] || e[6]->GetNextLeft() != e[10] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[9] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[6] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[13] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[10] || e[7]->GetPrevDest() != e[13] || e[7]->GetPrevLeft() != e[12] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[0] || e[9]->GetDest() != v[6] || e[9]->GetLeft() != f[3] || e[9]->GetRight() != f[4])
		return false;
	if (e[9]->GetNextOrig() != e[6] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[8] || e[9]->GetNextRight() != e[6])
		return false;
	if (e[9]->GetPrevOrig() != e[0] || e[9]->GetPrevDest() != e[8] || e[9]->GetPrevLeft() != e[2] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[6] || e[10]->GetDest() != v[5] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[7] || e[10]->GetNextLeft() != e[12] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[12] || e[10]->GetPrevLeft() != e[11] || e[10]->GetPrevRight() != e[6])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[10] || e[11]->GetRight() != f[5])
		return false;
	if (e[11]->GetNextOrig() != e[8] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[16] || e[11]->GetNextRight() != e[10])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[16] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[12])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[10] || e[12]->GetNextDest() != e[13] || e[12]->GetNextLeft() != e[11] || e[12]->GetNextRight() != e[7])
		return false;
	if (e[12]->GetPrevOrig() != e[6] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[10] || e[12]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[1] || e[13]->GetDest() != v[7] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[7] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[12] || e[13]->GetNextRight() != e[1])
		return false;
	if (e[13]->GetPrevOrig() != e[1] || e[13]->GetPrevDest() != e[12] || e[13]->GetPrevLeft() != e[7] || e[13]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[7] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[13] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[5] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[1] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[4] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[11] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_2(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "2_flipE13"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (e[10]->index != 10 || f[4]->index != 4)
		return false;
	if (e[13]->index != 13 || f[7]->index != 7 || f[6]->index != 6)
		return false;
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 4)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 5)
		return false;
	if (v[6]->GetEdge() != e[10] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 5)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[6] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[13] || f[7]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[9] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[7] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[13] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[7] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[15] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[13] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[4] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[10] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[9] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[7] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[1] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[10] || e[7]->GetPrevDest() != e[1] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[0] || e[9]->GetDest() != v[6] || e[9]->GetLeft() != f[3] || e[9]->GetRight() != f[4])
		return false;
	if (e[9]->GetNextOrig() != e[6] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[8] || e[9]->GetNextRight() != e[6])
		return false;
	if (e[9]->GetPrevOrig() != e[0] || e[9]->GetPrevDest() != e[8] || e[9]->GetPrevLeft() != e[2] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[6] || e[10]->GetDest() != v[5] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[7] || e[10]->GetNextLeft() != e[12] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[12] || e[10]->GetPrevLeft() != e[11] || e[10]->GetPrevRight() != e[6])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[10] || e[11]->GetRight() != f[5])
		return false;
	if (e[11]->GetNextOrig() != e[8] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[16] || e[11]->GetNextRight() != e[10])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[16] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[12])
		return false;
	if (e[12]->GetOrig() != v[7] || e[12]->GetDest() != v[5] || e[12]->GetLeft() != f[6] || e[12]->GetRight() != f[5])
		return false;
	if (e[12]->GetNextOrig() != e[14] || e[12]->GetNextDest() != e[10] || e[12]->GetNextLeft() != e[13] || e[12]->GetNextRight() != e[11])
		return false;
	if (e[12]->GetPrevOrig() != e[11] || e[12]->GetPrevDest() != e[13] || e[12]->GetPrevLeft() != e[14] || e[12]->GetPrevRight() != e[10])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[6] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[12] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[13] || e[14]->GetPrevDest() != e[12] || e[14]->GetPrevLeft() != e[13] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[4] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[11] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[5] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[1])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_3(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "3_flipE10"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 3)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 6)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[13] || f[7]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[7] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[13] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[7] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[15] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[13] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[7] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[1] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[1] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[7] || e[11]->GetDest() != v[6] || e[11]->GetLeft() != f[4] || e[11]->GetRight() != f[10])
		return false;
	if (e[11]->GetNextOrig() != e[12] || e[11]->GetNextDest() != e[8] || e[11]->GetNextLeft() != e[9] || e[11]->GetNextRight() != e[16])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[9] || e[11]->GetPrevLeft() != e[10] || e[11]->GetPrevRight() != e[8])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[6] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[12] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[13] || e[14]->GetPrevDest() != e[12] || e[14]->GetPrevLeft() != e[13] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[4] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[10] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[5] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[1])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[16] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_4(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "4_connectV2-V6"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 6)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[4] || v[3]->GetDegree() != 2)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[10] || v[7]->GetDegree() != 4)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[16] || f[8]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[13] || f[7]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[15] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[16] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[7] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[13] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[7] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[15] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[5] || e[4]->GetPrevDest() != e[15] || e[4]->GetPrevLeft() != e[5] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[9] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[16] || e[5]->GetNextDest() != e[4] || e[5]->GetNextLeft() != e[4] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[15] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[7] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[1] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[1] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[16] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[16] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[8] || e[11]->GetRight() != f[4])
		return false;
	if (e[11]->GetNextOrig() != e[16] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[14] || e[11]->GetNextRight() != e[9])
		return false;
	if (e[11]->GetPrevOrig() != e[9] || e[11]->GetPrevDest() != e[10] || e[11]->GetPrevLeft() != e[16] || e[11]->GetPrevRight() != e[10])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[7] || e[14]->GetDest() != v[2] || e[14]->GetLeft() != f[8] || e[14]->GetRight() != f[6])
		return false;
	if (e[14]->GetNextOrig() != e[10] || e[14]->GetNextDest() != e[1] || e[14]->GetNextLeft() != e[16] || e[14]->GetNextRight() != e[12])
		return false;
	if (e[14]->GetPrevOrig() != e[12] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[11] || e[14]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[15] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[1])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[14] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	if (e[15]->GetOrig() != v[2] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[10] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[13] || e[15]->GetNextDest() != e[4] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[5])
		return false;
	if (e[15]->GetPrevOrig() != e[16] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[16] || e[15]->GetPrevRight() != e[4])
		return false;
	if (e[16]->GetOrig() != v[2] || e[16]->GetDest() != v[6] || e[16]->GetLeft() != f[8] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[15] || e[16]->GetNextDest() != e[8] || e[16]->GetNextLeft() != e[11] || e[16]->GetNextRight() != e[15])
		return false;
	if (e[16]->GetPrevOrig() != e[5] || e[16]->GetPrevDest() != e[11] || e[16]->GetPrevLeft() != e[14] || e[16]->GetPrevRight() != e[8])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_5(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "5_removeE1"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 15 || T.GetMesh2D()->NumFaces() != 9)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[7] || v[1]->GetDegree() != 2)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 5)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[4] || v[3]->GetDegree() != 2)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[10] || v[7]->GetDegree() != 4)
		return false;
	if (f[1]->GetEdge() != e[0] || f[1]->GetDegree() != 6)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[16] || f[8]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[15] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[16] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[7] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[7])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[15] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[5] || e[4]->GetPrevDest() != e[15] || e[4]->GetPrevLeft() != e[5] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[9] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[16] || e[5]->GetNextDest() != e[4] || e[5]->GetNextLeft() != e[4] || e[5]->GetNextRight() != e[13])
		return false;
	if (e[5]->GetPrevOrig() != e[14] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[15] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[1] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[0] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[0] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[16] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[16] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[8] || e[11]->GetRight() != f[4])
		return false;
	if (e[11]->GetNextOrig() != e[16] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[14] || e[11]->GetNextRight() != e[9])
		return false;
	if (e[11]->GetPrevOrig() != e[9] || e[11]->GetPrevDest() != e[10] || e[11]->GetPrevLeft() != e[16] || e[11]->GetPrevRight() != e[10])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[7] || e[14]->GetDest() != v[2] || e[14]->GetLeft() != f[8] || e[14]->GetRight() != f[6])
		return false;
	if (e[14]->GetNextOrig() != e[10] || e[14]->GetNextDest() != e[5] || e[14]->GetNextLeft() != e[16] || e[14]->GetNextRight() != e[12])
		return false;
	if (e[14]->GetPrevOrig() != e[12] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[11] || e[14]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[1])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[15] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[5])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[14] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	if (e[15]->GetOrig() != v[2] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[10] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[13] || e[15]->GetNextDest() != e[4] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[5])
		return false;
	if (e[15]->GetPrevOrig() != e[16] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[16] || e[15]->GetPrevRight() != e[4])
		return false;
	if (e[16]->GetOrig() != v[2] || e[16]->GetDest() != v[6] || e[16]->GetLeft() != f[8] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[15] || e[16]->GetNextDest() != e[8] || e[16]->GetNextLeft() != e[11] || e[16]->GetNextRight() != e[15])
		return false;
	if (e[16]->GetPrevOrig() != e[5] || e[16]->GetPrevDest() != e[11] || e[16]->GetPrevLeft() != e[14] || e[16]->GetPrevRight() != e[8])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_6(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "6_removeE5");
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false;
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false;
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false;
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 7 || T.GetMesh2D()->NumEdges() != 13 || T.GetMesh2D()->NumFaces() != 8)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[7] || v[1]->GetDegree() != 2)
		return false;
	if (v[2]->GetEdge() != e[16] || v[2]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[10] || v[7]->GetDegree() != 4)
		return false;
	if (f[1]->GetEdge() != e[0] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[16] || f[8]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[16] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[7] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[7])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[15])
		return false;
	if (e[2]->GetPrevOrig() != e[15] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[1] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[0] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[0] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[16] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[16] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[8] || e[11]->GetRight() != f[4])
		return false;
	if (e[11]->GetNextOrig() != e[16] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[14] || e[11]->GetNextRight() != e[9])
		return false;
	if (e[11]->GetPrevOrig() != e[9] || e[11]->GetPrevDest() != e[10] || e[11]->GetPrevLeft() != e[16] || e[11]->GetPrevRight() != e[10])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[7] || e[14]->GetDest() != v[2] || e[14]->GetLeft() != f[8] || e[14]->GetRight() != f[6])
		return false;
	if (e[14]->GetNextOrig() != e[10] || e[14]->GetNextDest() != e[16] || e[14]->GetNextLeft() != e[16] || e[14]->GetNextRight() != e[12])
		return false;
	if (e[14]->GetPrevOrig() != e[12] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[11] || e[14]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[1])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[15] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[15])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[14] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	if (e[15]->GetOrig() != v[2] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[10] || e[15]->GetRight() != f[1])
		return false;
	if (e[15]->GetNextOrig() != e[13] || e[15]->GetNextDest() != e[2] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[13])
		return false;
	if (e[15]->GetPrevOrig() != e[16] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[16] || e[15]->GetPrevRight() != e[2])
		return false;
	if (e[16]->GetOrig() != v[2] || e[16]->GetDest() != v[6] || e[16]->GetLeft() != f[8] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[15] || e[16]->GetNextDest() != e[8] || e[16]->GetNextLeft() != e[11] || e[16]->GetNextRight() != e[15])
		return false;
	if (e[16]->GetPrevOrig() != e[14] || e[16]->GetPrevDest() != e[11] || e[16]->GetPrevLeft() != e[14] || e[16]->GetPrevRight() != e[8])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1(int& NumTests)
{
	//--------------0. Create Triangulation -----------------------------------
	Vector3D A, B(2), C(2, 2), D(0, 2), E(0.7, 1.3), F(1.6, 0.4), G(0.9, 1.15), H(1.81, 1.71);
	vector<Vector3D> input = { A,B,C,D,E,F,G,H };
	Triangulation T = DelaunayLifting::Triangulate(input);
	int Nv = 8;
	int h = 4;
	int t = 2 * Nv - h - 2;
	int Nf = t + 1;
	int Ne = Nv + Nf - 2;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_0(NumTests, T,input))
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		v[vv->index] = vv;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		f[ff->index] = ff;
		ff = itf.Next();
	}
	//--------------1. Remove E3  -------------------------------------------
	T.ShrinkBoundary(e[3]);
	--Ne;
	--Nf;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_1(NumTests, T,input))
		return false;
	//--------------2. Flip E13  ------------------------------------------
	T.CheckThenFlip(e[13]);//causes e[12] to reverse
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_2(NumTests, T, input))
		return false;
	//--------------3. Flip E10  ---------------------------------------------
	T.CheckThenFlip(e[10]);//causes e[9], e[11] and e[12] (again) to reverse
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_3(NumTests, T, input))
		return false;
	//--------------4. Connect v2-v6  -------------------------------------
	T.Connect(v[2], v[6]);
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_4(NumTests, T, input))
		return false;
	//--------------5. Remove E1  ------------------------------------------
	T.ShrinkBoundary(e[1]);
	--Ne;
	--Nf;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_5(NumTests, T, input))
		return false;
	//--------------6. Remove E5  ------------------------------------------
	T.ShrinkBoundary(e[5]);
	if (!tester_Triangulation_1_6(NumTests, T, input))
		return false;
	--Nv;
	Ne -= 2;
	--Nf;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	//--------------------------------------------------------------------------
	//T.PrintTriangulation();
	NumTests += 1;
	return true;
}
bool tester_Triangulation_2(int& NumTests)
{
	int Nx = 1;
	int Ny = 1;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	//T.PrintTriangulation();
	T.Draw("..\\Data\\tester_Triangulation_2.bmp");
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() !=  Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!T.TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_3(int& NumTests)
{
	int Nx = 2;
	int Ny = 2;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	T.Draw("..\\Data\\tester_Triangulation_3.bmp");
	//T.PrintTriangulation();
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!T.TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_4(int& NumTests)
{
	int Nx = 3;
	int Ny = 4;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	//T.PrintTriangulation();
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!T.TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_5(int& NumTests)
{
	int Nx = 5;
	int Ny = 7;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	//T.PrintTriangulation();
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	/*if (!T.TestIntegrity())
		return false;*///Passes okay but takes 5 sec to do it.
	NumTests += 1;
	return true;
}
bool tester_Triangulation_6(int& NumTests)
{
	Vector3D A, B(2), C(2, 2), D(0, 2), E(0.55, 0.54), F(1.55, 0.45), G(1.56, 1.55), H(0.54, 1.54);
	vector<Vector3D> input = { A,B,C,D,E,F,G,H };
	Triangulation T = DelaunayLifting::Triangulate(input);
	T.PrintTriangulation("..\\Data", "0");
	T.Draw("..\\Data\\tester_Triangulation_6.bmp");
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		v[vv->index] = vv;
		vv = itv.Next();
	}
	vector<Vertex*> boundary = { v[0] , v[1], v[5], v[7], v[4] };
	T.ImposeBoundary(boundary);
	//T.PrintTriangulation();
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv2(qe);
	Vertex* vv2 = itv2.Next();
	while (vv2)
	{
		if (v[vv2->index])
			return false;
		v[vv2->index] = vv2;
		if (vv2->GetPoint() != input[vv2->index])
			return false;
		vv2 = itv2.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false;
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false;
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 5 || T.GetMesh2D()->NumEdges() != 7 || T.GetMesh2D()->NumFaces() != 4)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 3)
		return false;
	if (v[1]->GetEdge() != e[7] || v[1]->GetDegree() != 2)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[8] || v[5]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[11] || v[7]->GetDegree() != 2)
		return false;
	if (f[1]->GetEdge() != e[0] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[11] || f[5]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[6] || e[0]->GetNextDest() != e[7] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[7])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[6] || e[2]->GetNextRight() != e[15])
		return false;
	if (e[2]->GetPrevOrig() != e[15] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[3] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[7] || e[6]->GetNextLeft() != e[8] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[0] || e[6]->GetPrevDest() != e[11] || e[6]->GetPrevLeft() != e[2] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[1] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[8] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[0] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[6] || e[7]->GetPrevDest() != e[0] || e[7]->GetPrevLeft() != e[11] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[5] || e[8]->GetLeft() != f[5] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[11] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[7] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[6])
		return false;
	if (e[15]->GetOrig() != v[7] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[5] || e[15]->GetRight() != f[1])
		return false;
	if (e[15]->GetNextOrig() != e[11] || e[15]->GetNextDest() != e[2] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[11])
		return false;
	if (e[15]->GetPrevOrig() != e[11] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[11] || e[15]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[5] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[5] || e[11]->GetRight() != f[1])
		return false;
	if (e[11]->GetNextOrig() != e[6] || e[11]->GetNextDest() != e[15] || e[11]->GetNextLeft() != e[15] || e[11]->GetNextRight() != e[7])
		return false;
	if (e[11]->GetPrevOrig() != e[8] || e[11]->GetPrevDest() != e[15] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[15])
		return false;
	NumTests += 1;
	return true;
}
double Vx(double x, double y)
{
	return 1.5;
}
double Vy(double x, double y)
{
	return 2.5;
}
bool tester_Triangulation_7(int& NumTests)
{
	{
		double Lx = 6, Ly = 6;
		Triangulation T = Triangulation::OffDiagonalGrid(1, 1, Lx, Ly);
		T.PrintTriangulation();
		EdgeIterator ite(T.GetBoundary());
		Edge* e = ite.Next();
		int counter = 0;
		while (e && counter < 10)
		{
			++counter;
			e = ite.Next();
		}
		if (counter != 4)
			return false;
	}
	{
		double Lx = 6, Ly = 6; 
		Triangulation T = Triangulation::OffDiagonalGrid(2, 2, Lx, Ly);
		EdgeIterator ite(T.GetBoundary());
		Edge* e = ite.Next();
		int counter = 0;
		while (e && counter < 100)
		{
			++counter;
			e = ite.Next();
		}
		if (counter != 8)
			return false; 
	}
	NumTests += 1;
	return true;
}
bool tester_Triangulation(int& NumTests)
{
	Vector3D A, B(1), C(0.5, 1), D(0.5, 0.5, 1);
	QuadEdge QE;
	QE.SetAsTetrahedron(A, B, C, D);
	if (!QE.TestIntegrity())
		return false;
	if (QE.NumVertices() != 4 || QE.NumEdges() != 6 || QE.NumFaces() != 4)
		return false;
	Triangulation T;
	T << QE;
	if (T.GetBoundary() != T.GetMesh2D()->GetFace(0))
		return false;
	if (T.GetMesh2D()->NumVertices() != 4 || T.GetMesh2D()->NumEdges() != 6 || T.GetMesh2D()->NumFaces() != 4)
		return false;
	if (!T.TestIntegrity())
		return false;
	Vector3D DD(0.5, 0.5), E(0.5,0.25), F(0.5,0.75);
	if (!T.isInside(DD) || !T.isOnOrInside(DD))
		return false;
	if (!T.isInside(E) || !T.isOnOrInside(E))
		return false;
	if (!T.isInside(F) || !T.isOnOrInside(F))
		return false;
	Vector3D H(0.5, 0);
	Vector3D J = B * 0.9 + C * 0.1;
	Vector3D K = C * 0.6 + A * 0.4;
	if (T.isInside(H) || !T.isOnOrInside(H))
		return false;
	if (T.isInside(J) || !T.isOnOrInside(J))
		return false;
	if (T.isInside(K) || !T.isOnOrInside(K))
		return false;
	Vector3D M = A * 0.5 + B * 1.5;
	Vector3D N(3, 3);
	if (T.isInside(M) || T.isOnOrInside(M))
		return false;
	if (T.isInside(N) || T.isOnOrInside(N))
		return false;
	//T.PrintTriangulation();
	if (!tester_Triangulation_1(NumTests))
		return false;
	if (!tester_Triangulation_2(NumTests))
		return false;
	if (!tester_Triangulation_3(NumTests))
		return false;
	if (!tester_Triangulation_4(NumTests))
		return false;
	if (!tester_Triangulation_5(NumTests))
		return false;
	if (!tester_Triangulation_6(NumTests))
		return false;
	if (!tester_Triangulation_7(NumTests))
		return false;
	NumTests += 1;
	return true;
}
#endif