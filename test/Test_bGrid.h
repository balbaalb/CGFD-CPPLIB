#ifndef Test_bGridH
#define Test_bGridH
#include "../src/bGrid.h"
#include "../src/QuadEdge.h"
bool tester_bGrid(int& NumTests)
{
	const int Nx = 5;
	const int Ny = 7;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 2 * Nx * Ny + Nx + Ny;
	int Nf = Nx * Ny + 1;
	bGrid g(Nx, Ny, 3.0, 4.0);
	if (g.GetMesh2D()->NumVertices() != Nv || g.GetMesh2D()->NumEdges() != Ne || g.GetMesh2D()->NumFaces() != Nf)
		return false;
	QuadEdge* qe = g.GetMesh2D();
	Vertex* v[48];
	for (int i = 0; i < 48; i++)
	{
		v[i] = qe->GetVertex(i);
	}
	Edge* e[82];
	for (int i = 0; i < 82; i++)
	{
		e[i] = qe->GetEdge(i);
	}
	Face* f[35];
	for (int i = 0; i < 35; i++)
	{
		f[i] = qe->GetFace(i);
	}
	if (!g.IsHorizontal(e[11]) || !g.IsHorizontal(e[12]) || !g.IsHorizontal(e[13]))
		return false;
	Edge* e12 = e[12];
	Edge* e13 = g.GetOppositeEdge(e12, f[11]);
	Edge* e11 = g.GetOppositeEdge(e12, f[10]);
	if (!e13 || e13 != e[13]  || !e11 || e11 != e[11])
		return false;
	if (e12 != g.GetOppositeEdge(e13, f[11]) || e12 != g.GetOppositeEdge(e11, f[10]))
		return false;
	Edge* e57 = e[57];
	Edge* e64 = e[64];
	Edge* e71 = e[71];
	if (!e57 || !e64 || !e71)
		return false;
	if (g.IsHorizontal(e57) || g.IsHorizontal(e64) || g.IsHorizontal(e71))
		return false;
	if (e57 != g.GetOppositeEdge(e64, f[17]) || e71 != g.GetOppositeEdge(e64, f[24]))
		return false;
	if (e64 != g.GetOppositeEdge(e57, f[17]) || e64 != g.GetOppositeEdge(e71, f[24]))
		return false;
	++NumTests;
	return true;
}
#endif