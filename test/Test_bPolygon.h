#ifndef Test_bPolygonH
#define Test_bPolygonH
#include "../src/bPolygon.h"
#include "../src/QuadEdge.h"
#include "../src/Edge.h"
#include "../src/VertexIterator.h"
#include "../src/EdgeIterator.h"
bool tester_bPolygon_1(int& NumTests)
{
	Vector3D A, B(1), C(2, 2), D(0, 1);
	bPolygon PPP;
	PPP.AddVertex(A);
	PPP.AddVertex(B);
	PPP.AddVertex(C);
	PPP.AddVertex(D);
	PPP.Close();
	if (!PPP.TestIntegrity())
		return false;
	bPolygon PP(PPP), P;
	P = PP;
	Face* inface = P.GetInsideFace();
	if (!inface)
		return false;
	VertexIterator itv(inface);
	EdgeIterator ite(inface);
	const int N = 4;
	Vertex* v[N];
	Edge* e[4];
	for (int i = 0; i < N; ++i){
		v[i] = 0;
		v[i] = itv.Next();
		e[i] = 0;
		e[i] = ite.Next();
		if (e[i]->GetLeft() != inface)
			return false;
	}
	if (itv.Next() || !v[N - 1])
		return false;
	if (ite.Next() || !e[N - 1])
		return false;
	if (v[0]->GetPoint() != D || v[1]->GetPoint() != A || v[2]->GetPoint() != B || v[3]->GetPoint() != C)
		return false;
	if (!P.isEar(v[0]) || !P.isEar(v[1]) || !P.isEar(v[2]) || !P.isEar(v[3])) 
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	Vector3D E(0.5, 0.5), F(3, 3), G(3, 0), H(0.5, 0);
	if (!P.isInside(E) || !P.isOnOrInside(E))
		return false;
	if (P.isInside(F) || P.isOnOrInside(F))
		return false;
	if (P.isInside(G) || P.isOnOrInside(G))
		return false;
	if (P.isInside(H) || !P.isOnOrInside(H))
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_2(int& NumTests)
{
	Vector3D A, B(1), C(2, 2), D(0, 1);
	bPolygon PPP;
	PPP.AddVertex(A);
	PPP.AddVertex(D);
	PPP.AddVertex(C);
	PPP.AddVertex(B);
	PPP.Close();
	if (!PPP.TestIntegrity())
		return false;
	bPolygon PP(PPP), P;
	P = PP;
	Face* inface = P.GetInsideFace();
	if (!inface)
		return false;
	VertexIterator itv(inface);
	EdgeIterator ite(inface);
	const int N = 4;
	Vertex* v[N];
	Edge* e[4];
	for (int i = 0; i < N; ++i){
		v[i] = 0;
		v[i] = itv.Next();
		e[i] = 0;
		e[i] = ite.Next();
		if (e[i]->GetLeft() != inface)
			return false;
	}
	if (itv.Next() || !v[N - 1])
		return false;
	if (ite.Next() || !e[N - 1])
		return false;
	if (v[0]->GetPoint() != A || v[1]->GetPoint() != B || v[2]->GetPoint() != C || v[3]->GetPoint() != D)
		return false;
	if (!P.isEar(v[0]) || !P.isEar(v[1]) || !P.isEar(v[2]) || !P.isEar(v[3]))
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_3(int& NumTests)
{
	Vector3D A, B(3,3), C(0, 1), D(-3, 3);
	bPolygon PPP;
	PPP.AddVertex(A);
	PPP.AddVertex(D);
	PPP.AddVertex(C);
	PPP.AddVertex(B);
	PPP.Close();
	if (!PPP.TestIntegrity())
		return false;
	bPolygon PP(PPP), P;
	P = PP;
	Face* inface = P.GetInsideFace();
	if (!inface)
		return false;
	VertexIterator itv(inface);
	EdgeIterator ite(inface);
	const int N = 4;
	Vertex* v[N];
	Edge* e[4];
	for (int i = 0; i < N; ++i){
		v[i] = 0;
		v[i] = itv.Next();
		e[i] = 0;
		e[i] = ite.Next();
		if (e[i]->GetLeft() != inface)
			return false;
	}
	if (itv.Next() || !v[N - 1])
		return false;
	if (ite.Next() || !e[N - 1])
		return false;
	if (v[0]->GetPoint() != A || v[1]->GetPoint() != B || v[2]->GetPoint() != C || v[3]->GetPoint() != D)
		return false;
	if (P.isEar(v[0]) || !P.isEar(v[1]/*B*/) || P.isEar(v[2]/*C*/) || !P.isEar(v[3]/*D*/))
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_4(int& NumTests)
{
	vector<Vector3D> pentagonPoints;
	pentagonPoints.resize(5);
	bPolygon Pentagon;
	for (int i = 0; i < 5; ++i){
		double theta = double(i) * 72.0 / 180.0 * pi + pi / 2.0;
		double r = 1.0;
		pentagonPoints[i](0) = r * cos(theta);
		pentagonPoints[i](1) = r * sin(theta);
		Pentagon.AddVertex(pentagonPoints[i]);
	}
	Pentagon.Close();
	if (!Pentagon.TestIntegrity())
		return false;
	Triangulation T = Pentagon.Triangulate();
	if (!T.TestIntegrity())
		return false;
	int N = 5;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_5(int& NumTests)
{
	vector<Vector3D> hexagonPoints;
	int N = 6;
	hexagonPoints.resize(N);
	bPolygon Hexagon;
	for (int i = 0; i < N; ++i){
		double theta = double(i) * (2 * pi) / double(N);
		double r = 1.0;
		hexagonPoints[i](0) = r * cos(theta);
		hexagonPoints[i](1) = r * sin(theta);
		Hexagon.AddVertex(hexagonPoints[i]);
	}
	Hexagon.Close();
	if (!Hexagon.TestIntegrity())
		return false;
	Triangulation T = Hexagon.Triangulate();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_6(int& NumTests)
{
	vector<Vector3D> points;
	points.resize(11);
	points[0](0) = 0; points[0](1) = 0;
	points[1](0) = 1; points[1](1) = 0;
	points[2](0) = 2; points[2](1) = 1;
	points[3](0) = 2; points[3](1) = 2;
	points[4](0) = 3; points[4](1) = 2;
	points[5](0) = 3; points[5](1) = 3;
	points[6](0) = 0.5; points[6](1) = 3;
	points[7](0) = 0; points[7](1) = 1.5;
	points[8](0) = -1; points[8](1) = 1.5;
	points[9](0) = -1; points[9](1) = 0.75;
	points[10](0) = 0; points[10](1) = 0.75;
	bPolygon P;
	for (int i = 0; i < 11; ++i){
		P.AddVertex(points[i]);
	}
	P.Close();
	if (!P.TestIntegrity())
		return false;
	Triangulation T = P.Triangulate();
	if (!T.TestIntegrity())
		return false;
	int N = 11;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_7(int& NumTests)
{
	vector<Vector3D> points;
	int Nx = 2;
	int Ny = 2;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	points.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		points[i](0) = hx * i; points[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		points[n](0) = Lx; points[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		points[n](0) = Lx - i * hx; points[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		points[n](0) = 0; points[n](1) = Ly - i * hy;
	}
	int N = points.size();
	bPolygon P;
	for (int i = 0; i < N; ++i){
		P.AddVertex(points[i]);
	}
	P.Close();
	if (!P.TestIntegrity())
		return false;
	Triangulation T = P.Triangulate();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_8(int& NumTests)
{
	double h = sqrt(3) / 2.0;
	vector<Vector3D> points;
	points.resize(27);
	points[0](0) = 6;		points[0](1) = 5;
	points[1](0) = 6;		points[1](1) = 6;
	points[2](0) = 5;		points[2](1) = 6;
	points[3](0) = 4;		points[3](1) = 6;
	points[4](0) = 3;		points[4](1) = 6;
	points[5](0) = 2;		points[5](1) = 6;
	points[6](0) = 1;		points[6](1) = 6;
	points[7](0) = 0;		points[7](1) = 6;
	points[8](0) = 0;		points[8](1) = 5;
	points[9](0) = 0;		points[9](1) = 4;
	points[10](0) = 0;		points[10](1) = 3;
	points[11](0) = 0;		points[11](1) = 2;
	points[12](0) = 0;		points[12](1) = 1;
	points[13](0) = 1;		points[13](1) = 0;
	points[14](0) = 1.5;	points[14](1) = h;
	points[15](0) = 2;		points[15](1) = 0;
	points[16](0) = 2.5;	points[16](1) = h;
	points[17](0) = 3;		points[17](1) = 0;
	points[18](0) = 3.5;	points[18](1) = h;
	points[19](0) = 4;		points[19](1) = 0;
	points[20](0) = 4.5;	points[20](1) = h;
	points[21](0) = 6;		points[21](1) = 2;
	points[22](0) = 6-h;	points[22](1) = 2.5;
	points[23](0) = 6;		points[23](1) = 3;
	points[24](0) = 6-h;	points[24](1) = 3.5;
	points[25](0) = 6;		points[25](1) = 4;
	points[26](0) = 6-h;	points[26](1) = 4.5;
	int N = points.size();
	bPolygon P;
	for (int i = 0; i < N; ++i){
		P.AddVertex(points[i]);
	}
	P.Close();
	if (!P.TestIntegrity())
		return false;
	Triangulation T = P.Triangulate();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon(int& NumTests)
{
	if (!tester_bPolygon_1(NumTests))
		return false;
	if (!tester_bPolygon_2(NumTests))
		return false;
	if (!tester_bPolygon_3(NumTests))
		return false;
	if (!tester_bPolygon_4(NumTests))
		return false;
	if (!tester_bPolygon_5(NumTests))
		return false;
	if (!tester_bPolygon_6(NumTests))
		return false;
	if (!tester_bPolygon_7(NumTests))
		return false;
	if (!tester_bPolygon_8(NumTests))
		return false;
	++NumTests;
	return true;
}
#endif