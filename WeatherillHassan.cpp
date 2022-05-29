#include <utility> //for std::pair
#include "WeatherillHassan.h"
#include "DelaunayLifting.h"
#include "QuadEdge.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "EdgeIterator.h"
#include "Edge.h"
#include "Face.h"
#include "Vertex.h"
#include "MathUtils.h"
double WeatherillHassan::alpha = 0.5;
double WeatherillHassan::beta = 0.5;
WeatherillHassan::LENGTH_SCALE_DEFINITION WeatherillHassan::length_scale_definition = AVERAGE_ADJ_EDGE_LENGTHS;
WeatherillHassan::WeatherillHassan()
{
	this->triang = 0;
	this->numInputPoints = 0;
	this->scale0 = 0;
}
WeatherillHassan::~WeatherillHassan()
{
	this->numInputPoints = 0;
	this->scale0 = 0;
}
void WeatherillHassan::UpdateBoundary()
{
	QuadEdge* QE = this->triang->GetMesh2D();
	this->boundary.rehash(0);
	Face* fb = this->triang->GetBoundary();
	VertexIterator itvb(fb);
	Vertex* vb = itvb.Next();
	while (vb){
		this->boundary.insert(vb);
		vb = itvb.Next();
	}
}
void WeatherillHassan::UpdateLengthScales()
{
	//unordered_set<Vertex*>::iterator its;
	QuadEdge* QE = this->triang->GetMesh2D();
	int Nv = QE->NumVertices();
	this->LenScale.rehash(0);
	this->LenScale.rehash(Nv);
	VertexIterator itv(QE);
	Vertex* v = itv.Next();
	while (v){
		double lenScale = 0;
		EdgeIterator ite(v);
		Edge* e = ite.Next();
		int Ne = 0;
		while (e){
			Vector3D E = e->GetVector();
			double len = E.abs();
			if (length_scale_definition == MAX_ADJ_EDGE_LENGTH){
				if (len > lenScale)
					lenScale = len;
			}
			else if (length_scale_definition == AVERAGE_ADJ_EDGE_LENGTHS){
				++Ne;
				lenScale += len;
			}
			e = ite.Next();
		}
		if (length_scale_definition == AVERAGE_ADJ_EDGE_LENGTHS)
			lenScale /= Ne;
		//its = this->boundary.find(v);
		this->LenScale[v] = lenScale;
		v = itv.Next();
	}
	if (!scale0){
		scale0 = 1e10;
		unordered_map<Vertex*, double>::iterator itm;
		for (itm = this->LenScale.begin(); itm != this->LenScale.end(); ++itm){
			if (itm->second < scale0)
				scale0 = itm->second;
		}
	}
}
int WeatherillHassan::InsertVertices(vector<Vector3D>& points)
{
	QuadEdge* QE = this->triang->GetMesh2D();
	Face* fb = this->triang->GetBoundary();
	FaceIterator itf(QE);
	Face* f = itf.Next();
	int NumNewPointInserted = 0;
	while (f){
		if (f != fb){
			VertexIterator itv(f);
			Vertex* a = itv.Next();
			Vertex* b = itv.Next();
			Vertex* c = itv.Next();
			if (itv.Next() || a == b || a == c || b == c || !a || !b || !c)
				throw "WeatherillHassan::InsertVertices()";
			Vector3D A = a->GetPoint();
			Vector3D B = b->GetPoint();
			Vector3D C = c->GetPoint();
			double aScale = this->LenScale[a];
			double bScale = this->LenScale[b];
			double cScale = this->LenScale[c];
			Vector3D M = (A + B + C) / 3.0;
			double MA = A.distance(M);
			double MB = B.distance(M);
			double MC = C.distance(M);
			double scale = scale0;// (aScale + bScale + cScale) / 3.0;
			if (MA > alpha * scale && MB > alpha * scale && MC > alpha * scale)
			{
				bool accepted = true;
				for (int i = this->numInputPoints; i < points.size(); ++i){
					//Check with previously accepted new points
					double MS = points[i].distance(M);
					if (MS < beta * scale){
						accepted = false;
						break;
					}
				}
				if (accepted){
					points.push_back(M);
					++NumNewPointInserted;
				}
			}
		}
		f = itf.Next();
	}
	return NumNewPointInserted;
}
void WeatherillHassan::Build(const vector<Vector3D>& BoundaryPoints, int maxIter)
{
	this->numInputPoints = BoundaryPoints.size();
	vector<Vector3D> points;
	vector<int> boundary;
	points.resize(BoundaryPoints.size());
	boundary.resize(BoundaryPoints.size());
	for (int i = 0; i < BoundaryPoints.size(); ++i)
	{
		points[i] = BoundaryPoints[i];
		boundary[i] = i;
	}
	int NumNewPointInserted = 0;
	int iter = 0;
	do {
		++iter;
		this->triang.reset();
		this->triang = make_unique<Triangulation>();
		*(this->triang) = DelaunayLifting::Triangulate(points,boundary);
		this->UpdateBoundary();
		this->UpdateLengthScales();
		NumNewPointInserted = this->InsertVertices(points);
	} while (NumNewPointInserted && iter < maxIter);
}
Triangulation WeatherillHassan::Tessellate(const vector<Vector3D>& BoundaryPoints, int maxIter)//static
{
	WeatherillHassan Tss;
	Tss.Build(BoundaryPoints,maxIter);
	return *(Tss.triang);
}
bool test_WeatherillHassan(int& NumTests)
{
	if (!test_WeatherillHassan_1(NumTests))
		return false;
	if (!test_WeatherillHassan_2(NumTests))
		return false;
	if (!test_WeatherillHassan_3(NumTests))
		return false;
	if (!test_WeatherillHassan_4(NumTests))
		return false;
	NumTests += 1;
	return true;
}
bool test_WeatherillHassan_1(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(8);
	input[0](0) = 0; input[0](1) = 0;
	input[1](0) = 1.0; input[1](1) = 0;
	input[2](0) = 2.0; input[2](1) = 0;
	input[3](0) = 2.0; input[3](1) = 1.0;
	input[4](0) = 2.0; input[4](1) = 2.0;
	input[5](0) = 1.0; input[5](1) = 2.0;
	input[6](0) = 0; input[6](1) = 2.0;
	input[7](0) = 0; input[7](1) = 1.0;
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation T = WeatherillHassan::Tessellate(input,100);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	++NumTests;
	return true;
}
bool test_WeatherillHassan_2(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 4;
	int Ny = 4;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2*(Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation T = WeatherillHassan::Tessellate(input, 100);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	++NumTests;
	return true;
}
bool test_WeatherillHassan_3(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 4;
	int Ny = 4;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i) {
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i) {
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i) {
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i) {
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	//Making the boundary non-convex:----
	int n1 = Nx / 2;
	input[n1](1) = hy;
	int n2 = Nx + Ny / 2;
	input[n2](0) = Lx + hx;
	//-----------------------------------
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation T = WeatherillHassan::Tessellate(input, 100);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;/**/
	++NumTests;
	return true;
}
bool test_WeatherillHassan_4(int& NumTests)
{
	int N = 4;
	vector<Vector3D> input;
	Vector3D A, B(1, 0), C(2, 1);
	Vector3D P = A;
	Vector3D dP = B - A;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	P = B;
	dP = C - B;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	P = C;
	dP = A - C;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation Triang = WeatherillHassan::Tessellate(input, 100);
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	if (fabs(alpha_max - 135) > 0.01)
		return false;
	int Nv = Triang.NumVertices();
	if (Nv != 12)
		return false;
	++NumTests;
	return true;
}
