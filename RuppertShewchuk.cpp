#include <list>
#include <utility> //for std::pair
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include "RuppertShewchuk.h"
#include "DelaunayLifting.h"
#include "QuadEdge.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "EdgeIterator.h"
#include "Edge.h"
#include "Face.h"
#include "Vertex.h"
#include "MathUtils.h"
#include "bPolygon.h"
double RuppertShewchuk::B = sqrt(2.0);
double RuppertShewchuk::minLengthToMinSegment = 0;
RuppertShewchuk::RuppertShewchuk()
{
	this->triang = 0;
	this->LengthScale = 0;
}
RuppertShewchuk::~RuppertShewchuk()
{
	this->triang = 0;
}
void RuppertShewchuk::RemovePoint(const Vector3D& point)
{
	int removed_i = -1;
	for (int i = 0; i < this->points.size(); ++i)
	{
		if (this->points[i] == point)
		{
			this->points.erase(this->points.begin() + i);
			removed_i = i;
			break;
		}
	}
	if (removed_i != -1)
	{
		for (int j = 0; j < this->constrainst.size(); ++j)
		{
			if (this->constrainst[j].point1 >= removed_i)
				--this->constrainst[j].point1;
			if (this->constrainst[j].point2 >= removed_i)
				--this->constrainst[j].point2;
		}
	}
}
bool RuppertShewchuk::IsEncroched(Edge* e, Face* fb)
{
	bool encroched = false;
	Vertex* vO = e->GetOrig();
	Vertex* vD = e->GetDest();
	Vector3D pO = vO->GetPoint();
	Vector3D pD = vD->GetPoint();
	Vector3D pM = (pO + pD) / 2.0;
	double R = pO.distance(pD) / 2.0;
	if (e->GetRight() != fb) {
		Vertex* vOR = this->triang->GetRightFaceOpposite(e);
		Vector3D pOR = vOR->GetPoint();
		double r = pOR.distance(pM);
		if (r <= R /*- 1.0e-10*/)
			encroched = true;
	}
	if (!encroched && e->GetLeft() != fb) {
		Vertex* vOL = this->triang->GetLeftFaceOpposite(e);
		Vector3D pOL = vOL->GetPoint();
		double r = pOL.distance(pM);
		if (r <= R/* - 1.0e-10*/)
			encroched = true;
	}
	return encroched;
}
void RuppertShewchuk::ReBuildInput(QuadEdge* QE, Face* fb, list<Edge*>& newConstraints)
{
	int Nv = QE->NumVertices();
	this->points.resize(0);
	this->segmentPoints.resize(0);
	unordered_set<Vertex*> boundaryVertices;
	boundaryVertices.rehash(Nv);
	unordered_map<Vertex*, int> VertexLocation;
	VertexLocation.rehash(Nv);
	VertexIterator itvb(fb);
	Vertex* vb = itvb.Next();
	int counter = 0;
	while (vb)
	{
		boundaryVertices.insert(vb);
		this->points.push_back(vb->GetPoint());
		this->segmentPoints.push_back(counter);
		VertexLocation[vb] = counter;
		++counter;
		vb = itvb.Next();
	}
	VertexIterator itv(QE);
	Vertex* v = itv.Next();
	while (v)
	{
		if (boundaryVertices.find(v) == boundaryVertices.end())
		{
			this->points.push_back(v->GetPoint());
			VertexLocation[v] = counter;
			++counter;
		}
		v = itv.Next();
	}
	this->constrainst.resize(0);
	for (auto it = newConstraints.begin(); it != newConstraints.end(); ++it) {
		Edge* e = *it;
		Vertex* vO = e->GetOrig();
		Vertex* vD = e->GetDest();
		TriangulationConstraint constraint;
		constraint.point1 = VertexLocation[vO];
		constraint.point2 = VertexLocation[vD];
		this->constrainst.push_back(constraint);
	}
	this->triang.reset();
}
Vertex* RuppertShewchuk::findVertex(Vector3D* newPoint)
{
	if (newPoint && this->triang)
	{
		VertexIterator itv(this->triang->GetMesh2D());
		Vertex* v = itv.Next();
		while (v)
		{
			if (v->GetPoint() == *newPoint)
				return v;
			v = itv.Next();
		}
	}
	return nullptr;
}
bool RuppertShewchuk::IsAcceptable(Vector3D* P)
{
	//Check to see if the new point is inside boundary
	if (P && this->triang)
	{
		return this->triang->isInside(*P);
	}
	return false;
}
int RuppertShewchuk::UnEncrochEdges(Vector3D* newPoint)
{
	int NumNewPointInserted = 0;
	if (!this->triang)
	{
		this->triang = make_unique<Triangulation>();
		*(this->triang) = DelaunayLifting::Triangulate(this->points, this->segmentPoints, this->constrainst);
	}
	Vertex* newVertex = this->findVertex(newPoint);
	QuadEdge* QE = this->triang->GetMesh2D();
	Face* fb = this->triang->GetBoundary();
	int Ne = QE->NumEdges();
	int Nv = QE->NumVertices();
	int Nf = QE->NumFaces();
	EdgeIterator ite(fb);
	Edge* e = ite.Next();
	list<Edge*> encrochedEdges;
	while (e){
		double eL = e->GetVector().abs();
		if (eL > this->LengthScale * RuppertShewchuk::minLengthToMinSegment)
		{
			if (!newVertex || this->triang->GetLeftFaceOpposite(e) == newVertex
				|| this->triang->GetRightFaceOpposite(e) == newVertex)
			{
				bool encroched = this->IsEncroched(e, fb);
				if (encroched) {
					++NumNewPointInserted;
					encrochedEdges.push_back(e);
				}
			}
		}
		e = ite.Next();
	}
	list<Edge*> encrochedConstraints;
	list<Edge*> newConstraints;
	for (int i = 0; i < this->constrainst.size(); ++i)//need O(n) optimization
	{
		int index0 = this->constrainst[i].point1;
		int index1 = this->constrainst[i].point2;
		Vertex* v0 = QE->GetVertexOfIndex(index0);//O(n)
		Vertex* v1 = QE->GetVertexOfIndex(index1);//O(n)
		if (v0 && v1)
		{
			Edge* e = QE->GetEdgeConnecting(v0, v1);
			if (e)
			{
				if (!newVertex || this->triang->GetLeftFaceOpposite(e) == newVertex
					|| this->triang->GetRightFaceOpposite(e) == newVertex)
				{
					bool encroched = this->IsEncroched(e, fb);
					if (encroched) {
						++NumNewPointInserted;
						encrochedConstraints.push_back(e);
					}
					else
					{
						newConstraints.push_back(e);
					}
				}
				else
				{
					newConstraints.push_back(e);
				}
			}
		}
	}
	if (NumNewPointInserted)
	{
		for (auto it = encrochedEdges.begin(); it != encrochedEdges.end(); ++it) {
			Edge* encrochingEdge = *it;
			QE->InsertVertexAtMid(encrochingEdge);
		}
		for (auto it = encrochedConstraints.begin(); it != encrochedConstraints.end(); ++it) {
			Edge* encrochingEdge = *it;
			Vertex* vMid = QE->InsertVertexAtMid(encrochingEdge);
			Edge* eNew = encrochingEdge->GetNext(vMid);
			newConstraints.push_back(encrochingEdge);
			newConstraints.push_back(eNew);
		}
		this->ReBuildInput(QE, fb, newConstraints);
	}
	return NumNewPointInserted;
}
Vector3D* RuppertShewchuk::SplitOneSkinnyTriangle()
{
	Vector3D* newPoint = 0;
	if (!this->triang)
	{
		this->triang = make_unique<Triangulation>();
		*(this->triang) = DelaunayLifting::Triangulate(this->points, this->segmentPoints, this->constrainst);
	}
	QuadEdge* QE = this->triang->GetMesh2D();
	Face* fb = this->triang->GetBoundary();
	FaceIterator itf(QE);
	Face* f = itf.Next();
	Face* MostSkinny = 0;
	Circle MostSkinnyCircle;
	double max_R_over_minL = 0;
	while (f){
		if (f != fb){
			Circle Cir = this->triang->GetCircumCircle(f);
			double R = Cir.GetRadius();
			double minL = this->triang->GetMinLength(f);
			double R_over_minL = R / minL;
			if (R_over_minL > RuppertShewchuk::B && R_over_minL > max_R_over_minL){
				Vector3D* newPointCandidate = new Vector3D(Cir.GetCenter());
				if (this->IsAcceptable(newPointCandidate))
				{
					max_R_over_minL = R_over_minL;
					MostSkinny = f;
					MostSkinnyCircle = Cir;
				}
			}
		}
		f = itf.Next();
	}
	if (MostSkinny)
	{
		newPoint = new Vector3D(MostSkinnyCircle.GetCenter());
		this->triang.reset();
	}
	return newPoint;
}
void RuppertShewchuk::UpdateLengthScale()
{
	this->LengthScale = -1;
	for (int i = 0; i < this->segmentPoints.size(); ++i) {
		int im1 = i ? i - 1 : this->segmentPoints.size() - 1;
		int n0 = this->segmentPoints[im1];
		int n1 = this->segmentPoints[i];
		Vector3D P0 = this->points[n0];
		Vector3D P1 = this->points[n1];
		Vector3D P01 = P1 - P0;
		double L = P01.abs();
		if (L < this->LengthScale || this->LengthScale < 0)
			this->LengthScale = L;
	}
}
void RuppertShewchuk::Build(const vector<Vector3D>& BoundaryPoints, int maxIter)
{
	this->UpdateLengthScale();
	this->triang.reset();
	this->triang = make_unique<Triangulation>();
	*(this->triang) = DelaunayLifting::Triangulate(this->points, this->segmentPoints, this->constrainst);
	int NumNewPointInserted = 0;
	int iter = 1; 
	do
	{
		++iter; 
		NumNewPointInserted = this->UnEncrochEdges();
		//Check for skinny triangles
		bool go_on;
		do{
			if (iter > maxIter)
				break;
			++iter; 
			go_on = false;
			Vector3D* newPoint = this->SplitOneSkinnyTriangle();
			if (newPoint){
				go_on = true;
				this->points.push_back(*newPoint);
				int NumEncrochedEdges = this->UnEncrochEdges(newPoint);
				if (!NumEncrochedEdges){//point is accepted
					++NumNewPointInserted;
				}
				else
				{//point is rejected
					this->RemovePoint(*newPoint);
					this->triang.reset();
				}
				NumNewPointInserted += NumEncrochedEdges;
				delete newPoint;
			}
		} while (go_on && iter < maxIter);
	} while (NumNewPointInserted && iter < maxIter);
	//Build the triangulation based on final points
	this->triang.reset();
	this->triang = make_unique<Triangulation>();
	*(this->triang) = DelaunayLifting::Triangulate(this->points, this->segmentPoints, this->constrainst);
}
Triangulation RuppertShewchuk::Tessellate(const vector<Vector3D>& BoundaryPoints, int maxIter)//static
{
	RuppertShewchuk Tss;
	Tss.points.resize(0);
	Tss.segmentPoints.resize(0);
	Tss.points.resize(BoundaryPoints.size());
	Tss.segmentPoints.resize(BoundaryPoints.size());
	for (int i = 0; i < BoundaryPoints.size(); ++i)
	{
		Tss.points[i] = BoundaryPoints[i];
		Tss.segmentPoints[i] = i;
	}
	Tss.Build(BoundaryPoints, maxIter);
	return *(Tss.triang);
}
Triangulation RuppertShewchuk::Tessellate(const TriangulationInput& input, int maxIter)
{
	
	RuppertShewchuk Tss;
	Tss.points.resize(0);
	Tss.segmentPoints.resize(0);
	Tss.points.resize(input.points.size());
	for (int i = 0; i < input.points.size(); ++i)
	{
		Tss.points[i] = input.points[i];
	}
	Tss.segmentPoints.resize(input.boundary.size());
	for (int i = 0; i < input.boundary.size(); ++i)
	{
		Tss.segmentPoints[i] = input.boundary[i];
	}
	Tss.constrainst.resize(input.constraints.size());
	for (int i = 0; i < input.constraints.size(); ++i)
	{
		Tss.constrainst[i].point1 = input.constraints[i].point1;
		Tss.constrainst[i].point2 = input.constraints[i].point2;
	}
	Tss.Build(input.points, maxIter);
	return *(Tss.triang);
}
bool tester_RuppertShewchuk(int& NumTests)
{
	if (!tester_RuppertShewchuk_1(NumTests))
		return false;
	if (!tester_RuppertShewchuk_2(NumTests))
		return false;
	if (!tester_RuppertShewchuk_3(NumTests))
		return false;
	if (!tester_RuppertShewchuk_4(NumTests))
		return false;
	if (!tester_RuppertShewchuk_5(NumTests))
		return false;
	if (!tester_RuppertShewchuk_6(NumTests))
		return false;
	if (!tester_RuppertShewchuk_7(NumTests))
		return false;
	if (!tester_RuppertShewchuk_8(NumTests))
		return false;
	if (!tester_RuppertShewchuk_10(NumTests))
		return false;
	if (!tester_RuppertShewchuk_11(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_RuppertShewchuk_1(int& NumTests)
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
	Triangulation T = RuppertShewchuk::Tessellate(input, 10);
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
bool tester_RuppertShewchuk_2(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 3;
	int Ny = 3;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = RuppertShewchuk::Tessellate(input, 10);
	T.PrintTriangulation();
	ofstream file;
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
bool tester_RuppertShewchuk_3(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 5;
	int Ny = 5;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = RuppertShewchuk::Tessellate(input, 10);
	T.PrintTriangulation();
	ofstream file;
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
bool tester_RuppertShewchuk_4(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 5;
	int Ny = 5;
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
	input.erase(input.begin() + 19);
	input.erase(input.begin() + 18);
	input.erase(input.begin() + 17);
	input.erase(input.begin() + 16);
	input.erase(input.begin() + 14);
	input.erase(input.begin() + 13);
	input.erase(input.begin() + 12);
	input.erase(input.begin() + 11);
	input.erase(input.begin() + 9);
	input.erase(input.begin() + 8);
	input.erase(input.begin() + 7);
	input.erase(input.begin() + 6);
	Triangulation T = RuppertShewchuk::Tessellate(input, 100);
	T.PrintTriangulation();
	ofstream file;
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
bool tester_RuppertShewchuk_5(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 5;
	int Ny = 5;
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
	input.erase(input.begin() + 19);
	input.erase(input.begin() + 18);
	input.erase(input.begin() + 17);
	input.erase(input.begin() + 16);
	input[13](1) *= 2.0;
	input.erase(input.begin() + 11);
	input.erase(input.begin() + 9);
	input.erase(input.begin() + 8);
	input.erase(input.begin() + 7);
	input.erase(input.begin() + 6);
	Triangulation T = RuppertShewchuk::Tessellate(input, 100);
	T.PrintTriangulation();
	ofstream file;
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
bool tester_RuppertShewchuk_6(int& NumTests)
{
	//Like test 5 but  input[2:3] are shifted upward
	vector<Vector3D> input;
	int Nx = 5;
	int Ny = 5;
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
	input.erase(input.begin() + 19);
	input.erase(input.begin() + 18);
	input.erase(input.begin() + 17);
	input.erase(input.begin() + 16);
	input[13](1) *= 2.0;
	input.erase(input.begin() + 11);
	input.erase(input.begin() + 9);
	input.erase(input.begin() + 8);
	input.erase(input.begin() + 7);
	input.erase(input.begin() + 6);
	input[2](1) = 0.5;
	input[3](1) = 0.5;
	Triangulation T = RuppertShewchuk::Tessellate(input, 10);
	T.PrintTriangulation();
	ofstream file;
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	++NumTests;
	return true;
}
bool tester_RuppertShewchuk_7(int& NumTests)
{
	TriangulationInput input;
	input.points.resize(8);
	input.points[0](0) = 0;			input.points[0](1) = 0;
	input.points[1](0) = 1;			input.points[1](1) = 0;
	input.points[2](0) = 2;			input.points[2](1) = 0;
	input.points[3](0) = 2;			input.points[3](1) = 1;
	input.points[4](0) = 2;			input.points[4](1) = 2;
	input.points[5](0) = 1;			input.points[5](1) = 2;
	input.points[6](0) = 0;			input.points[6](1) = 2;
	input.points[7](0) = 0;			input.points[7](1) = 1;
	input.boundary = { 0,1,2,3,4,5,6,7 };
	TriangulationConstraint c;
	c.point1 = 4; c.point2 = 7; input.constraints.push_back(c);
	Triangulation T = RuppertShewchuk::Tessellate(input, 0);//UNDER DEVELOPMENT !!!!!!!!!!!!!!!!!!!!!!!
	T.PrintTriangulation();
	ofstream file;
	if (!T.TestIntegrity())
		return false;
	/*if (!T.TestDelaunay())
		return false;*/
	++NumTests;
	return true;
}
bool tester_RuppertShewchuk_8(int& NumTests)
{//Shewchuk 2002 examole on page  748, same input as in tester_DelaunayLifting_10(...)
	TriangulationInput input;
	input.points.resize(10);
	input.points[0](0) = 3.02;		input.points[0](1) = 0;
	input.points[1](0) = 10;		input.points[1](1) = 0;
	input.points[2](0) = 10;		input.points[2](1) = 10;
	input.points[3](0) = 3.04;		input.points[3](1) = 10;
	input.points[4](0) = 0.82;		input.points[4](1) = 7.52;
	input.points[5](0) = 0;			input.points[5](1) = 6.01;
	input.points[6](0) = 0.52;		input.points[6](1) = 3.53;
	input.points[7](0) = 1.52;		input.points[7](1) = 4.77;
	input.points[8](0) = 1.52;		input.points[8](1) = 7.19;
	input.points[9](0) = 2.02;		input.points[9](1) = 5.49;
	input.boundary = { 0,1,2,3,4,5,6 };
	TriangulationConstraint c;
	c.point1 = 6; c.point2 = 7; input.constraints.push_back(c);
	c.point1 = 7; c.point2 = 8; input.constraints.push_back(c);
	c.point1 = 8; c.point2 = 4; input.constraints.push_back(c);
	Triangulation T = RuppertShewchuk::Tessellate(input, 1000);
	T.PrintTriangulation();
	ofstream file;
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	double alpha_max = T.GetMaxAngle() / pi * 180.0;
	if (fabs(alpha_max - 120.9729) > 0.01)
		return false;
	int Nv = T.NumVertices();
	if (Nv != 20)
		return false;
	++NumTests;
	return true;
}
bool tester_RuppertShewchuk_10(int& NumTests)
{
	int N = 3;
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
	RuppertShewchuk::minLengthToMinSegment = 0.1;
	Triangulation Triang = RuppertShewchuk::Tessellate(input, 1000);
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	if (fabs(alpha_max - 116.56) > 0.01)
		return false;
	int Nv = Triang.NumVertices();
	if (Nv != 36)
		return false;
	++NumTests;
	return true;
}
bool tester_RuppertShewchuk_11(int& NumTests)//shows an unfixed bug
{
	int N = 1;
	vector<Vector3D> input;
	Vector3D A, B(1, 0), C;
	double A_angle = 5.0;
	double B_angle = 180.0 - 2.0 * A_angle;
	double B_prime = 2.0 * A_angle;
	double B_prime_Rad = 2.0 * A_angle / 180.0 * pi;
	C(0) = 1.0 + cos(B_prime_Rad); C(1) = sin(B_prime_Rad);
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
	RuppertShewchuk::minLengthToMinSegment = 0.1;
	Triangulation Triang = RuppertShewchuk::Tessellate(input, 10);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//It will triangulate if iterations are increased to 1000. 
	//However, the shape of the resulting triangulation is suspicious.
	//it might have very skinny triangles
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	/*if (fabs(alpha_max - 116.56) > 0.01)
		return false;*/
	int Nv = Triang.NumVertices();
	/*if (Nv != 36)
		return false;*/
	++NumTests;
	return true;
}
bool tester_RuppertShewchuk_12(int& NumTests)//Gosper Island of order 2
{
	int N = 2;//number of point on a side
	vector<Vector3D> input; 
	double a = 1.0;
	double perturbationRed = a / 2.0;
	int rm = 100;//Random number generation rX[n+1] = (ra * rX[n] + rb) % rm; 
	int ra = 25;
	int rb = 47;
	int rX = 0;
	bMatrix M(3, 3);
	M(0, 0) = cos(pi / 3.0); M(0, 1) = -sin(pi / 3.0);
	M(1, 0) = sin(pi / 3.0); M(1, 1) = cos(pi / 3.0);
	M(2, 2) = 1.0;
	Vector3D T[3];
	T[0](0) = a;
	T[1] = M * (T[0]);
	T[2] = M * (T[1]);
	Vector3D P(3.5 * a, sin(pi / 3.0) * a);
	for (int i = 0; i < 6; ++i)
	{
		for (int k = 0; k < 3; ++k)
		{
			Vector3D dP = T[k];
			dP(0) += (double(rX) / double(rm) - 0.5) * perturbationRed;
			rX = (ra * rX + rb) % rm;
			dP(1) += (double(rX) / double(rm) - 0.5) * perturbationRed;
			rX = (ra * rX + rb) % rm;
			dP = dP / double(N);
			for (int j = 0; j < N; ++j)
			{
				input.push_back(P);
				P = P + dP;
			}
			T[k] = M * T[k];
		}
	}

	Triangulation Triang = RuppertShewchuk::Tessellate(input, 1000);
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	/*if (fabs(alpha_max - 120.9729) > 0.01)
		return false;*/
	int Nv = Triang.NumVertices();
	/*if (Nv != 20)
		return false;*/
	++NumTests;
	return true;
}
bool tester_RuppertShewchuk_13(int& NumTests)
{
	int N = 1;//number of point on a side
	double R1 = 1.0, R1a = 1.25, R2 = 2.0, R2a = 1.75;
	vector<Vector3D> input;
	Vector3D P, A, B, dP;
	for (int i = 0; i < N; ++i) {
		P(0) = R1 + (R2 - R1) / double(N) * double(i);
		P(1) = 0;
		input.push_back(P);
	}
	int N0 = input.size();
	for (int i = 0; i < N; ++i) {
		double theta = 3.0 * pi / 16.0 / double(N) * double(i);
		P(0) = R2 * cos(theta);
		P(1) = R2 * sin(theta);
		input.push_back(P);
	}
	A(0) = R2 * cos(3.0 * pi / 16.0); A(1) =  R2 * sin(3.0 * pi / 16.0);
	B(0) = R2a * cos(3.5 * pi / 16.0); B(1) = R2a * sin(3.5 * pi / 16.0);
	dP = (B - A);
	dP = dP / double(N);
	P = A;
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	A = B;
	B(0) = R2 * cos(pi / 4.0); B(1) = R2 * sin(pi / 4.0);
	dP = (B - A);
	dP = dP / double(N);
	P = A;
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	int N1 = input.size();
	bMatrix M(3, 3);
	M(0, 0) = cos(pi / 4.0); M(0, 1) = -sin(pi / 4.0);
	M(1, 0) = sin(pi / 4.0); M(1, 1) = cos(pi / 4.0);
	M(2, 2) = 1.0;
	for (int i = N0; i < N1; ++i) {
		P = M * input[i];
		input.push_back(P);
	}
	for (int i = 0; i < N; ++i) {
		P(0) = 0;
		P(1) = R2 - (R2 - R1) / double(N) * double(i);
		input.push_back(P);
	}
	int N2 = input.size();
	for (int i = 0; i < N; ++i) {
		double theta = pi / 2.0 - (pi / 2.0 - 6.0 * pi / 16.0) / double(N) * double(i);
		P(0) = R1 * cos(theta);
		P(1) = R1 * sin(theta);
		input.push_back(P);
	}
	A(0) = R1 * cos(6.0 * pi / 16.0); A(1) = R1 * sin(6.0 * pi / 16.0);
	B(0) = R1a * cos(5.0 * pi / 16.0); B(1) = R1a * sin(5.0 * pi / 16.0);
	dP = (B - A);
	dP = dP / double(N);
	P = A;
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	A = B;
	B(0) = R1 * cos(pi / 4.0); B(1) = R1 * sin(pi / 4.0);
	dP = (B - A);
	dP = dP / double(N);
	P = A;
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	bMatrix M2(3, 3);
	M2(0, 0) = cos(pi / 4.0); M2(0, 1) = sin(pi / 4.0);
	M2(1, 0) = -sin(pi / 4.0); M2(1, 1) = cos(pi / 4.0);
	M2(2, 2) = 1.0;
	int N3 = input.size();
	for (int i = N2; i < N3; ++i) {
		P = M2 * input[i];
		input.push_back(P);
	}
	Triangulation Triang = RuppertShewchuk::Tessellate(input, 1000);
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	/*if (fabs(alpha_max - 120.9729) > 0.01)
		return false;*/
	int Nv = Triang.NumVertices();
	/*if (Nv != 20)
		return false;*/
	++NumTests;
	return true;
}

