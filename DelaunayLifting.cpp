#include <unordered_set>
#include "DelaunayLifting.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include "QuadEdge.h"
#include "EdgeList.h"
#include "EdgeIterator.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "ConvexHull3D.h"
#include "Triangulation.h"
#include "LineSegment2D.h"
#include "bPolygon.h"
void DelaunayLifting::CopyBody(const DelaunayLifting& rhs)
{
	this->lifted = 0;
	if (rhs.lifted)
		this->lifted = new ConvexHull3D(*rhs.lifted);
	this->DelTriangulation = 0;
	if (rhs.DelTriangulation)
	{
		this->DelTriangulation = new Triangulation(*rhs.DelTriangulation);
	}
}
void DelaunayLifting::DeleteBody()
{
	delete this->lifted;
	this->lifted = 0;
	delete this->DelTriangulation;
	this->DelTriangulation = 0;
}
void DelaunayLifting::Lift(const vector<Vector3D>& input, vector<Vector3D>& liftedPoints)
{
	liftedPoints.resize(0);
	for (int i = 0; i < input.size(); ++i)
	{
		liftedPoints.push_back(input[i]);
		double x = liftedPoints[i](0);
		double y = liftedPoints[i](1);
		double lift = x * x + y * y;
		double perturbation = 1.0e-8 * x;
		liftedPoints[i](2) = -(lift + perturbation);
	}
}
void DelaunayLifting::UpdateFaceVisibility()//Mark faces that are visible from the top.
{
	QuadEdge* hull = this->lifted->GetHullPointer();
	FaceIterator itf(hull);
	FaceVisible* f = (FaceVisible*)itf.Next();
	while (f)
	{
		Vector3D n = hull->GetNormal(f);
		double nz = n(2);
		if (fabs(nz) < 1.0e-10)
			nz = 0;
		f->SetVisible(nz <= 0);
		f = (FaceVisible*)itf.Next();
	}
}
void DelaunayLifting::BuildMesh2D()
{
	delete this->DelTriangulation;
	this->DelTriangulation = new Triangulation;
	*(this->DelTriangulation) << *(this->lifted->GetHull());
	this->DelTriangulation->GetMesh2D()->UpdateFacesDegree();
	this->DelTriangulation->GetMesh2D()->UpdateVerticesDegree();
}
DelaunayLifting::DelaunayLifting()
{
	this->lifted = 0;
	this->DelTriangulation = 0;
}
DelaunayLifting::DelaunayLifting(const DelaunayLifting& rhs)
{
	this->CopyBody(rhs);
}
DelaunayLifting::~DelaunayLifting()
{
	this->DeleteBody();
}
DelaunayLifting& DelaunayLifting::operator=(const DelaunayLifting& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
	return *this;
}
void DelaunayLifting::Build(const vector<Vector3D>& input)
{
	this->DeleteBody();
	if (input.size() > 3)
	{
		this->lifted = new ConvexHull3D;
		vector<Vector3D> liftedPoints;
		DelaunayLifting::Lift(input, liftedPoints);
		this->lifted->Build(liftedPoints);
		this->UpdateFaceVisibility();
		this->lifted->UpdateEdgesVisibility();
		this->lifted->DeleteEdges();
		this->BuildMesh2D();
	}
	else if (input.size() == 3)
	{
		QuadEdge QE;
		QE.SetAsDoubleTriangle(input[0], input[1], input[2]);
		this->DelTriangulation = new Triangulation;
		*(this->DelTriangulation) << QE;
	}
}
void DelaunayLifting::Build(const vector<Vector3D>& input, const vector<int>& boundary)
{
	//it assumes that all points in input are inside boundary 
	this->Build(input);
	QuadEdge* qe = this->DelTriangulation->GetMesh2D();
	if (input.size() > 3 && boundary.size() > 3)
	{
		this->VertexHash.resize(0);
		this->VertexHash.resize(qe->NumVertices());
		VertexIterator itv(qe);
		Vertex* v = itv.Next();
		while (v) {
			this->VertexHash[v->index] = v;//Relies on the fact that the Vertex::index created from input[i] is equal to i. 
			v = itv.Next();
		}
		vector<Vertex*> boundaryVertices;
		for (int i = 0; i < boundary.size(); ++i)
		{
			int index = boundary[i];
			boundaryVertices.push_back(VertexHash[index]);
		}
		this->DelTriangulation->ImposeBoundary(boundaryVertices);
	}
}
void DelaunayLifting::Build(const vector<Vector3D>& input, const vector<int>& Boundary, 
	const vector<TriangulationConstraint>& Constraints)
{
	this->Build(input, Boundary);//Builds this->VertexHash
	if (input.size() > 3)
	{
		for (int i = 0; i < Constraints.size(); ++i)
		{
			int index0 = Constraints[i].point1;
			int index1 = Constraints[i].point2;
			Vertex* v0 = this->VertexHash[index0];
			Vertex* v1 = this->VertexHash[index1];
			this->DelTriangulation->Connect(v0,v1);
		}
	}
}
Triangulation DelaunayLifting::GetTriangulation() const
{
	return *(this->DelTriangulation);
}
Triangulation DelaunayLifting::Triangulate(const vector<Vector3D>& inputPoints)//static
{
	DelaunayLifting meshBuilder;
	meshBuilder.Build(inputPoints);
	return *(meshBuilder.DelTriangulation);
}
Triangulation DelaunayLifting::Triangulate(const vector<Vector3D>& inputPoints, const vector<int>& Boundary)
{
	DelaunayLifting meshBuilder;
	meshBuilder.Build(inputPoints, Boundary);
	return *(meshBuilder.DelTriangulation);
}
Triangulation DelaunayLifting::Triangulate(const vector<Vector3D>& inputPoints, const vector<int>& Boundary,
	const vector<TriangulationConstraint>& Constraints)
{
	DelaunayLifting meshBuilder;
	meshBuilder.Build(inputPoints, Boundary, Constraints);
	return *(meshBuilder.DelTriangulation);
}
Triangulation DelaunayLifting::Triangulate(const TriangulationInput& input)
{
	DelaunayLifting meshBuilder;
	meshBuilder.Build(input.points, input.boundary, input.constraints);
	return *(meshBuilder.DelTriangulation);
}
bool tester_DelaunayLifting(int& NumTests)
{
	if (!tester_DelaunayLifting_1(NumTests))
		return false;
	if (!tester_DelaunayLifting_2(NumTests))
		return false;
	if (!tester_DelaunayLifting_3(NumTests))
		return false;
	if (!tester_DelaunayLifting_4(NumTests))
		return false;
	if (!tester_DelaunayLifting_5(NumTests))
		return false;
	if (!tester_DelaunayLifting_6(NumTests))
		return false;
	if (!tester_DelaunayLifting_7(NumTests))
		return false;
	if (!tester_DelaunayLifting_8(NumTests))
		return false;
	if (!tester_DelaunayLifting_9(NumTests))
		return false;
	if (!tester_DelaunayLifting_10(NumTests))
		return false;
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_1(int& NumTests)
{//SINGLE TRIANGLE
	vector<Vector3D> input;
	input.resize(3);
	input[0](0) = 0; input[0](1) = 0;
	input[1](0) = 1; input[1](1) = 0;
	input[2](0) = 0; input[2](1) = 1;
	Triangulation triangulation = DelaunayLifting::Triangulate(input);
	//triangulation.PrintTriangulation();
	if (!triangulation.TestIntegrity())
		return false;
	if (!triangulation.TestDelaunay())
		return false;
	QuadEdge* mesh = triangulation.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	if (mesh->NumVertices() != 3 || mesh->NumEdges() != 3 || mesh->NumFaces() != 2)
		return false;
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_2(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(4);
	input[1](0) = 1.0;
	input[2](0) = input[2](1) = 1.0 / sqrt(2.0);
	input[3](0) = 1.0 + 1.0 / sqrt(2.0); input[3](1) = 1.0;
	Triangulation triangulation = DelaunayLifting::Triangulate(input);
	//triangulation.PrintTriangulation();
	if (!triangulation.TestIntegrity())
		return false;
	if (!triangulation.TestDelaunay())
		return false;
	QuadEdge* mesh = triangulation.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	if (mesh->NumVertices() != 4 || mesh->NumEdges() != 5 || mesh->NumFaces() != 3)
		return false;
	EdgeIterator ite(mesh);
	Edge* e = ite.Next();
	bool found_e01 = false;
	bool found_e02 = false;
	bool found_e12 = false;
	bool found_e13 = false;
	bool found_e23 = false;
	while (e)
	{
		int n0 = e->GetOrig()->index;
		int n1 = e->GetDest()->index;
		if ((n0 == 0 && n1 == 3) || (n0 == 3 && n1 == 0))
			return false;
		if (!found_e01)
			if ((n0 == 0 && n1 == 1) || (n0 == 1 && n1 == 0))
				found_e01 = true;
		if (!found_e02)
			if ((n0 == 0 && n1 == 2) || (n0 == 2 && n1 == 0))
				found_e02 = true;
		if (!found_e12)
			if ((n0 == 1 && n1 == 2) || (n0 == 2 && n1 == 1))
				found_e12 = true;
		if (!found_e13)
			if ((n0 == 1 && n1 == 3) || (n0 == 3 && n1 == 1))
				found_e13 = true;
		if (!found_e23)
			if ((n0 == 2 && n1 == 3) || (n0 == 3 && n1 == 2))
				found_e23 =  true;
		e = ite.Next();
	}
	if (!found_e01 || !found_e02 || !found_e12 || !found_e13 || !found_e23)
		return false;
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_3(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(4);
	input[1](0) = 1;
	input[2](1) = 1;
	input[3](0) = 1; input[3](1) = 1;
	Triangulation triangulation = DelaunayLifting::Triangulate(input);
	//triangulation.PrintTriangulation();
	if (!triangulation.TestIntegrity())
		return false;
	if (!triangulation.TestDelaunay())
		return false;
	QuadEdge* mesh = triangulation.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	if (mesh->NumVertices() != 4 || mesh->NumEdges() != 5 || mesh->NumFaces() != 3)
		return false;
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_4(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(8);
	input[1](0) = 1;
	input[2](0) = 2;
	input[3](0) = 1.5; input[3](1) = 1;
	input[4](0) = 2; input[4](1) = 2;
	input[5](0) = 1; input[5](1) = 2.5;
	input[6](1) = 2;
	input[7](1) = 1;
	Triangulation triangulation = DelaunayLifting::Triangulate(input);
	//triangulation.PrintTriangulation();
	if (!triangulation.TestIntegrity())
		return false;
	if (!triangulation.TestDelaunay())
		return false;
	QuadEdge* mesh = triangulation.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	if (mesh->NumVertices() != 8 || mesh->NumEdges() != 14 || mesh->NumFaces() != 8)
		return false;
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_5(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(16);
	input[0](0) = 0; input[0](1) = 0;
	input[1](0) = 1; input[1](1) = 0;
	input[2](0) = 2; input[2](1) = 0;
	input[3](0) = 3; input[3](1) = 0;
	input[4](0) = 0; input[4](1) = 1;
	input[5](0) = 1; input[5](1) = 1;
	input[6](0) = 2; input[6](1) = 1;
	input[7](0) = 3; input[7](1) = 1;
	input[8](0) = 0; input[8](1) = 2;
	input[9](0) = 1; input[9](1) = 2;
	input[10](0) = 2; input[10](1) = 2;
	input[11](0) = 3; input[11](1) = 2;
	input[12](0) = 0; input[12](1) = 3;
	input[13](0) = 1; input[13](1) = 3;
	input[14](0) = 2; input[14](1) = 3;
	input[15](0) = 3; input[15](1) = 3;
	Triangulation triangulation = DelaunayLifting::Triangulate(input);
	//triangulation.PrintTriangulation();
	if (!triangulation.TestIntegrity())
		return false;
	if (!triangulation.TestDelaunay())
		return false;
	QuadEdge* mesh = triangulation.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	if (mesh->NumVertices() != 16 || mesh->NumEdges() != 33 || mesh->NumFaces() != 19)
		return false;
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_6(int& NumTests)
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
	Triangulation triangulation = DelaunayLifting::Triangulate(input);
	//triangulation.PrintTriangulation();
	if (!triangulation.TestIntegrity())
		return false;
	if (!triangulation.TestDelaunay())
		return false;
	QuadEdge* mesh = triangulation.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	if (mesh->NumVertices() != 8 || mesh->NumEdges() != 13 || mesh->NumFaces() != 7)
		return false;
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_7(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(16);
	input[0](0) = 0; input[0](1) = 0;
	input[1](0) = 1.3; input[1](1) = 0;
	input[2](0) = 1.9; input[2](1) = 0;
	input[3](0) = 3; input[3](1) = 0;
	input[4](0) = 0; input[4](1) = 1.1;
	input[5](0) = 1.5; input[5](1) = 0.6;
	input[6](0) = 1.9; input[6](1) = 1.5;
	input[7](0) = 3; input[7](1) = 0.9;
	input[8](0) = 0; input[8](1) = 1.8;
	input[9](0) = 0.9; input[9](1) = 2.3;
	input[10](0) = 2; input[10](1) = 2;
	input[11](0) = 3; input[11](1) = 2.05;
	input[12](0) = 0; input[12](1) = 3;
	input[13](0) = 1.1; input[13](1) = 3;
	input[14](0) = 1.75; input[14](1) = 3;
	input[15](0) = 3; input[15](1) = 3;
	//16 points, 12 boundary points. This should result in T = 2 * 16 - 12 - 2 = 18 triangles or 19 faces.
	Triangulation triangulation = DelaunayLifting::Triangulate(input);
	//triangulation.PrintTriangulation();
	if (!triangulation.TestIntegrity())
		return false;
	if (!triangulation.TestDelaunay())
		return false;
	QuadEdge* mesh = triangulation.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	const int Nv = 16;
	const int Ne = 33;
	const int Nf = 19;
	if (mesh->NumVertices() != Nv || mesh->NumEdges() != Ne || mesh->NumFaces() != Nf)
		return false;
	Vertex* v[Nv];
	for (int i = 0; i < Nv; ++i) v[i] = 0;
	Edge* e[Ne];
	for (int i = 0; i < Ne; ++i) e[i] = 0;
	Face* f[Nf];
	for (int i = 0; i < Nf; ++i) f[i] = 0;
	VertexIterator itv(mesh);
	Vertex* vv = itv.Next();
	while (vv){
		for (int i = 0; i < Nv; ++i){
			if (vv->GetPoint() == input[i]){
				v[i] = vv;
				break;
			}
		}
		vv = itv.Next();
	}
	for (int i = 0; i < Nv; ++i)
	{
		if (!v[i])
			return false;
	}
	EdgeIterator ite(mesh);
	Edge* ee = ite.Next();
	while (ee){
		Vertex* O = ee->GetOrig();
		Vertex* D = ee->GetDest();
		for (int i = 0; i < Ne; ++i){
			int j = -1, k = -1;
			if (i == 0){ j = 0; k = 1;}
			else if (i == 0){ j = 0; k = 1; }
			else if (i == 1){ j = 1; k = 2; }
			else if (i == 2){ j = 0; k = 4; }
			else if (i == 3){ j = 2; k = 3; }
			else if (i == 4){ j = 4; k = 5; }
			else if (i == 5){ j = 2; k = 5; }
			else if (i == 6){ j = 1; k = 5; }
			else if (i == 7){ j = 0; k = 5; }
			else if (i == 8){ j = 6; k = 5; }
			else if (i == 9){ j = 3; k = 7; }
			else if (i == 10){ j = 2; k = 7; }
			else if (i == 11){ j = 5; k = 7; }
			else if (i == 12){ j = 7; k = 6; }
			else if (i == 13){ j = 8; k = 4; }
			else if (i == 14){ j = 6; k = 9; }
			else if (i == 15){ j = 5; k = 9; }
			else if (i == 16){ j = 4; k = 9; }
			else if (i == 17){ j = 9; k = 8; }
			else if (i == 18){ j = 6; k = 10; }
			else if (i == 19){ j = 10; k = 9; }
			else if (i == 20){ j = 11; k = 10; }
			else if (i == 21){ j = 7; k = 11; }
			else if (i == 22){ j = 6; k = 11; }
			else if (i == 23){ j = 9; k = 12; }
			else if (i == 24){ j = 8; k = 12; }
			else if (i == 25){ j = 9; k = 13; }
			else if (i == 26){ j = 13; k = 12; }
			else if (i == 27){ j = 10; k = 14; }
			else if (i == 28){ j = 9; k = 14; }
			else if (i == 29){ j = 14; k = 13; }
			else if (i == 30){ j = 11; k = 15; }
			else if (i == 31){ j = 10; k = 15; }
			else if (i == 32){ j = 14; k = 15; }
			if ((O == v[j] && D == v[k]) || (O == v[k] && D == v[j])){
				e[i] = ee;
				break;
			}
		}
		ee = ite.Next();
	}
	for (int i = 0; i < Ne; ++i){
		if (!e[i])
			return false;
	}
	//Next check faces
	Face* fb = 0;
	FaceIterator itf(mesh);
	Face* ff = itf.Next();
	Vertex* bv[12];
	while(ff){
		VertexIterator itv(ff);
		Vertex* v0 = itv.Next();
		Vertex* v1 = itv.Next();
		Vertex* v2 = itv.Next();
		if (itv.Next()){
			if (!fb){
				fb = ff;
			}
			else
				return false;
		}
		for (int i = 0; i < Nf; ++i)
		{
			unordered_set<Vertex*> vf;
			vf.rehash(3);
			if (i == 0)     { vf.insert(v[2]); vf.insert(v[1]); vf.insert(v[5]); }
			else if (i == 1) { vf.insert(v[1]); vf.insert(v[0]); vf.insert(v[5]); }
			else if (i == 2){ vf.insert(v[0]); vf.insert(v[4]); vf.insert(v[5]); }
			else if (i == 3){ vf.insert(v[3]); vf.insert(v[2]); vf.insert(v[7]); }
			else if (i == 4){ vf.insert(v[2]); vf.insert(v[5]); vf.insert(v[7]); }
			else if (i == 5){ vf.insert(v[5]); vf.insert(v[6]); vf.insert(v[7]); }
			else if (i == 6){ vf.insert(v[6]); vf.insert(v[5]); vf.insert(v[9]); }
			else if (i == 7){ vf.insert(v[5]); vf.insert(v[4]); vf.insert(v[9]); }
			else if (i == 8){ vf.insert(v[4]); vf.insert(v[8]); vf.insert(v[9]); }
			else if (i == 9){ vf.insert(v[6]); vf.insert(v[9]); vf.insert(v[10]); }
			else if (i == 10){ vf.insert(v[7]); vf.insert(v[6]); vf.insert(v[11]); }
			else if (i == 11){ vf.insert(v[6]); vf.insert(v[10]); vf.insert(v[11]); }
			else if (i == 12){ vf.insert(v[9]); vf.insert(v[8]); vf.insert(v[12]); }
			else if (i == 13){ vf.insert(v[9]); vf.insert(v[12]); vf.insert(v[13]); }
			else if (i == 14){ vf.insert(v[10]); vf.insert(v[9]); vf.insert(v[14]); }
			else if (i == 15){ vf.insert(v[9]); vf.insert(v[13]); vf.insert(v[14]); }
			else if (i == 16){ vf.insert(v[11]); vf.insert(v[10]); vf.insert(v[15]); }
			else if (i == 17){ vf.insert(v[10]); vf.insert(v[14]); vf.insert(v[15]); }
			else if (i == 18){
				vf.insert(v[14]); vf.insert(v[13]); vf.insert(v[12]);
				vf.insert(v[8]); vf.insert(v[4]); vf.insert(v[0]);
				vf.insert(v[1]); vf.insert(v[2]); vf.insert(v[3]);
				vf.insert(v[7]); vf.insert(v[11]); vf.insert(v[15]);
			}
			unordered_set<Vertex*>::iterator its0,its1,its2,itsb;
			its0 = vf.find(v0);
			its1 = vf.find(v1);
			its2 = vf.find(v2);
			if (its0 != vf.end() && its1 != vf.end() && its2 != vf.end()){
				f[i] = ff;
			}
			if (ff == fb && i == 18){
				VertexIterator itvb(ff);
				Vertex* vb = itvb.Next();
				while (vb){
					itsb = vf.find(vb);
					if (itsb == vf.end())
						return false;
					vb = itvb.Next();
				}
			}
		}
		ff = itf.Next();
	}
	for (int i = 0; i < 19; ++i){
		if (!f[i])
			return false;
	}
	NumTests += 1;
	return true;
}
bool tester_DelaunayLifting_8(int& NumTests)
{//constrained Delaunay Triangulation	
	vector<Vector3D> input;
	input.resize(7);
	input[0](0) = 0; input[0](1) = 0;
	input[1](0) = 2; input[1](1) = 0;
	input[2](0) = 2; input[2](1) = 2;
	input[3](0) = 0; input[3](1) = 2;
	input[4](0) = 2.0/3.0; input[4](1) = 1.0/3.0;
	input[5](0) = 4.0 / 3.0; input[5](1) = 1.0 / 3.0;
	input[6](0) = 2.0 / 3.0; input[6](1) = 1.0;
	//Boundary: P0-P1-P2-P3-P4
	vector<int> boundary = { 0,1,2,3,4 };
	Triangulation T = DelaunayLifting::Triangulate(input,boundary);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[20];
	Edge* e[20];
	Face* f[20];
	for (int i = 0; i < 20; ++i) v[i] = 0;
	for (int i = 0; i < 20; ++i) e[i] = 0;
	for (int i = 0; i < 20; ++i) f[i] = 0;
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
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 3)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[4] || v[2]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[6] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[10] || v[3]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[8] || v[5]->GetDegree() != 5)
		return false;
	if (v[6]->GetEdge() != e[13] || v[6]->GetDegree() != 4)
		return false;
	if (f[0]->GetEdge() != e[1] || f[0]->GetDegree() != 5)
		return false;
	if (f[1]->GetEdge() != e[2] || f[1]->GetDegree() != 3)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[1] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[4] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[11] || f[5]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[6] || f[7]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[0])
		return false;
	if (e[0]->GetNextOrig() != e[7] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[8] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[8] || e[0]->GetPrevLeft() != e[7] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[3] || e[1]->GetRight() != f[0])
		return false;
	if (e[1]->GetNextOrig() != e[8] || e[1]->GetNextDest() != e[4] || e[1]->GetNextLeft() != e[5] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[5] || e[1]->GetPrevLeft() != e[8] || e[1]->GetPrevRight() != e[4])
		return false;
	if (e[2]->GetOrig() != v[0] || e[2]->GetDest() != v[4] || e[2]->GetLeft() != f[0] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[0] || e[2]->GetNextDest() != e[6] || e[2]->GetNextLeft() != e[11] || e[2]->GetNextRight() != e[7])
		return false;
	if (e[2]->GetPrevOrig() != e[7] || e[2]->GetPrevDest() != e[12] || e[2]->GetPrevLeft() != e[0] || e[2]->GetPrevRight() != e[6])
		return false;
	if (e[4]->GetOrig() != v[2] || e[4]->GetDest() != v[3] || e[4]->GetLeft() != f[4] || e[4]->GetRight() != f[0])
		return false;
	if (e[4]->GetNextOrig() != e[9] || e[4]->GetNextDest() != e[10] || e[4]->GetNextLeft() != e[10] || e[4]->GetNextRight() != e[1])
		return false;
	if (e[4]->GetPrevOrig() != e[1] || e[4]->GetPrevDest() != e[11] || e[4]->GetPrevLeft() != e[9] || e[4]->GetPrevRight() != e[11])
		return false;
	if (e[5]->GetOrig() != v[5] || e[5]->GetDest() != v[2] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[3])
		return false;
	if (e[5]->GetNextOrig() != e[6] || e[5]->GetNextDest() != e[1] || e[5]->GetNextLeft() != e[9] || e[5]->GetNextRight() != e[8])
		return false;
	if (e[5]->GetPrevOrig() != e[13] || e[5]->GetPrevDest() != e[9] || e[5]->GetPrevLeft() != e[13] || e[5]->GetPrevRight() != e[1])
		return false;
	if (e[6]->GetOrig() != v[4] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[7] || e[6]->GetRight() != f[1])
		return false;
	if (e[6]->GetNextOrig() != e[11] || e[6]->GetNextDest() != e[7] || e[6]->GetNextLeft() != e[13] || e[6]->GetNextRight() != e[2])
		return false;
	if (e[6]->GetPrevOrig() != e[2] || e[6]->GetPrevDest() != e[5] || e[6]->GetPrevLeft() != e[12] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[0] || e[7]->GetDest() != v[5] || e[7]->GetLeft() != f[1] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[2] || e[7]->GetNextDest() != e[8] || e[7]->GetNextLeft() != e[6] || e[7]->GetNextRight() != e[0])
		return false;
	if (e[7]->GetPrevOrig() != e[0] || e[7]->GetPrevDest() != e[6] || e[7]->GetPrevLeft() != e[2] || e[7]->GetPrevRight() != e[8])
		return false;
	if (e[8]->GetOrig() != v[1] || e[8]->GetDest() != v[5] || e[8]->GetLeft() != f[2] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[0] || e[8]->GetNextDest() != e[13] || e[8]->GetNextLeft() != e[7] || e[8]->GetNextRight() != e[1])
		return false;
	if (e[8]->GetPrevOrig() != e[1] || e[8]->GetPrevDest() != e[7] || e[8]->GetPrevLeft() != e[0] || e[8]->GetPrevRight() != e[5])
		return false;
	if (e[9]->GetOrig() != v[2] || e[9]->GetDest() != v[6] || e[9]->GetLeft() != f[8] || e[9]->GetRight() != f[4])
		return false;
	if (e[9]->GetNextOrig() != e[5] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[13] || e[9]->GetNextRight() != e[4])
		return false;
	if (e[9]->GetPrevOrig() != e[4] || e[9]->GetPrevDest() != e[13] || e[9]->GetPrevLeft() != e[5] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[6] || e[10]->GetDest() != v[3] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[12] || e[10]->GetNextDest() != e[11] || e[10]->GetNextLeft() != e[11] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[4] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[4])
		return false;
	if (e[12]->GetOrig() != v[4] || e[12]->GetDest() != v[6] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[7])
		return false;
	if (e[12]->GetNextOrig() != e[2] || e[12]->GetNextDest() != e[13] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[6])
		return false;
	if (e[12]->GetPrevOrig() != e[11] || e[12]->GetPrevDest() != e[10] || e[12]->GetPrevLeft() != e[11] || e[12]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[6] || e[13]->GetLeft() != f[7] || e[13]->GetRight() != f[8])
		return false;
	if (e[13]->GetNextOrig() != e[5] || e[13]->GetNextDest() != e[9] || e[13]->GetNextLeft() != e[12] || e[13]->GetNextRight() != e[5])
		return false;
	if (e[13]->GetPrevOrig() != e[8] || e[13]->GetPrevDest() != e[12] || e[13]->GetPrevLeft() != e[6] || e[13]->GetPrevRight() != e[9])
		return false;
	if (e[11]->GetOrig() != v[3] || e[11]->GetDest() != v[4] || e[11]->GetLeft() != f[5] || e[11]->GetRight() != f[0])
		return false;
	if (e[11]->GetNextOrig() != e[4] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[12] || e[11]->GetNextRight() != e[4])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[6] || e[11]->GetPrevLeft() != e[10] || e[11]->GetPrevRight() != e[2])
		return false;
	++NumTests;
	return true;
}
bool tester_DelaunayLifting_9(int& NumTests)//--------ERRROR---------------------------
{
	int Nx = 7;
	int Ny = Nx;
	double Lx = 10;
	double Ly = Lx;
	vector<Vector3D> input;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	int counter = -1;
	for (int i = 2; i < Ny; ++i) {
		input[++counter](0) = 0; input[counter](1) = double(i) * hy;
	}
	for (int i = 0; i < Nx; ++i) {
		input[++counter](0) = hx * double(i); input[counter](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i) {
		input[++counter](0) = Lx; input[counter](1) = Ly - double(i) * hy;
	}
	for (int i = 0; i < Nx; ++i) {
		input[++counter](0) = Lx - double(i) * hx; input[counter](1) = 0;
	}
	for (int i = 0; i < 2; ++i) {
		input[++counter](0) = 0; input[counter](1) = double(i) * hy;
	}
	//input.resize(23);
	Triangulation T = DelaunayLifting::Triangulate(input); 
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	++NumTests;
	return true;
}
bool tester_DelaunayLifting_10(int& NumTests)
{
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
	Triangulation T = DelaunayLifting::Triangulate(input);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
