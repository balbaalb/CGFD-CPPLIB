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