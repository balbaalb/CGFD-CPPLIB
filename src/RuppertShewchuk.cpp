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
