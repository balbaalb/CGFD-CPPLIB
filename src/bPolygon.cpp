#include "bPolygon.h"
#include "QuadEdge.h"
#include "EdgeIterator.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "Vertex.h"
#include "Edge.h"
#include "LineSegment2D.h"
#include "MathUtils.h"
void bPolygon::CopyBody(const bPolygon& rhs)
{
	this->qe = 0;
	this->inFace = 0;
	this->outFace = 0;
	if (rhs.qe){
		this->qe = new QuadEdge(*rhs.qe);
		if (rhs.inFace == rhs.qe->GetFace(0))
		{
			this->inFace = this->qe->GetFace(0);
			this->outFace = this->qe->GetFace(1);
		}
		else
		{
			this->inFace = this->qe->GetFace(1);
			this->outFace = this->qe->GetFace(0);
		}
	}
	this->closed = rhs.closed;
	this->latestVertex = 0;
	this->firstVertex = 0;
	if (this->closed)
	{
		this->UpdateArrows();
	}
}
void bPolygon::DeleteBody()
{
	delete this->qe;
	this->qe = 0;
	this->closed = false;
	this->latestVertex = 0;
	this->firstVertex = 0;
	this->inFace = 0;
	this->outFace = 0;
	this->Arrows.empty();
}
void bPolygon::UpdateArrows()
{
	if (this->closed){
		int Df = this->qe->UpdateFaceDegree(this->inFace);
		this->Arrows.empty();
		this->Arrows.rehash(Df);
		EdgeIterator ite(this->inFace);
		Edge* e = ite.Next();
		while (e){
			Vertex* v = e->GetDest();
			this->Arrows[v] = e;
			e = ite.Next();
		}
	}
}
bool bPolygon::isCCW() const
{
	double sumTheta = 0;
	Vertex* a = 0;
	Vertex* b = 0;
	Vertex* c = 0;
	Vertex* v0 = 0;
	Vertex* v1 = 0;
	VertexIterator itv(this->qe);
	Vertex* v = itv.Next();
	while (v){
		a = b;
		b = c;
		c = v;
		if (a && b && c)
			sumTheta += Vector3D::RotationAngle_z(a->GetPoint(), b->GetPoint(), c->GetPoint());
		if (!v0)
			v0 = v;
		else if (!v1)
			v1 = v;
		v = itv.Next();
	}
	a = b;
	b = c;
	c = v0;
	sumTheta += Vector3D::RotationAngle_z(a->GetPoint(), b->GetPoint(), c->GetPoint());
	a = b;
	b = c;
	c = v1;
	sumTheta += Vector3D::RotationAngle_z(a->GetPoint(), b->GetPoint(), c->GetPoint());
	return (sumTheta > 0);
}
void bPolygon::Reverse()
{
	if (this->closed)
	{
		EdgeIterator ite(this->qe);
		Edge* e = ite.Next();
		while (e){
			e->Reverse();
			e = ite.Next();
		}
		Vertex* v = this->firstVertex;
		this->firstVertex = this->latestVertex;
		this->latestVertex = v;
		this->qe->ReverseFace(this->inFace);
		this->qe->ReverseFace(this->outFace);
	}
}
void bPolygon::ReArrangeVertices()
{
	int Nv = this->qe->NumVertices();
	vector<Vertex*> vArray;
	vArray.resize(Nv, 0);
	VertexIterator itv(this->inFace);
	Vertex* v = itv.Next();
	int i = 0;
	while (v){
		if (i >= Nv)
			throw "bPolygon::ReArrangeVertices()";
		vArray[i] = v;
		++i;
		v = itv.Next();
	}
	if (i != Nv)
		throw "bPolygon::ReArrangeVertices()";
	for (int i = 0; i < Nv; ++i){
		int ip1 = (i != Nv - 1) ? i + 1 : 0;
		int im1 = i ? i - 1 : Nv - 1;
		vArray[i]->next = vArray[ip1];
		vArray[i]->prev = vArray[im1];
	}
}
bool bPolygon::isEar(Edge* e)//O(n)
{
	Vertex* v = e->GetDest();
	Vertex* vp = e->GetOrig();
	Vertex* vn = e->GetNextLeft()->GetDest();
	Vector3D A = vp->GetPoint();
	Vector3D B = v->GetPoint();
	Vector3D C = vn->GetPoint();
	double theta = Vector3D::RotationAngle_z(A, B, C);
	if (theta <= 0)
		return false;
	LineSegment2D AC(A, C);
	Vertex* vpp = e->GetPrevLeft()->GetOrig();
	Vector3D D = vpp->GetPoint();
	if (AC.GetRelationTo(D) != ON_LEFT)
		return false;
	EdgeIterator ite(this->inFace);
	Edge* ee = ite.Next();
	while (ee){
		Vertex* vO = ee->GetOrig();
		Vertex* vD = ee->GetDest();
		Vector3D pO = vO->GetPoint();
		Vector3D pD = vD->GetPoint();
		LineSegment2D LE(pO, pD);
		if (AC.isIntersecting(LE))
			return false;
		ee = ite.Next();
	}
	return true;
}
Edge* bPolygon::FindAnEar()//O(n^2)
{
	EdgeIterator ite(this->inFace);
	Edge* e = ite.Next();
	while (e){
		bool is_Ear = this->isEar(e);//O(n)
		if (is_Ear)
			return e;
		e = ite.Next();
	}
	return 0;
}
bPolygon::POINT_POLYGON_RELATION bPolygon::GetRelationTo(const Vector3D& p) const//O(n)
{
	EdgeIterator ite(this->inFace);
	Edge* e = ite.Next();
	int numIntersecting = 0;
	while (e)
	{
		Vertex* vO = e->GetOrig();
		Vertex* vD = e->GetDest();
		Vector3D pO = vO->GetPoint();
		Vector3D pD = vD->GetPoint();
		LineSegment2D LE(pO, pD);
		if (LE.isOn(p))
			return POINT_POLYGON_RELATION::ON_POLYGON;
		if (LE.isIntersectingRayFrom(p,1.356844))
		{
			numIntersecting += 1;
			if (fabs(p(1) - pO(1)) < 1.0e-15 || fabs(p(1) - pD(1)) < 1.0e-15)
			{
				Vector3D mOD = (pO + pD) / 2.0;
				if(mOD(1) > p(1) - 1.0e-15)
					numIntersecting -= 1;
			}
		}
		e = ite.Next();
	}
	if (numIntersecting % 2 == 1)
		return POINT_POLYGON_RELATION::INSIDE_POLYGON;
	return POINT_POLYGON_RELATION::OUTSIDE_POLYGON;
}
bPolygon::bPolygon()
{
	this->qe = 0;
	this->closed = false;
	this->latestVertex = 0;
	this->firstVertex = 0;
	this->inFace = 0;
	this->outFace = 0;
}
bPolygon::bPolygon(const bPolygon& rhs)
{
	this->CopyBody(rhs);
}
bPolygon::~bPolygon()
{
	this->DeleteBody();
}
bPolygon& bPolygon::operator=(const bPolygon& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
	return *this;
}
QuadEdge* bPolygon::GetQuadEdgePtr() const
{
	QuadEdge* qeCopy = 0;
	if(this->qe)
		qeCopy = new QuadEdge(*this->qe);
	return qeCopy;
}
void bPolygon::AddVertex(const Vector3D& p)
{
	if (this->closed)
		throw "bPolygon::AddVertex(): Addding vertex to a closed Polygon.";
	if (p(2))
		throw "bPolygon::AddVertex(): 3D points are not accepted.";
	if (!this->qe)
		this->qe = new QuadEdge;
	Vertex* v = this->qe->AddVertex(p);
	if (!this->firstVertex)
		this->firstVertex = v;
	if (this->latestVertex)
	{
		Edge* e = this->qe->AddEdge();
		this->qe->AddEdgeToVertex(v, e, true);
		this->qe->AddEdgeToVertex(this->latestVertex, e, false);
	}
	this->latestVertex = v;
}
void bPolygon::AddVerticesAndClose(const vector<Vector3D>& points)
{
	for (int i = 0; i < points.size(); ++i)
	{
		this->AddVertex(points[i]);
	}
	this->Close();
}
void bPolygon::Close()
{
	if (!this->qe || this->qe->NumVertices() < 3)
		throw "Polygon::Close(): Not Enough Vertices";
	if (!this->closed)
	{
		Edge* eLast = this->qe->AddEdge(); 
		this->qe->AddEdgeToVertex(this->latestVertex, eLast, false);
		this->qe->AddEdgeToVertex(this->firstVertex, eLast, true);
		this->closed = true;
		this->inFace = this->qe->AddFace();
		this->outFace = this->qe->AddFace();
		bool ccw = this->isCCW();
		EdgeIterator ite(this->qe);
		Edge* e = ite.Next();
		while (e){
			this->qe->AddEdgeToFace(this->inFace, e, ccw);//assuming ccw
			this->qe->AddEdgeToFace(this->outFace, e, !ccw);//assuming ccw
			e = ite.Next();
		}
		if (!ccw)
			this->Reverse();
		this->ReArrangeVertices();
		this->UpdateArrows();
	}
}
Face* bPolygon::GetInsideFace()
{
	return this->inFace;
}
bool bPolygon::isEar(Vertex* v)
{
	Edge* e = this->Arrows[v];
	return this->isEar(e);
}
Triangulation bPolygon::Triangulate() const//O(n^3)
{
	Triangulation T;
	if (this->closed){
		int Nv = this->qe->NumVertices();
		bPolygon P(*this);
		while (Nv > 3)
		{
			Edge* e = P.FindAnEar();
			if (!e)
				throw "bPolygon::Triangulate() : CANNOT FIND AN EAR!!";
			Vertex* v = e->GetDest();
			Vertex* vp = e->GetOrig();
			Vertex* vn = e->GetNextLeft()->GetDest();
			Edge* eNew = P.qe->DivideFaceOnRight(P.inFace, vp, vn);
			--Nv;
		}
		T << *P.qe;
	}
	return T;
}
bool bPolygon::isInside(const Vector3D& p) const//O(n)
{
	return (this->GetRelationTo(p) == POINT_POLYGON_RELATION::INSIDE_POLYGON);
}
bool bPolygon::isOn(const Vector3D& p) const//O(n)
{
	return (this->GetRelationTo(p) == POINT_POLYGON_RELATION::ON_POLYGON);
}
bool bPolygon::isOnOrInside(const Vector3D& p) const//O(n)
{
	POINT_POLYGON_RELATION relationToP = this->GetRelationTo(p);
	return (relationToP == POINT_POLYGON_RELATION::ON_POLYGON 
		|| relationToP == POINT_POLYGON_RELATION::INSIDE_POLYGON);
}
bool bPolygon::TestIntegrity()
{
	if (this->closed)
	{
		EdgeIterator ite(this->inFace);
		Edge* e = ite.Next();
		while (e)
		{
			LineSegment2D L(e->GetOrig()->GetPoint(), e->GetDest()->GetPoint());
			EdgeIterator ite2(this->inFace);
			Edge* e2 = ite2.Next();
			while (e2)
			{
				LineSegment2D L2(e2->GetOrig()->GetPoint(), e2->GetDest()->GetPoint()); 
				if (e2 != e && L.isIntersecting(L2))
					return false;
				e2 = ite2.Next();
			}
			e = ite.Next();
		}
		if (this->qe->NumVertices() > 3)
		{
			Edge* ear = this->FindAnEar();
			if (!ear)
				return false;
		}
	}
	return true;
}