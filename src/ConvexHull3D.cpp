#include "math.h"
#include <cmath>
#include "ConvexHull3D.h"
#include "QuadEdge.h"
#include "Vector3D.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "EdgeIterator.h"
#include "EdgeList.h"
FaceVisible::FaceVisible() : Face()
{
	this->type = 2;
}
bool FaceVisible::GetVisible() const
{
	return this->visible;
}
void FaceVisible::SetVisible(bool v)
{
	int t = this->type;
	this->visible = v;
}
Face* FaceVisible::Clone() const
{
	return (new FaceVisible);
}
EdgeVisible::EdgeVisible()
{
	this->type = 2;
	this->visibility = EDGE_STATUS_UNKNOWN;
}
int EdgeVisible::GetVisiblity() const
{
	return this->visibility;
}
void EdgeVisible::SetVisiblity(EDGE_VISIBILITY v)
{
	if (v < 0 || v > 4)
		throw "EdgeVisible::SetVisiblity()";
	this->visibility = v;
}
Edge* EdgeVisible::Clone() const
{
	return (new EdgeVisible);
}
void ConvexHull3D::CopyBody(const ConvexHull3D& rhs)
{
	this->hull = 0;
	if (rhs.hull)
		this->hull = new QuadEdge(*rhs.hull);
	this->inputCopy.resize(0);
	for (int i = 0; i < rhs.inputCopy.size(); ++i)
		this->inputCopy.push_back(rhs.inputCopy[i]);
}
void ConvexHull3D::DeleteBody()
{
	delete this->hull;
	this->hull = 0;
}
bool ConvexHull3D::RemoveEdge(EdgeContainer* E)
{
	if (!this->hull || !this->hull->isRemovable(E->edge))
		return false;
	this->hull->DeleteEdgeAndRightFace(E);
	return true;
}
bool ConvexHull3D::UpdateFacesVisibility(const Vector3D& p)
{
	bool oneVisible = false;
	if (!this->hull)
		throw "ConvexHull3D::UpdateFacesVisibility()";
	FaceIterator itf(this->hull);
	FaceVisible* f = (FaceVisible*) itf.Next();
	while (f)
	{
		POINT_FACE_VISIBILITY visiblity = this->hull->isVisible(f, p);
		bool v = (visiblity == COPLANAR || visiblity == VISIBLE);
		oneVisible = oneVisible || v;
		f->SetVisible(v);
		f = (FaceVisible*) itf.Next();
	}
	return oneVisible;
}
Edge* ConvexHull3D::FindOneBoundaryEdge() const
{
	EdgeIterator ite(this->hull);
	Edge* e = ite.Next();
	while (e)
	{
		FaceVisible* fR = (FaceVisible*)e->GetRight();
		FaceVisible* fL = (FaceVisible*)e->GetLeft();
		bool vR = fR->GetVisible();
		bool vL = fL->GetVisible();
		if (vR != vL)
			return e;
		e = ite.Next();
	}
	return 0;
}
void ConvexHull3D::UpdateEdgesVisibility()
{
	EdgeIterator ite(this->hull);
	EdgeVisible* e = (EdgeVisible*) ite.Next();
	while (e)
	{
		FaceVisible* fR = (FaceVisible*)e->GetRight();
		FaceVisible* fL = (FaceVisible*)e->GetLeft();
		bool vR = fR->GetVisible();
		bool vL = fL->GetVisible();
		if (!vR && !vL)
			e->SetVisiblity(BOTH_FACES_INVISIBLE);
		else if (vR && vL)
			e->SetVisiblity(BOTH_FACES_VISIBLE);
		else if (vR)
			e->SetVisiblity(RIGHT_FACE_VISIBLE);
		else
			e->SetVisiblity(LEFT_FACE_VISIBLE);
		e = (EdgeVisible*) ite.Next();
	}
}
void ConvexHull3D::DeleteEdges()
{
	EdgeIterator ite(this->hull);
	EdgeVisible* e = (EdgeVisible*)ite.Next();
	while (e)
	{
		EdgeContainer* E = ite.GetCurrentContainer();
		if (e->GetVisiblity() == BOTH_FACES_VISIBLE)
		{
			bool edgeRemoved = this->RemoveEdge(E);
			if (!edgeRemoved)
				this->hull->GetEdgeList()->MoveToEnd(E);
		}
		e = (EdgeVisible*)ite.Next();
	}
}
Vertex* ConvexHull3D::AddPyramid(const Vector3D& p, Edge* eBoundary)
{
	FaceVisible* fR = (FaceVisible*)eBoundary->GetRight();
	FaceVisible* fL = (FaceVisible*)eBoundary->GetLeft();
	Vertex* v = 0;
	if (fR->GetVisible())
		v = this->hull->AddPyramidAndDeleteFace(p, fR);
	else if (fL->GetVisible())
		v = this->hull->AddPyramidAndDeleteFace(p, fL);
	else
		throw "ConvexHull3D::AddPyramid()";
	return v;
}
bool ConvexHull3D::PointAcceptable(const Vector3D& p)
{
	if (!this->hull)
		return false;
	VertexIterator itv(this->hull);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D q = v->GetPoint();
		if (p == q)
			return false;
		v = itv.Next();
	}
	return true;
}
void ConvexHull3D::BuildTetrahedron(const vector<Vector3D>& input, vector<int>& usedPoints)
{
	int iA = -1, iB = -1, iC = -1, iD = -1;
	Vector3D A, B, C, D;
	Vector3D AB, AC, AD, n;
	for (int i = 0; i < input.size(); ++i)
	{
		if (iA == -1)
		{
			iA = i;
			A = input[i];
		}
		else if (iB == -1)
		{
			B = input[i];
			AB = B - A;
			if (AB.abs() > 1e-10)
				iB = i;
		}
		else if (iC == -1)
		{
			C = input[i];
			AC = C - A;
			n = AB && AC;
			if (n.abs() > 1e-10)
			{
				n = n / n.abs();
				iC = i;
			}
		}
		else if (iD == -1)
		{
			D = input[i];
			AD = D - A;
			if (fabs(AD || n) > 1e-10)
			{
				iD = i;
				this->hull = new QuadEdge;
				this->hull->SetPrototype(new FaceVisible);
				this->hull->SetPrototype(new EdgeVisible);
				this->hull->SetAsTetrahedron(A, B, C, D);
				break;
			}
		}
	}
	if (iA != -1 && iB != -1 && iC != -1 && iD == -1)
	{
		for (int i = 0; i < input.size(); ++i)
		{
			if (i != iA && i != iB && i != iC)
			{
				D = input[i];
				if (D != A && D != B && D != C)
				{
					iD = i;
					this->hull = new QuadEdge;
					this->hull->SetPrototype(new FaceVisible);
					this->hull->SetPrototype(new EdgeVisible);
					this->hull->SetAsDoubleTriangulatedSquare(A, B, C, D);
					break;
				}
			}
		}
	}
	if (iD != -1)
	{
		usedPoints.push_back(iA);
		usedPoints.push_back(iB);
		usedPoints.push_back(iC);
		usedPoints.push_back(iD);
		VertexIterator itv(this->hull);
		itv.Next()->index = iA;
		itv.Next()->index = iB;
		itv.Next()->index = iC;
		itv.Next()->index = iD;
	}
}
QuadEdge* ConvexHull3D::GetHullPointer()
{
	return this->hull;
}
ConvexHull3D::ConvexHull3D()
{
	this->hull = 0;
}
ConvexHull3D::ConvexHull3D(const ConvexHull3D& rhs)
{
	this->CopyBody(rhs);
}
ConvexHull3D::~ConvexHull3D()
{
	this->DeleteBody();
}
void ConvexHull3D::AddPoint(const Vector3D& p, int index)
{
	if (!this->hull)
	{
		this->inputCopy.push_back(p);
		if (this->inputCopy.size() >= 4)
		{
			this->Build(this->inputCopy);
		}
	}
	if (this->hull)
	{
		bool pointAcceptable = this->PointAcceptable(p);
		if (!pointAcceptable)
			return;
		bool oneVisible = this->UpdateFacesVisibility(p);
		if (!oneVisible)
			return;//point is inside the convex hull
		Edge* eBoundary = this->FindOneBoundaryEdge();
		if (!eBoundary)
			throw "ConvexHull3D: Boundary edge is Null.";
		this->UpdateEdgesVisibility();
		this->DeleteEdges();
		this->hull->SetNextVetexIndex(index);
		this->AddPyramid(p, eBoundary);
	}
}
ConvexHull3D& ConvexHull3D::operator=(const ConvexHull3D& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
	return *this;
}
void ConvexHull3D::Build(const vector<Vector3D>& input)
{
	vector<int> usedPoints;
	this->BuildTetrahedron(input, usedPoints);
	if (this->hull)
	{
		for (int i = 0; i < input.size(); ++i)
		{
			if (i != usedPoints[0] && i != usedPoints[1] && i != usedPoints[2] && i != usedPoints[3])
			{
				this->AddPoint(input[i],i);
			}
		}
	}
}
QuadEdge* ConvexHull3D::GetHull() const
{
	QuadEdge* QE = 0;
	if (this->hull)
		QE = new QuadEdge(*this->hull);
	return QE;
}
bool ConvexHull3D::GetFaceVisibiliy(int i) const
{
	FaceVisible* f = (FaceVisible*) this->hull->GetFace(i);
	return f->GetVisible();
}
int ConvexHull3D::GetEdgeVisibiliy(int i) const
{
	EdgeVisible* e = (EdgeVisible*) this->hull->GetEdge(i);
	return e->GetVisiblity();
}
bool ConvexHull3D::TestConvexity() const
{
	if (this->hull)
	{
		VertexIterator itv(this->hull);
		Vertex* v = itv.Next();
		while (v)
		{
			vector<Face*> adjacentFaces;
			FaceIterator itf0(v);
			Face* f = itf0.Next();
			while (f)
			{
				adjacentFaces.push_back(f);
				f = itf0.Next();
			}
			FaceIterator itf(this->hull);
			f = itf.Next();
			Vector3D p = v->GetPoint();
			while (f)
			{
				bool ignore = false;
				for (int i = 0; i < adjacentFaces.size(); ++i)
				{
					if (adjacentFaces[i] == f)
					{
						ignore = true;
						break;
					}
				}
				if (!ignore)
				{
					POINT_FACE_VISIBILITY visible = this->hull->isVisible(f, p);
					if (visible == VISIBLE)
						return false;
				}
				f = itf.Next();
			}
			v = itv.Next();
		}
	}
	return true;
}