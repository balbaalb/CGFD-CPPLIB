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
Vertex* ConvexHull3D::AddPoint(const Vector3D& p, int index)
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
			return 0;
		bool oneVisible = this->UpdateFacesVisibility(p);
		if (!oneVisible)
			return 0;//point is inside the convex hull
		Edge* eBoundary = this->FindOneBoundaryEdge();
		if (!eBoundary)
			throw "ConvexHull3D: Boundary edge is Null.";
		this->UpdateEdgesVisibility();
		this->DeleteEdges();
		this->hull->SetNextVetexIndex(index);
		Vertex* v = this->AddPyramid(p, eBoundary);
		return v;
	}
	return nullptr;
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
				Vertex* v = this->AddPoint(input[i],i);
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
bool tester_ConvexHull3D(int& NumTests)
{
	//Tetrahedron and a point inside it
	vector<Vector3D> input;
	input.resize(5);
	input[1](0) = 1;
	input[2](1) = 1;
	input[3](2) = 1;
	input[4](0) = input[4](1) = input[4](2) = 1.0 / 4.0;
	ConvexHull3D h1;
	h1.Build(input);
	QuadEdge* h1qe = h1.GetHull();
	if (!h1qe->TestIntegrity())
		return false;
	if (h1qe->NumVertices() != 4 || h1qe->NumFaces() != 4 || h1qe->NumEdges() != 6)
		return false;
	if (!h1.TestConvexity())
		return false;
	input[4](0) = input[4](1) = input[4](2) = 1.0 / 3.0;
	h1.Build(input);
	h1qe = h1.GetHull();
	if (!h1qe->TestIntegrity())
		return false;
	if (h1qe->NumVertices() != 5 || h1qe->NumFaces() != 6 || h1qe->NumEdges() != 9)
		return false;
	if (!h1.TestConvexity())
		return false;
	if (!tester_ConvexHull3D_1(NumTests))
		return false;
	if (!tester_ConvexHull3D_2(NumTests))
		return false;
	if (!tester_ConvexHull3D_PyramidToCube(NumTests))
		return false;
	if (!tester_ConvexHull3D_3(NumTests))
		return false;
	if (!tester_ConvexHull3D_Icosahedron(NumTests))
		return false;
	if (!tester_ConvexHull3D_Dodecahedron(NumTests))
		return false;
	NumTests += 1;
	return true;
}
bool tester_ConvexHull3D_1(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(5);
	input[1](0) = 1;
	input[2](1) = 1;
	input[3](0) = input[3](1) = input[3](2) = 1.0 / 3.0;
	input[4](2) = 1;
	ConvexHull3D h1;
	vector<int> usedPoints;
	h1.BuildTetrahedron(input, usedPoints);
	if (!h1.TestConvexity())
		return false;
	if (usedPoints.size() != 4 || usedPoints[0] != 0 || usedPoints[1] != 1 || usedPoints[2] != 2 || usedPoints[3] != 3)
		return false;
	QuadEdge* hull = h1.GetHull();
	if (!hull)
		return false;
	if (!hull->TestIntegrity())
		return false;
	Vertex* v[4];
	Face* f[4];
	Edge* e[6];
	EdgeContainer* E[6];
	if (!tester_QuadEdge_Tetrahedron(*hull, v, f, e, input))
		return false;
	EdgeVisible* ev[6];
	for (int i = 0; i < 6; ++i)
		ev[i] = (EdgeVisible*)e[i];
	for (int i = 0; i < 6; ++i)
		E[i] = hull->GetEdgeList()->GetContainer(i);
	bool oneVisible = h1.UpdateFacesVisibility(input[4]);
	if (!oneVisible)
		return false;
	if (h1.GetFaceVisibiliy(0) || !h1.GetFaceVisibiliy(1) || !h1.GetFaceVisibiliy(2) || !h1.GetFaceVisibiliy(3))
		return false;
	Edge* eBoundary = h1.FindOneBoundaryEdge();
	if (!eBoundary)
		return false;
	h1.UpdateEdgesVisibility();
	int n = ev[0]->GetVisiblity();
	if (h1.GetEdgeVisibiliy(0) != LEFT_FACE_VISIBLE || h1.GetEdgeVisibiliy(1) != LEFT_FACE_VISIBLE
		|| h1.GetEdgeVisibiliy(2) != LEFT_FACE_VISIBLE)
		return false;
	if (h1.GetEdgeVisibiliy(3) != BOTH_FACES_VISIBLE || h1.GetEdgeVisibiliy(4) != BOTH_FACES_VISIBLE
		|| h1.GetEdgeVisibiliy(5) != BOTH_FACES_VISIBLE)
		return false;
	h1.DeleteEdges();
	hull = h1.GetHull();
	if (!hull)
		return false;
	if (!hull->TestIntegrity())
		return false;
	if (!hull->TestIntegrity())
		return false;
	if (hull->NumVertices() != 3 || hull->NumFaces() != 2 || hull->NumEdges() != 3)
		return false;
	for (int i = 0; i < 3; ++i)
	{
		v[i] = hull->GetVertex(i);
		e[i] = hull->GetEdge(i);
		f[i] = i < 2 ? hull->GetFace(i) : 0;
	}

	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2])
		return false;
	if (f[0]->GetEdge() != e[0] || f[1]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[1])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[1])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[2] || e[0]->GetPrevOrig() != e[2] || e[0]->GetNextDest() != e[1] || e[0]->GetPrevDest() != e[1])
		return false;
	if (e[1]->GetNextOrig() != e[0] || e[1]->GetPrevOrig() != e[0] || e[1]->GetNextDest() != e[2] || e[1]->GetPrevDest() != e[2])
		return false;
	if (e[2]->GetNextOrig() != e[1] || e[2]->GetPrevOrig() != e[1] || e[2]->GetNextDest() != e[0] || e[2]->GetPrevDest() != e[0])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[2] || e[0]->GetPrevRight() != e[1] || e[0]->GetNextLeft() != e[1] || e[0]->GetPrevLeft() != e[2])
		return false;
	if (e[1]->GetNextRight() != e[0] || e[1]->GetPrevRight() != e[2] || e[1]->GetNextLeft() != e[2] || e[1]->GetPrevLeft() != e[0])
		return false;
	if (e[2]->GetNextRight() != e[1] || e[2]->GetPrevRight() != e[0] || e[2]->GetNextLeft() != e[0] || e[2]->GetPrevLeft() != e[1])
		return false;
	NumTests += 1;
	return true;
}
bool tester_ConvexHull3D_2(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(5);
	input[1](0) = 1;
	input[2](1) = 1;
	input[3](0) = input[3](1) = input[3](2) = 1.0 / 3.0;
	input[4](2) = 1;
	ConvexHull3D h1;
	h1.Build(input);
	if (!h1.TestConvexity())
		return false;
	QuadEdge* hull = h1.GetHull();
	if (!hull)
		return false;
	if (!hull->TestIntegrity())
		return false;
	Vector3D p[4];
	for (int i = 0; i < 4; ++i)
		p[i] = hull->GetVertex(i)->GetPoint();
	if (p[0] != input[0] || p[1] != input[1] || p[2] != input[2] || p[3] != input[4])
		return false;
	Vertex* v[4];
	Face* f[4];
	Edge* e[6];
	EdgeContainer* E[6];
	if (hull->NumVertices() != 4 || hull->NumEdges() != 6 || hull->NumFaces() != 4)
		return false;
	for (int i = 0; i < 6; ++i)
	{
		e[i] = hull->GetEdge(i);
		E[i] = hull->GetEdgeList()->GetContainer(i);
	}
	for (int i = 0; i < 4; ++i)
		f[i] = hull->GetFace(i);
	for (int i = 0; i < 4; ++i)
		v[i] = hull->GetVertex(i);
	for (int i = 0; i < 6; ++i)
	{
		if (E[i]->edge != e[i])
			return false;
	}
	vector<Edge*> edgeVerts[4];
	for (int i = 0; i < 4; ++i)
	{
		if (v[i]->GetDegree() != 3)
			return false;
	}
	hull->UpdateVerticesDegree();
	for (int i = 0; i < 4; ++i)
	{
		if (v[i]->GetDegree() != 3)
			return false;
	}
	for (int i = 0; i < 4; ++i)
	{
		hull->GetVertexEdges(v[i], edgeVerts[i]);
		if (edgeVerts[i].size() != 3)
			return false;
		if (v[i]->GetDegree() != 3)
			return false;
		if (hull->UpdateVertexDegree(v[i]) != 3)
			return false;
	}
	vector<Edge*> edgeFaces[4];
	for (int i = 0; i < 4; ++i)
	{
		hull->GetFaceEdges(hull->GetFace(i), edgeFaces[i]);
		if (edgeFaces[i].size() != 3)
			return false;
	}
	for (int i = 0; i < 6; ++i)
	{
		if (!hull->isRemovable(e[i]))
			return false;
	}
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2] || v[3]->GetEdge() != e[5])
		return false;
	if (f[0]->GetEdge() != e[0] || f[1]->GetEdge() != e[2] || f[2]->GetEdge() != e[0] || f[3]->GetEdge() != e[1])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[2])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[3])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[1])
		return false;
	if (e[3]->GetOrig() != v[2] || e[3]->GetDest() != v[3] || e[3]->GetRight() != f[1] || e[3]->GetLeft() != f[3])
		return false;
	if (e[4]->GetOrig() != v[0] || e[4]->GetDest() != v[3] || e[4]->GetRight() != f[2] || e[4]->GetLeft() != f[1])
		return false;
	if (e[5]->GetOrig() != v[1] || e[5]->GetDest() != v[3] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[4] || e[0]->GetPrevOrig() != e[2] || e[0]->GetNextDest() != e[1] || e[0]->GetPrevDest() != e[5])
		return false;	
	if (e[1]->GetNextOrig() != e[5] || e[1]->GetPrevOrig() != e[0] || e[1]->GetNextDest() != e[2] || e[1]->GetPrevDest() != e[3])
		return false;
	if (e[2]->GetNextOrig() != e[3] || e[2]->GetPrevOrig() != e[1] || e[2]->GetNextDest() != e[0] || e[2]->GetPrevDest() != e[4])
		return false;
	if (e[3]->GetNextOrig() != e[1] || e[3]->GetPrevOrig() != e[2] || e[3]->GetNextDest() != e[4] || e[3]->GetPrevDest() != e[5])
		return false;
	if (e[4]->GetNextOrig() != e[2] || e[4]->GetPrevOrig() != e[0] || e[4]->GetNextDest() != e[5] || e[4]->GetPrevDest() != e[3])
		return false;
	if (e[5]->GetNextOrig() != e[0] || e[5]->GetPrevOrig() != e[1] || e[5]->GetNextDest() != e[3] || e[5]->GetPrevDest() != e[4])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[2] || e[0]->GetPrevRight() != e[1] || e[0]->GetNextLeft() != e[5] || e[0]->GetPrevLeft() != e[4])
		return false;
	if (e[1]->GetNextRight() != e[0] || e[1]->GetPrevRight() != e[2] || e[1]->GetNextLeft() != e[3] || e[1]->GetPrevLeft() != e[5])
		return false;
	if (e[2]->GetNextRight() != e[1] || e[2]->GetPrevRight() != e[0] || e[2]->GetNextLeft() != e[4] || e[2]->GetPrevLeft() != e[3])
		return false;
	if (e[3]->GetNextRight() != e[2] || e[3]->GetPrevRight() != e[4] || e[3]->GetNextLeft() != e[5] || e[3]->GetPrevLeft() != e[1])
		return false;
	if (e[4]->GetNextRight() != e[0] || e[4]->GetPrevRight() != e[5] || e[4]->GetNextLeft() != e[3] || e[4]->GetPrevLeft() != e[2])
		return false;
	if (e[5]->GetNextRight() != e[1] || e[5]->GetPrevRight() != e[3] || e[5]->GetNextLeft() != e[4] || e[5]->GetPrevLeft() != e[0])
		return false;
	Vector3D Pv[4],Center;
	for (int i = 0; i < 4; ++i)
	{
		Pv[i] = v[i]->GetPoint();
		Center = Center + Pv[i] / 4.0;
	}
	if (hull->isVisible(f[0], Pv[3]))
		return false;
	if (hull->isVisible(f[1], Pv[1]))
		return false;
	if (hull->isVisible(f[2], Pv[2]))
		return false;
	if (hull->isVisible(f[3], Pv[0]))
		return false;
	if (hull->isVisible(f[0], Center))
		return false;
	if (hull->isVisible(f[1], Center))
		return false;
	if (hull->isVisible(f[2], Center))
		return false;
	if (hull->isVisible(f[3], Center))
		return false;
	if (!hull->TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
bool tester_ConvexHull3D_PyramidToCube(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(5);
	input[1](0) = 1;
	input[2](1) = 1;
	input[3](0) = input[3](1) = 1;
	input[4](0) = input[4](1) = input[4](2) = 0.4;
	ConvexHull3D h1;
	h1.Build(input);
	if (!h1.TestConvexity())
		return false;
	QuadEdge* hull = h1.GetHull();
	if (!hull)
		return false;
	if (!hull->TestIntegrity())
		return false;
	if (hull->NumVertices() != 5 || hull->NumEdges() != 9 || hull->NumFaces() != 6)
		return false;
	Vertex* v[8];
	Face* f[12];
	Edge* e[18];
	EdgeContainer* E[18];
	for (int i = 0; i < 5; ++i)
	{
		v[i] = hull->GetVertex(i);
	}
	for (int i = 0; i < 6; ++i)
	{
		f[i] = hull->GetFace(i);
	}
	for (int i = 0; i < 9; ++i)
	{
		e[i] = hull->GetEdge(i);
		E[i] = hull->GetEdgeList()->GetContainer(i);
	}
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[3] || v[2]->GetEdge() != e[1]
		|| v[3]->GetEdge() != e[4] || v[4]->GetEdge() != e[8])
		return false;
	if (f[0]->GetEdge() != e[0] || f[1]->GetEdge() != e[1] || f[2]->GetEdge() != e[4]
		|| f[3]->GetEdge() != e[3] || f[4]->GetEdge() != e[0] || f[5]->GetEdge() != e[1])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[1] || e[0]->GetDest() != v[0] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[4])
		return false;
	if (e[1]->GetOrig() != v[0] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[1] || e[1]->GetLeft() != f[5])
		return false;
	if (e[2]->GetOrig() != v[0] || e[2]->GetDest() != v[3] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[1])
		return false;
	if (e[3]->GetOrig() != v[3] || e[3]->GetDest() != v[1] || e[3]->GetRight() != f[0] || e[3]->GetLeft() != f[3])
		return false;
	if (e[4]->GetOrig() != v[2] || e[4]->GetDest() != v[3] || e[4]->GetRight() != f[1] || e[4]->GetLeft() != f[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[4] || e[5]->GetRight() != f[2] || e[5]->GetLeft() != f[5])
		return false;
	if (e[6]->GetOrig() != v[3] || e[6]->GetDest() != v[4] || e[6]->GetRight() != f[3] || e[6]->GetLeft() != f[2])
		return false;
	if (e[7]->GetOrig() != v[1] || e[7]->GetDest() != v[4] || e[7]->GetRight() != f[4] || e[7]->GetLeft() != f[3])
		return false;
	if (e[8]->GetOrig() != v[0] || e[8]->GetDest() != v[4] || e[8]->GetRight() != f[5] || e[8]->GetLeft() != f[4])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[3] || e[0]->GetPrevOrig() != e[7] || e[0]->GetNextDest() != e[8] || e[0]->GetPrevDest() != e[1])
		return false;
	if (e[1]->GetNextOrig() != e[0] || e[1]->GetPrevOrig() != e[2] || e[1]->GetNextDest() != e[5] || e[1]->GetPrevDest() != e[4])
		return false;
	if (e[2]->GetNextOrig() != e[1] || e[2]->GetPrevOrig() != e[8] || e[2]->GetNextDest() != e[3] || e[2]->GetPrevDest() != e[6])
		return false;
	if (e[3]->GetNextOrig() != e[4] || e[3]->GetPrevOrig() != e[2] || e[3]->GetNextDest() != e[7] || e[3]->GetPrevDest() != e[0])
		return false;
	if (e[4]->GetNextOrig() != e[1] || e[4]->GetPrevOrig() != e[5] || e[4]->GetNextDest() != e[6] || e[4]->GetPrevDest() != e[3])
		return false;
	if (e[5]->GetNextOrig() != e[4] || e[5]->GetPrevOrig() != e[1] || e[5]->GetNextDest() != e[6] || e[5]->GetPrevDest() != e[8])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetPrevOrig() != e[4] || e[6]->GetNextDest() != e[7] || e[6]->GetPrevDest() != e[5])
		return false;
	if (e[7]->GetNextOrig() != e[0] || e[7]->GetPrevOrig() != e[3] || e[7]->GetNextDest() != e[8] || e[7]->GetPrevDest() != e[6])
		return false;
	if (e[8]->GetNextOrig() != e[2] || e[8]->GetPrevOrig() != e[0] || e[8]->GetNextDest() != e[5] || e[8]->GetPrevDest() != e[7])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[3] || e[0]->GetPrevRight() != e[2] || e[0]->GetNextLeft() != e[8] || e[0]->GetPrevLeft() != e[7])
		return false;
	if (e[1]->GetNextRight() != e[2] || e[1]->GetPrevRight() != e[4] || e[1]->GetNextLeft() != e[5] || e[1]->GetPrevLeft() != e[8])
		return false;
	if (e[2]->GetNextRight() != e[0] || e[2]->GetPrevRight() != e[3] || e[2]->GetNextLeft() != e[4] || e[2]->GetPrevLeft() != e[1])
		return false;
	if (e[3]->GetNextRight() != e[2] || e[3]->GetPrevRight() != e[0] || e[3]->GetNextLeft() != e[7] || e[3]->GetPrevLeft() != e[6])
		return false;
	if (e[4]->GetNextRight() != e[1] || e[4]->GetPrevRight() != e[2] || e[4]->GetNextLeft() != e[6] || e[4]->GetPrevLeft() != e[5])
		return false;
	if (e[5]->GetNextRight() != e[4] || e[5]->GetPrevRight() != e[6] || e[5]->GetNextLeft() != e[8] || e[5]->GetPrevLeft() != e[1])
		return false;
	if (e[6]->GetNextRight() != e[3] || e[6]->GetPrevRight() != e[7] || e[6]->GetNextLeft() != e[5] || e[6]->GetPrevLeft() != e[4])
		return false;
	if (e[7]->GetNextRight() != e[0] || e[7]->GetPrevRight() != e[8] || e[7]->GetNextLeft() != e[6] || e[7]->GetPrevLeft() != e[3])
		return false;
	if (e[8]->GetNextRight() != e[1] || e[8]->GetPrevRight() != e[5] || e[8]->GetNextLeft() != e[7] || e[8]->GetPrevLeft() != e[0])
		return false;
	//---------------------------------------------------------------------------------
	for (int i = 0; i < 9; ++i)
	{
		if (E[i]->edge != e[i])
			return false;
	}
	//====================================================================================
	Vector3D dz(0, 0, 1);
	Vector3D p5 = input[0] + dz;
	h1.AddPoint(p5);
	hull = h1.GetHull();
	if (!hull->TestIntegrity())
		return false;
	if (!h1.TestConvexity())
		return false;
	if (hull->NumVertices() != 5 || hull->NumEdges() != 9 || hull->NumFaces() != 6)
		return false;
	Vector3D p6 = input[1] + dz;
	h1.AddPoint(p6);
	hull = h1.GetHull();
	if (!hull->TestIntegrity())
		return false;
	if (!h1.TestConvexity())
		return false;
	if (hull->NumVertices() != 6 || hull->NumEdges() != 12 || hull->NumFaces() != 8)
		return false;
	Vector3D p7 = input[2] + dz;
	h1.AddPoint(p7);
	hull = h1.GetHull();
	if (!hull->TestIntegrity())
		return false;
	if (!h1.TestConvexity())
		return false;
	if (hull->NumVertices() != 7 || hull->NumEdges() != 15 || hull->NumFaces() != 10)
		return false;
	Vector3D p8 = input[3] + dz;
	h1.AddPoint(p8);
	hull = h1.GetHull();
	if (!hull->TestIntegrity())
		return false;
	if (!h1.TestConvexity())
		return false;
	if (hull->NumVertices() != 8 || hull->NumEdges() != 18 || hull->NumFaces() != 12)
		return false;
	for (int i = 0; i < 8; ++i)
	{
		v[i] = hull->GetVertex(i);
	}
	for (int i = 0; i < 12; ++i)
	{
		f[i] = hull->GetFace(i);
	}
	for (int i = 0; i < 18; ++i)
	{
		e[i] = hull->GetEdge(i);
		E[i] = hull->GetEdgeList()->GetContainer(i);
	}
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[3] || v[2]->GetEdge() != e[1] || v[3]->GetEdge() != e[4] 
		|| v[4]->GetEdge() != e[8] || v[5]->GetEdge() != e[8] || v[6]->GetEdge() != e[9] || v[7]->GetEdge() != e[17])
		return false;
	if (f[0]->GetEdge() != e[0] || f[1]->GetEdge() != e[1] || f[2]->GetEdge() != e[0] || f[3]->GetEdge() != e[5] 
		|| f[4]->GetEdge() != e[1] || f[5]->GetEdge() != e[5] || f[6]->GetEdge() != e[8] || f[7]->GetEdge() != e[7]
		|| f[8]->GetEdge() != e[3] || f[9]->GetEdge() != e[2] || f[10]->GetEdge() != e[9] || f[11]->GetEdge() != e[11])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[2])
		return false;
	if (e[1]->GetOrig() != v[2] || e[1]->GetDest() != v[0] || e[1]->GetRight() != f[1] || e[1]->GetLeft() != f[4])
		return false;
	if (e[2]->GetOrig() != v[3] || e[2]->GetDest() != v[2] || e[2]->GetRight() != f[1] || e[2]->GetLeft() != f[9])
		return false;
	if (e[3]->GetOrig() != v[1] || e[3]->GetDest() != v[3] || e[3]->GetRight() != f[0] || e[3]->GetLeft() != f[8])
		return false;
	if (e[4]->GetOrig() != v[0] || e[4]->GetDest() != v[3] || e[4]->GetRight() != f[1] || e[4]->GetLeft() != f[0])
		return false;
	if (e[5]->GetOrig() != v[0] || e[5]->GetDest() != v[4] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[5])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetRight() != f[2] || e[6]->GetLeft() != f[3])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetRight() != f[2] || e[7]->GetLeft() != f[7])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[5] || e[8]->GetRight() != f[3] || e[8]->GetLeft() != f[6])
		return false;
	if (e[9]->GetOrig() != v[2] || e[9]->GetDest() != v[6] || e[9]->GetRight() != f[4] || e[9]->GetLeft() != f[10])
		return false;
	if (e[10]->GetOrig() != v[0] || e[10]->GetDest() != v[6] || e[10]->GetRight() != f[5] || e[10]->GetLeft() != f[4])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[4] || e[11]->GetRight() != f[5] || e[11]->GetLeft() != f[11])
		return false;
	if (e[12]->GetOrig() != v[4] || e[12]->GetDest() != v[7] || e[12]->GetRight() != f[6] || e[12]->GetLeft() != f[11])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[7] || e[13]->GetRight() != f[7] || e[13]->GetLeft() != f[6])
		return false;
	if (e[14]->GetOrig() != v[1] || e[14]->GetDest() != v[7] || e[14]->GetRight() != f[8] || e[14]->GetLeft() != f[7])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetRight() != f[9] || e[15]->GetLeft() != f[8])
		return false;
	if (e[16]->GetOrig() != v[2] || e[16]->GetDest() != v[7] || e[16]->GetRight() != f[10] || e[16]->GetLeft() != f[9])
		return false;
	if (e[17]->GetOrig() != v[6] || e[17]->GetDest() != v[7] || e[17]->GetRight() != f[11] || e[17]->GetLeft() != f[10])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetPrevOrig() != e[1] || e[0]->GetNextDest() != e[3] || e[0]->GetPrevDest() != e[7])
		return false;
	if (e[1]->GetNextOrig() != e[16] || e[1]->GetPrevOrig() != e[2] || e[1]->GetNextDest() != e[0] || e[1]->GetPrevDest() != e[4])
		return false;
	if (e[2]->GetNextOrig() != e[3] || e[2]->GetPrevOrig() != e[15] || e[2]->GetNextDest() != e[1] || e[2]->GetPrevDest() != e[9])
		return false;
	if (e[3]->GetNextOrig() != e[14] || e[3]->GetPrevOrig() != e[0] || e[3]->GetNextDest() != e[4] || e[3]->GetPrevDest() != e[2])
		return false;
	if (e[4]->GetNextOrig() != e[1] || e[4]->GetPrevOrig() != e[5] || e[4]->GetNextDest() != e[15] || e[4]->GetPrevDest() != e[3])
		return false;
	if (e[5]->GetNextOrig() != e[4] || e[5]->GetPrevOrig() != e[6] || e[5]->GetNextDest() != e[8] || e[5]->GetPrevDest() != e[11])
		return false;
	if (e[6]->GetNextOrig() != e[5] || e[6]->GetPrevOrig() != e[10] || e[6]->GetNextDest() != e[7] || e[6]->GetPrevDest() != e[13])
		return false;
	if (e[7]->GetNextOrig() != e[8] || e[7]->GetPrevOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetPrevDest() != e[14])
		return false;
	if (e[8]->GetNextOrig() != e[12] || e[8]->GetPrevOrig() != e[5] || e[8]->GetNextDest() != e[13] || e[8]->GetPrevDest() != e[7])
		return false;
	if (e[9]->GetNextOrig() != e[2] || e[9]->GetPrevOrig() != e[16] || e[9]->GetNextDest() != e[17] || e[9]->GetPrevDest() != e[11])
		return false;
	if (e[10]->GetNextOrig() != e[6] || e[10]->GetPrevOrig() != e[0] || e[10]->GetNextDest() != e[11] || e[10]->GetPrevDest() != e[17])
		return false;
	if (e[11]->GetNextOrig() != e[9] || e[11]->GetPrevOrig() != e[10] || e[11]->GetNextDest() != e[5] || e[11]->GetPrevDest() != e[12])
		return false;
	if (e[12]->GetNextOrig() != e[11] || e[12]->GetPrevOrig() != e[8] || e[12]->GetNextDest() != e[13] || e[12]->GetPrevDest() != e[17])
		return false;
	if (e[13]->GetNextOrig() != e[6] || e[13]->GetPrevOrig() != e[8] || e[13]->GetNextDest() != e[14] || e[13]->GetPrevDest() != e[12])
		return false;
	if (e[14]->GetNextOrig() != e[7] || e[14]->GetPrevOrig() != e[3] || e[14]->GetNextDest() != e[15] || e[14]->GetPrevDest() != e[13])
		return false;
	if (e[15]->GetNextOrig() != e[2] || e[15]->GetPrevOrig() != e[4] || e[15]->GetNextDest() != e[16] || e[15]->GetPrevDest() != e[14])
		return false;
	if (e[16]->GetNextOrig() != e[9] || e[16]->GetPrevOrig() != e[1] || e[16]->GetNextDest() != e[17] || e[16]->GetPrevDest() != e[15])
		return false;
	if (e[17]->GetNextOrig() != e[10] || e[17]->GetPrevOrig() != e[9] || e[17]->GetNextDest() != e[12] || e[17]->GetPrevDest() != e[16])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[4] || e[0]->GetPrevRight() != e[3] || e[0]->GetNextLeft() != e[7] || e[0]->GetPrevLeft() != e[6])
		return false;
	if (e[1]->GetNextRight() != e[2] || e[1]->GetPrevRight() != e[4] || e[1]->GetNextLeft() != e[10] || e[1]->GetPrevLeft() != e[9])
		return false;
	if (e[2]->GetNextRight() != e[4] || e[2]->GetPrevRight() != e[1] || e[2]->GetNextLeft() != e[16] || e[2]->GetPrevLeft() != e[15])
		return false;
	if (e[3]->GetNextRight() != e[0] || e[3]->GetPrevRight() != e[4] || e[3]->GetNextLeft() != e[15] || e[3]->GetPrevLeft() != e[14])
		return false;
	if (e[4]->GetNextRight() != e[1] || e[4]->GetPrevRight() != e[2] || e[4]->GetNextLeft() != e[3] || e[4]->GetPrevLeft() != e[0])
		return false;
	if (e[5]->GetNextRight() != e[6] || e[5]->GetPrevRight() != e[8] || e[5]->GetNextLeft() != e[11] || e[5]->GetPrevLeft() != e[10])
		return false;
	if (e[6]->GetNextRight() != e[0] || e[6]->GetPrevRight() != e[7] || e[6]->GetNextLeft() != e[8] || e[6]->GetPrevLeft() != e[5])
		return false;
	if (e[7]->GetNextRight() != e[6] || e[7]->GetPrevRight() != e[0] || e[7]->GetNextLeft() != e[14] || e[7]->GetPrevLeft() != e[13])
		return false;
	if (e[8]->GetNextRight() != e[5] || e[8]->GetPrevRight() != e[6] || e[8]->GetNextLeft() != e[13] || e[8]->GetPrevLeft() != e[12])
		return false;
	if (e[9]->GetNextRight() != e[1] || e[9]->GetPrevRight() != e[10] || e[9]->GetNextLeft() != e[17] || e[9]->GetPrevLeft() != e[16])
		return false;
	if (e[10]->GetNextRight() != e[5] || e[10]->GetPrevRight() != e[11] || e[10]->GetNextLeft() != e[9] || e[10]->GetPrevLeft() != e[1])
		return false;
	if (e[11]->GetNextRight() != e[10] || e[11]->GetPrevRight() != e[5] || e[11]->GetNextLeft() != e[12] || e[11]->GetPrevLeft() != e[17])
		return false;
	if (e[12]->GetNextRight() != e[8] || e[12]->GetPrevRight() != e[13] || e[12]->GetNextLeft() != e[17] || e[12]->GetPrevLeft() != e[11])
		return false;
	if (e[13]->GetNextRight() != e[7] || e[13]->GetPrevRight() != e[14] || e[13]->GetNextLeft() != e[12] || e[13]->GetPrevLeft() != e[8])
		return false;
	if (e[14]->GetNextRight() != e[3] || e[14]->GetPrevRight() != e[15] || e[14]->GetNextLeft() != e[13] || e[14]->GetPrevLeft() != e[7])
		return false;
	if (e[15]->GetNextRight() != e[2] || e[15]->GetPrevRight() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetPrevLeft() != e[3])
		return false;
	if (e[16]->GetNextRight() != e[9] || e[16]->GetPrevRight() != e[17] || e[16]->GetNextLeft() != e[15] || e[16]->GetPrevLeft() != e[2])
		return false;
	if (e[17]->GetNextRight() != e[11] || e[17]->GetPrevRight() != e[12] || e[17]->GetNextLeft() != e[16] || e[17]->GetPrevLeft() != e[9])
		return false;
	//---------------------------------------------------------------------------------
	for (int i = 0; i < 18; ++i)
	{
		if (E[i]->edge != e[i])
			return false;
	}
	NumTests += 1;
	return true;
}
bool tester_ConvexHull3D_3(int& NumTests)
{
	//test convex hull performance when coincidental points are given in input
	vector<Vector3D> input;
	input.resize(20);
	input[2](0) = 1;
	input[3](1) = 1;
	input[4](0) = input[4](1) = 1;
	input[5](0) = input[5](1) = input[5](2) = 0.4;
	Vector3D dz(0, 0, 1);
	input[10] = input[0] + dz;
	input[12] = input[2] + dz;
	input[13] = input[3] + dz;
	input[14] = input[4] + dz;
	input[15] = input[5];
	input[19] = input[13];
	ConvexHull3D h1;
	h1.Build(input);
	QuadEdge* hull = h1.GetHull();
	if (!hull)
		return false;
	if (!hull->TestIntegrity())
		return false;
	if (!h1.TestConvexity())
		return false;
	if (hull->NumVertices() != 8 || hull->NumEdges() != 18 || hull->NumFaces() != 12)
		return false;
	NumTests += 1;
	return true;
}
bool tester_ConvexHull3D_Icosahedron(int& NumTests)
{
	double phi = (1.0 + sqrt(5.0)) / 2.0;
	vector<Vector3D> input;
	input.resize(30);
	input[1](0) = 0; input[1](1) = +1; input[1](2) = +phi;
	input[3](0) = 0; input[3](1) = +1; input[3](2) = -phi;
	input[5](0) = 0; input[5](1) = -1; input[5](2) = +phi;
	input[7](0) = 0; input[7](1) = -1; input[7](2) = -phi;
	
	input[9](0) = +phi; input[9](1) = 0; input[9](2) = +1;
	input[11](0) = -phi; input[11](1) = 0; input[11](2) = +1;
	input[13](0) = +phi; input[13](1) = 0; input[13](2) = -1;
	input[15](0) = -phi; input[15](1) = 0; input[15](2) = -1;

	input[17](0) = +1; input[17](1) = +phi; input[17](2) = 0;
	input[19](0) = +1; input[19](1) = -phi; input[19](2) = 0;
	input[21](0) = -1; input[21](1) = +phi; input[21](2) = 0;
	input[23](0) = -1; input[23](1) = -phi; input[23](2) = 0;

	input[2] = input[1];

	input[4](0) = -0.5; input[4](1) = -0.5; input[4](2) = -0.5;
	input[6](0) = -0.5; input[6](1) = -0.5; input[6](2) = +0.5;
	input[8](0) = -0.5; input[8](1) = +0.5; input[8](2) = -0.5;
	input[10](0) = -0.5; input[10](1) = +0.5; input[10](2) = +0.5;
	input[12](0) = +0.5; input[12](1) = -0.5; input[12](2) = -0.5;
	input[14](0) = +0.5; input[14](1) = -0.5; input[14](2) = +0.5;
	input[16](0) = +0.5; input[16](1) = +0.5; input[16](2) = -0.5;
	input[18](0) = +0.5; input[18](1) = +0.5; input[18](2) = +0.5;

	ConvexHull3D h1;
	h1.Build(input);
	QuadEdge* hull = h1.GetHull();
	if (!hull)
		return false;
	if (!hull->TestIntegrity())
		return false;
	if (!h1.TestConvexity())
		return false;
	if (hull->NumVertices() != 12 || hull->NumEdges() != 30 || hull->NumFaces() != 20)
		return false;
	NumTests += 1;
	return true;
}
bool tester_ConvexHull3D_Dodecahedron(int& NumTests)
{
	double phi = (1.0 + sqrt(5.0)) / 2.0;
	vector<Vector3D> input;
	input.resize(50);

	input[1](0) = -1; input[1](1) = -1; input[1](2) = -1;
	input[2](0) = -1; input[2](1) = -1; input[2](2) = +1;
	input[3](0) = -1; input[3](1) = +1; input[3](2) = -1;
	input[4](0) = -1; input[4](1) = +1; input[4](2) = +1;
	input[5](0) = +1; input[5](1) = -1; input[5](2) = -1;
	input[6](0) = +1; input[6](1) = -1; input[6](2) = +1;
	input[7](0) = +1; input[7](1) = +1; input[7](2) = -1;
	input[8](0) = +1; input[8](1) = +1; input[8](2) = +1;

	input[9](0)  = 0; input[9](1)  = +phi; input[9](2)  = +1.0 / phi;
	input[10](0) = 0; input[10](1) = +phi; input[10](2) = -1.0 / phi;
	input[11](0) = 0; input[11](1) = -phi; input[11](2) = +1.0 / phi;
	input[24](0) = 0; input[24](1) = -phi; input[24](2) = -1.0 / phi;

	input[26](0) = +1.0 / phi; input[26](1) = 0; input[26](2) = +phi;
	input[28](0) = -1.0 / phi; input[28](1) = 0; input[28](2) = +phi;
	input[30](0) = +1.0 / phi; input[30](1) = 0; input[30](2) = -phi;
	input[32](0) = -1.0 / phi; input[32](1) = 0; input[32](2) = -phi;

	input[34](0) = +phi; input[34](1) = +1.0 / phi; input[34](2) = 0;
	input[36](0) = +phi; input[36](1) = -1.0 / phi; input[36](2) = 0;
	input[38](0) = -phi; input[38](1) = +1.0 / phi; input[38](2) = 0;
	input[40](0) = -phi; input[40](1) = -1.0 / phi; input[40](2) = 0;

	input[12] = input[1];

	input[13](0) = -0.2; input[13](1) = -0.2; input[13](2) = -0.2;
	input[14](0) = -0.2; input[14](1) = -0.2; input[14](2) = +0.2;
	input[15](0) = -0.2; input[15](1) = +0.2; input[15](2) = -0.2;
	input[16](0) = -0.2; input[16](1) = +0.2; input[16](2) = +0.2;
	input[17](0) = +0.2; input[17](1) = -0.2; input[17](2) = -0.2;
	input[25](0) = +0.2; input[25](1) = -0.2; input[25](2) = +0.2;
	input[27](0) = +0.2; input[27](1) = +0.2; input[27](2) = -0.2;
	input[39](0) = +0.2; input[39](1) = +0.2; input[39](2) = +0.2;

	ConvexHull3D h1;
	h1.Build(input);
	QuadEdge* hull = h1.GetHull();
	if (!hull)
		return false;
	if (!hull->TestIntegrity())
		return false;
	if (!h1.TestConvexity())
		return false;
	if (hull->NumVertices() != 20 || hull->NumEdges() != 54 || hull->NumFaces() != 36)
		return false;
	//--------------------------------------------------------------------------------
	ConvexHull3D h2;
	for (int i = 0; i < input.size(); ++i)
		h2.AddPoint(input[i]);
	QuadEdge* hull2 = h2.GetHull();
	if (!hull2)
		return false;
	if (!hull2->TestIntegrity())
		return false;
	if (!h2.TestConvexity())
		return false;
	if (hull2->NumVertices() != 20 || hull2->NumEdges() != 54 || hull2->NumFaces() != 36)
		return false;
	//--------------------------------------------------------------------------------
	ConvexHull3D h3;
	for (int i = 0; i < 3; ++i)
		h3.AddPoint(input[i]);
	ConvexHull3D h4(h3);
	for (int i = 3; i < input.size(); ++i)
		h4.AddPoint(input[i]);
	QuadEdge* hull4 = h4.GetHull();
	if (!hull4)
		return false;
	if (!hull4->TestIntegrity())
		return false;
	if (!h4.TestConvexity())
		return false;
	if (hull4->NumVertices() != 20 || hull4->NumEdges() != 54 || hull4->NumFaces() != 36)
		return false;
	//--------------------------------------------------------------------------------
	ConvexHull3D h5;
	for (int i = 0; i < 17; ++i)
	{
		h5.AddPoint(input[i]);
	}
	if (!h5.TestConvexity())
		return false;
	QuadEdge* hull5 = h5.GetHull();
	hull5->IssueIndices();
	ConvexHull3D h6(h5);
	QuadEdge* hull6 = h6.GetHull();
	hull6->IssueIndices();
	if (*hull5 != *hull6)
		return false;
	if (!hull6)
		return false;
	if (!hull6->TestIntegrity())
		return false;
	if (!h6.TestConvexity())
		return false;
	//------Test copying quadEqge
	int Nv5 = hull5->NumVertices();
	int Nf5 = hull5->NumFaces();
	int Ne5 = hull5->NumEdges();
	int Nv6 = hull6->NumVertices();
	int Nf6 = hull6->NumFaces();
	int Ne6 = hull6->NumEdges();
	Vertex* v5[11];
	Face* f5[18];
	Edge* e5[27];
	EdgeContainer* E5[27];
	Vertex* v6[11];
	Face* f6[18];
	Edge* e6[27];
	EdgeContainer* E6[27];
	Vector3D p5[11];
	Vector3D p6[11];
	for (int i = 0; i < 11; ++i)
	{
		v5[i] = hull5->GetVertex(i);
		p5[i] = v5[i]->GetPoint();
		v6[i] = hull6->GetVertex(i);
		p6[i] = v6[i]->GetPoint();
		if (p5[i] != p6[i])
			return false;
		if (v5[i]->index == -1)
			return false;
		if (v6[i]->index == -1)
			return false;
	}
	for (int i = 0; i < 18; ++i)
	{
		f5[i] = hull5->GetFace(i);
		f6[i] = hull6->GetFace(i);
		if (f5[i]->index == -1)
			return false;
		if (f6[i]->index == -1)
			return false;
	}
	for (int i = 0; i < 27; ++i)
	{
		e5[i] = hull5->GetEdge(i);
		E5[i] = hull5->GetEdgeList()->GetContainer(i);
		if (E5[i]->edge != e5[i])
			return false;
		e6[i] = hull6->GetEdge(i);
		E6[i] = hull6->GetEdgeList()->GetContainer(i);
		if (E6[i]->edge != e6[i])
			return false;
		if (e5[i]->index == -1)
			return false;
		if (e6[i]->index == -1)
			return false;
		if (E5[i]->index == -1)
			return false;
		if (E6[i]->index == -1)
			return false;
	}
	for (int i = 0; i < 11; ++i)
	{
		int ip1 = i < 10 ? i + 1 : 0;
		int im1 = i > 0 ? i - 1 : 10;
		if (v5[i]->next != v5[ip1] || v5[i]->prev != v5[im1])
			return false;
		if (v6[i]->next != v6[ip1] || v6[i]->prev != v6[im1])
			return false;
		if (v5[i]->next->index != v6[i]->next->index)
			return false;
		if (v5[i]->prev->index != v6[i]->prev->index)
			return false;
		if (v5[i]->GetEdge()->index != v6[i]->GetEdge()->index)
			return false;
	}
	for (int i = 0; i < 18; ++i)
	{
		int ip1 = i < 17 ? i + 1 : 0;
		int im1 = i > 0 ? i - 1 : 17;
		if (f5[i]->next != f5[ip1] || f5[i]->prev != f5[im1])
			return false;
		if (f6[i]->next != f6[ip1] || f6[i]->prev != f6[im1])
			return false;
		if (f5[i]->next->index != f6[i]->next->index)
			return false;
		if (f5[i]->prev->index != f6[i]->prev->index)
			return false;
		if (f5[i]->GetEdge()->index != f6[i]->GetEdge()->index)
			return false;
	}
	for (int i = 0; i < 27; ++i)
	{
		int ip1 = i < 26 ? i + 1 : 0;
		int im1 = i > 0 ? i - 1 : 26;
		if (E5[i]->next != E5[ip1] || E5[i]->prev != E5[im1])
			return false;
		if (E6[i]->next != E6[ip1] || E6[i]->prev != E6[im1])
			return false;
		if (E5[i]->next->index != E6[i]->next->index)
			return false;
		if (E5[i]->prev->index != E6[i]->prev->index)
			return false;
		if (E5[i]->edge->index != E6[i]->edge->index)
			return false;
		if (e5[i]->GetOrig()->index != e6[i]->GetOrig()->index)
			return false;
		if (e5[i]->GetDest()->index != e6[i]->GetDest()->index)
			return false;
		if (e5[i]->GetRight()->index != e6[i]->GetRight()->index)
			return false;
		if (e5[i]->GetLeft()->index != e6[i]->GetLeft()->index)
			return false;
		if (e5[i]->GetNextOrig()->index != e6[i]->GetNextOrig()->index)
			return false;
		if (e5[i]->GetNextDest()->index != e6[i]->GetNextDest()->index)
			return false;
		if (e5[i]->GetNextRight()->index != e6[i]->GetNextRight()->index)
			return false;
		if (e5[i]->GetNextLeft()->index != e6[i]->GetNextLeft()->index)
			return false;
		if (e5[i]->GetPrevOrig()->index != e6[i]->GetPrevOrig()->index)
			return false;
		if (e5[i]->GetPrevDest()->index != e6[i]->GetPrevDest()->index)
			return false;
		if (e5[i]->GetPrevRight()->index != e6[i]->GetPrevRight()->index)
			return false;
		if (e5[i]->GetPrevLeft()->index != e6[i]->GetPrevLeft()->index)
			return false;
	}
	//------End test copying quadEdge
	for (int i = 17; i < input.size(); ++i)
	{
		h6.AddPoint(input[i]);
	}
	hull6 = h6.GetHull();
	if (!hull6)
		return false;
	if (!hull6->TestIntegrity())
		return false;
	if (hull6->NumVertices() != 20 || hull6->NumEdges() != 54 || hull6->NumFaces() != 36)
		return false;
	NumTests += 1;
	return true;//cannot get to here, crashes.
}