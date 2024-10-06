#include "Edge.h"
#include "Vertex.h"
#include "Face.h"
#include "MathUtils.h"
void Edge::Reset()
{
	this->Orig = 0;
	this->Dest = 0;
	this->Right = 0;
	this->Left = 0;
	this->nextOrig = 0;
	this->prevOrig = 0;
	this->nextDest = 0;
	this->prevDest = 0;
	this->nextRight = 0;
	this->prevRight = 0;
	this->nextLeft = 0;
	this->prevLeft = 0;
	this->index = -1;
}
void Edge::CopyBody(const Edge& e)
{
	this->Reset();
}
Edge::Edge()
{
	this->Reset();
	this->type = 0;
}
Edge::Edge(const Edge& e)
{
	this->CopyBody(e);
	this->type = 0;
}
Edge& Edge::operator=(const Edge& e)
{
	this->CopyBody(e);
	return *this;
}
void Edge::SetOrig(Vertex* v)
{
	this->Orig = v;
}
void Edge::SetDest(Vertex* v)
{
	this->Dest = v;
}
void Edge::SetRight(Face* f)
{
	this->Right = f;
}
void Edge::SetLeft(Face* f)
{
	this->Left = f;
}
void Edge::SetNextOrig(Edge* e)
{
	this->nextOrig = e;
}
void Edge::SetPrevOrig(Edge* e)
{
	this->prevOrig = e;
}
void Edge::SetNextDest(Edge* e)
{
	this->nextDest = e;
}
void Edge::SetPrevDest(Edge* e)
{
	this->prevDest = e;
}
void Edge::SetNextRight(Edge* e)
{
	this->nextRight = e;
}
void Edge::SetPrevRight(Edge* e)
{
	this->prevRight = e;
}
void Edge::SetNextLeft(Edge* e)
{
	this->nextLeft = e;
}
void Edge::SetPrevLeft(Edge* e)
{
	this->prevLeft = e;
}
Vertex* Edge::GetOrig() const
{
	return this->Orig;
}
Vertex* Edge::GetDest() const
{
	return this->Dest;
}
Face* Edge::GetRight() const
{
	return this->Right;
}
Face* Edge::GetLeft() const
{
	return this->Left;
}
Edge* Edge::GetNextOrig() const
{
	return this->nextOrig;
}
Edge* Edge::GetPrevOrig() const
{
	return this->prevOrig;
}
Edge* Edge::GetNextDest() const
{
	return this->nextDest;
}
Edge* Edge::GetPrevDest() const
{
	return this->prevDest;
}
Edge* Edge::GetNextRight() const
{
	return this->nextRight;
}
Edge* Edge::GetPrevRight() const
{
	return this->prevRight;
}
Edge* Edge::GetNextLeft() const
{
	return this->nextLeft;
}
Edge* Edge::GetPrevLeft() const
{
	return this->prevLeft;
}
Edge* Edge::GetNext(Face* f) const
{//needs testing
	//Will not work if Right = Left
	if (this->Left == f)
		return this->nextLeft;
	if (this->Right == f)
		return this->nextRight;
	return 0;
}
Edge* Edge::GetNext(Vertex* v) const
{//needs testing
	//Will not Work if Origin = Dest
	if (this->Orig == v)
		return this->nextOrig;
	if (this->Dest == v)
		return this->nextDest;
	return 0;
}
void Edge::ReverseFaces()
{
	Face* f = this->Right;
	this->Right = this->Left;
	this->Left = f;

	Edge* e = this->nextRight;
	this->nextRight = this->nextLeft;
	this->nextLeft = e;

	e = this->prevRight;
	this->prevRight = this->prevLeft;
	this->prevLeft = e;
}
void Edge::ReverseVertices()
{
	Vertex* v = this->Orig;
	this->Orig = this->Dest;
	this->Dest = v;

	Edge* e = this->nextOrig;
	this->nextOrig = this->nextDest;
	this->nextDest = e;

	e = this->prevOrig;
	this->prevOrig = this->prevDest;
	this->prevDest = e;
}
void Edge::Reverse()
{
	this->ReverseFaces();
	this->ReverseVertices();	
}
Vector3D Edge::GetVector() const
{
	if (!this->Orig || !this->Dest)
		throw "Edge::GetVector(): Origin or Dist is NULL";
	Vector3D v = this->Dest->GetPoint() - this->Orig->GetPoint();
	return v;
}
Vector3D Edge::GetNormal_2D() const //Only works for 2D edges
{
	if (!this->Orig || !this->Dest)
		throw "Edge::GetNormal(): Origin or Dist is NULL";
	Vector3D ez(0, 0, 1);
	Vector3D eta = this->GetVector();
	double etaL = eta.abs();
	if (!etaL)
		throw "Edge::GetNormal(): Edge Length is zero.";
	Vector3D e_eta = eta / etaL;
	Vector3D n = e_eta && ez;
	n(2) = 0;
	return n;
}
double Edge::GetLength() const
{
	if (this->Orig && this->Dest)
	{
		double dx = this->Dest->GetX() - this->Orig->GetX();
		double dy = this->Dest->GetY() - this->Orig->GetY();
		double dz = this->Dest->GetZ() - this->Orig->GetZ();
		return sqrt(dx * dx + dy * dy + dz * dz);
	}
	return 0.0;
}
Vector3D Edge::GetMidPoint() const
{
	return ((this->Orig->GetPoint() + this->Dest->GetPoint()) / 2.0);
}
Vector3D Edge::GetPoint()
{
	return this->GetMidPoint();
}
LineSegment2D Edge::GetLineSegment2D() const
{//needs testing
	if (!this->Orig || !this->Dest)
		throw "Edge::GetVector(): Origin or Dist is NULL.";
	LineSegment2D L(this->Orig->GetPoint(),this->Dest->GetPoint());
	return L;
}
EDGE_RELATION Edge::RelationTo(Edge* e1) const
{
	
	if (!e1)
		return EDGES_DISJOINT;
	if (this->Orig == e1->Dest && this->Dest == e1->Orig)
		return EDGES_CYCLIC;
	if (this->Orig == e1->Dest || this->Dest == e1->Orig)
		return EDGES_SEQUENTIAL;
	if (this->Orig == e1->Orig && this->Dest == e1->Dest)
		return EDGES_CONVERGENT_DIVERGENT;
	if (this->Orig == e1->Orig)
		return EDGES_DIVERGENT;
	if (this->Dest == e1->Dest)
		return EDGES_CONVERGENT;
	return EDGES_DISJOINT;
}
EDGE_RELATION Edge::DirectionTo(Edge* e1) const
{
	if (!e1)
		return EDGES_DISJOINT;
	if (this->Orig == e1->Dest || this->Dest == e1->Orig)
		return EDGES_SEQUENTIAL;
	if (this->Orig == e1->Orig || this->Dest == e1->Dest)
		return EDGES_NONSEQUENTIAL;
	return EDGES_DISJOINT;
}
void Edge::UpdateNextRight(Edge* e)
{
	this->SetNextRight(e);
	if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetPrevRight(this);
	else if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetPrevLeft(this);
}
void Edge::UpdateNextLeft(Edge* e)
{
	this->SetNextLeft(e);
	if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetPrevLeft(this);
	else if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetPrevRight(this);
}
void Edge::UpdateNextOrig(Edge* e)
{
	this->SetNextOrig(e);
	if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetPrevOrig(this);
	else if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetPrevDest(this);
}
void Edge::UpdateNextDest(Edge* e)
{
	this->SetNextDest(e);
	if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetPrevDest(this);
	else if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetPrevOrig(this);
}
void Edge::UpdatePrevRight(Edge* e)
{
	this->SetPrevRight(e);
	if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetNextRight(this);
	else if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetNextLeft(this);
}
void Edge::UpdatePrevLeft(Edge* e)
{
	this->SetPrevLeft(e);
	if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetNextLeft(this);
	else if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetNextRight(this);
}
void Edge::UpdatePrevOrig(Edge* e)
{
	this->SetPrevOrig(e);
	if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetNextOrig(this);
	else if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetNextDest(this);
}
void Edge::UpdatePrevDest(Edge* e)
{
	this->SetPrevDest(e);
	if (this->DirectionTo(e) == EDGES_NONSEQUENTIAL)
		e->SetNextDest(this);
	else if (this->DirectionTo(e) == EDGES_SEQUENTIAL)
		e->SetNextOrig(this);
}
Face* Edge::GetOtherFace(Face* f)
{
	if (f == this->Left)
		return this->Right;
	if (f == this->Right)
		return this->Left;
	return 0;
}
void Edge::CloneAdjacency(Edge* e)
{
	this->Orig = e->Orig; this->Dest = e->Dest;
	this->Left = e->Left; this->Right = e->Right;
	this->nextOrig = e->nextOrig; this->prevOrig = e->prevOrig;
	this->nextDest = e->nextDest; this->prevDest = e->prevDest;
	this->nextLeft = e->nextLeft; this->prevLeft = e->prevLeft;
	this->nextRight = e->nextRight; this->prevRight = e->prevRight;
}
EdgeBasic::EdgeBasic() : Edge()
{
	this->type = 1;
}
Edge* EdgeBasic::Clone() const
{
	return (new EdgeBasic);
}
bool tester_Edge(int& NumTests)
{
	//A C3 graph
	VertexBasic v1, v2;
	v1.SetX(11.5); v1.SetY(-21.8); v1.SetZ(14.36);
	v2.SetX(v1.GetX() + 3.0); v2.SetY(v1.GetY() + 4.0); v2.SetZ(v1.GetZ());
	EdgeBasic e0, e1, e2, e3, e4, e5, e6, e7, e8;
	FaceBasic f1, f2;
	EdgeBasic e00(e0), e000;
	e000 = e00;
	e000.SetOrig(&v1);
	e000.SetDest(&v2);
	if (fabs(e000.GetLength() - 5.0) > 1.0e-10)
		return false;
	e000.SetRight(&f1);
	e000.SetLeft(&f2);
	e000.SetNextOrig(&e1);
	e000.SetPrevOrig(&e2);
	e000.SetNextDest(&e3);
	e000.SetPrevDest(&e4);
	e000.SetNextRight(&e5);
	e000.SetPrevRight(&e6);
	e000.SetNextLeft(&e7);
	e000.SetPrevLeft(&e8);
	if (!e000.GetOrig() || e000.GetOrig() != &v1)
		return false;
	if (!e000.GetDest() || e000.GetDest() != &v2)
		return false;
	if (!e000.GetRight() || e000.GetRight() != &f1)
		return false;
	if (!e000.GetLeft() || e000.GetLeft() != &f2)
		return false;
	if (!e000.GetNextOrig() || e000.GetNextOrig() != &e1)
		return false;
	if (!e000.GetPrevOrig() || e000.GetPrevOrig() != &e2)
		return false;
	if (!e000.GetNextDest() || e000.GetNextDest() != &e3)
		return false;
	if (!e000.GetPrevDest() || e000.GetPrevDest() != &e4)
		return false;
	if (!e000.GetNextRight() || e000.GetNextRight() != &e5)
		return false;
	if (!e000.GetPrevRight() || e000.GetPrevRight() != &e6)
		return false;
	if (!e000.GetNextLeft() || e000.GetNextLeft() != &e7)
		return false;
	if (!e000.GetPrevLeft() || e000.GetPrevLeft() != &e8)
		return false;
	e000.Reverse();
	if (!e000.GetOrig() || e000.GetOrig() != &v2)
		return false;
	if (!e000.GetDest() || e000.GetDest() != &v1)
		return false;
	if (!e000.GetRight() || e000.GetRight() != &f2)
		return false;
	if (!e000.GetLeft() || e000.GetLeft() != &f1)
		return false;
	if (!e000.GetNextOrig() || e000.GetNextOrig() != &e3)
		return false;
	if (!e000.GetPrevOrig() || e000.GetPrevOrig() != &e4)
		return false;
	if (!e000.GetNextDest() || e000.GetNextDest() != &e1)
		return false;
	if (!e000.GetPrevDest() || e000.GetPrevDest() != &e2)
		return false;
	if (!e000.GetNextRight() || e000.GetNextRight() != &e7)
		return false;
	if (!e000.GetPrevRight() || e000.GetPrevRight() != &e8)
		return false;
	if (!e000.GetNextLeft() || e000.GetNextLeft() != &e5)
		return false;
	if (!e000.GetPrevLeft() || e000.GetPrevLeft() != &e6)
		return false;
	NumTests += 1;
	return true;
}