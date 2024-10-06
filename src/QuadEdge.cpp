#include <unordered_set>
#include <iomanip>
#include <sstream>
#include "QuadEdge.h"
#include "Circle.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "Vector3D.h"
#include "EdgeList.h"
#include "EdgeIterator.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "QuadEdgeIndex.h"
#include "Triangulation.h" //for testing only
#include "RuppertShewchuk.h"//for testing only
void QuadEdge::CopyBody(const QuadEdge& g)
{
	this->Reset(); 
	delete this->facePrototype;
	this->facePrototype = g.facePrototype->Clone();
	delete this->edgePrototype;
	this->edgePrototype = g.edgePrototype->Clone();
	delete this->vertexPrototype;
	this->vertexPrototype = g.vertexPrototype->Clone(); 
	EdgeIterator ite(g.edges);
	Edge* e = ite.Next();
	vector<Edge*> EdgeAssociation[2];
	int Ne = g.edges ? g.edges->GetNumEdgeContainers() : 0;
	while (e)
	{
		Edge* e1 = this->AddEdge();
		e1->index = e->index;
		EdgeAssociation[0].push_back(e);
		EdgeAssociation[1].push_back(e1);
		e = ite.Next();
	}
	EdgeIterator itee0(g.edges);
	Edge* ee0 = itee0.Next();
	EdgeIterator itee1(this);
	Edge* ee1 = itee1.Next();
	while (ee0 && ee1)
	{
		itee1.GetCurrentContainer()->index = itee0.GetCurrentContainer()->index;
		ee0 = itee0.Next();
		ee1 = itee1.Next();
	}
	Vertex* v = g.v0;
	vector<Vertex*> VertexAssociation[2];
	for (int i = 0; i < g.Nv; ++i)
	{
		if (!v)
			throw "QuadEdge::CopyBody()";
		if (i > 0 && v == g.v0)
			throw "QuadEdge::CopyBody()";
		this->nextVertexIndex = v->index;
		Vertex* v1 = this->AddVertex(v);
		v1->SetDegree(v->GetDegree());
		VertexAssociation[0].push_back(v);
		VertexAssociation[1].push_back(v1);
		Edge* e0 = v->GetEdge();
		for (int j = 0; j < Ne; ++j)
		{
			if (EdgeAssociation[0][j] == e0)
			{
				v1->SetEdge(EdgeAssociation[1][j]);
				break;
			}
		}
		v = v->next;
	}
	Face* f = g.face0;
	vector<Face*> FaceAssociation[2];
	for (int i = 0; i < g.Nf; ++i)
	{
		if (!f)
			throw "QuadEdge::CopyBody()";
		if (i > 0 && f == g.face0)
			throw "QuadEdge::CopyBody()";
		Face* f1 = this->AddFace();
		f1->index = f->index;
		FaceAssociation[0].push_back(f);
		FaceAssociation[1].push_back(f1);
		Edge* e0 = f->GetEdge();
		for (int j = 0; j < Ne; ++j)
		{
			if (EdgeAssociation[0][j] == e0)
			{
				f1->SetEdge(EdgeAssociation[1][j]);
				break;
			}
		}
		f = f->next;
	}
	EdgeIterator ite0(g.edges);
	EdgeIterator ite1(this->edges);
	Edge* e0 = ite0.Next();
	Edge* e1 = ite1.Next();
	while (e0 && e1)
	{
		Vertex* O0 = e0->GetOrig();
		Vertex* D0 = e0->GetDest();
		bool found1 = false, found2 = false;
		for (int j = 0; j < g.Nv; ++j)
		{
			if (!found1 && VertexAssociation[0][j] == O0)
			{
				e1->SetOrig(VertexAssociation[1][j]);
				found1 = true;
			}
			if (!found2 && VertexAssociation[0][j] == D0)
			{
				e1->SetDest(VertexAssociation[1][j]);
				found2 = true;
			}
			if (found1 && found2)
				break;
		}
		Face* R0 = e0->GetRight();
		Face* L0 = e0->GetLeft();
		found1 = false; found2 = false;
		for (int j = 0; j < g.Nf; ++j)
		{
			if (!found1 && FaceAssociation[0][j] == R0)
			{
				e1->SetRight(FaceAssociation[1][j]);
				found1 = true;
			}
			if (!found2 && FaceAssociation[0][j] == L0)
			{
				e1->SetLeft(FaceAssociation[1][j]);
				found2 = true;
			}
			if (found1 && found2)
				break;
		}
		Edge* nO0 = e0->GetNextOrig();
		Edge* nD0 = e0->GetNextDest();
		Edge* nR0 = e0->GetNextRight();
		Edge* nL0 = e0->GetNextLeft();
		found1 = false; found2 = false;
		bool found3 = false, found4 = false;
		for (int j = 0; j < Ne; ++j)
		{
			if (!found1 && EdgeAssociation[0][j] == nO0)
			{
				e1->SetNextOrig(EdgeAssociation[1][j]);
				found1 = true;
			}
			if (!found2 && EdgeAssociation[0][j] == nD0)
			{
				e1->SetNextDest(EdgeAssociation[1][j]);
				found2 = true;
			}
			if (!found3 && EdgeAssociation[0][j] == nR0)
			{
				e1->SetNextRight(EdgeAssociation[1][j]);
				found3 = true;
			}
			if (!found4 && EdgeAssociation[0][j] == nL0)
			{
				e1->SetNextLeft(EdgeAssociation[1][j]);
				found4 = true;
			}
			if (found1 && found2 && found3 && found4)
				break;
		}
		e0 = ite0.Next();
		e1 = ite1.Next();
	}
	this->UpdatePrevs();
	this->nextVertexIndex = -1;
	this->isPlanar = g.isPlanar;
	this->UpdateFacesDegree();
	this->UpdateVerticesDegree();
}
void QuadEdge::AddEdgeToVertex(Vertex* v, Edge* e1, bool isDest)
{//Does NOT update Prev relations
	if (!v || !e1)
		throw "QuadEdge::AddEdgeToVertex()";
	Edge* e0 = v->GetEdge();
	if (!e0)
	{
		v->SetEdge(e1);
		e0 = e1;
	}
	if (isDest)
	{
		e1->SetDest(v);
		if (e0 == e1)
			e1->SetNextDest(e1);
	}
	else
	{
		e1->SetOrig(v);
		if (e0 == e1)
			e1->SetNextOrig(e1);
	}
	Edge* e2;
	if (e0->GetOrig() == v)
	{
		e2 = e0->GetNextOrig();
		e0->UpdateNextOrig(e1);
	}
	else if (e0->GetDest() == v)
	{
		e2 = e0->GetNextDest();
		e0->UpdateNextDest(e1);
	}
	else
		throw "QuadEdge::AddEdgeToVertex()";
	if (isDest)
		e1->UpdateNextDest(e2);
	else
		e1->UpdateNextOrig(e2);
	v->IncreaseDegree();
}
void QuadEdge::AddEdgeToFace(Face* f, Edge* e1, bool isLeft)
{//Does NOT update Prev relations
	if (!f || !e1)
		throw "QuadEdge::AddFaceToVertex()";
	Edge* e0 = f->GetEdge();
	if (!e0)
	{
		f->SetEdge(e1);
		e0 = e1;
	}
	if (isLeft)
	{
		e1->SetLeft(f);
		if (e0 == e1)
			e1->SetNextLeft(e1);
	}
	else
	{
		e1->SetRight(f);
		if (e0 == e1)
			e1->SetNextRight(e1);
	}
	Edge* e2;
	if (e0->GetRight() == f)
	{
		e2 = e0->GetNextRight();
		e0->UpdateNextRight(e1);
	}
	else if (e0->GetLeft() == f)
	{
		e2 = e0->GetNextLeft();
		e0->UpdateNextLeft(e1);
	}
	else
		throw "QuadEdge::AddFaceToVertex()";
	if (isLeft)
		e1->UpdateNextLeft(e2);
	else
		e1->UpdateNextRight(e2);
	f->SetEdge(e1);
	f->IncrementDegree();
}
void QuadEdge::DeleteEdgeAndRightFace_arrangeEdgeFaces(Edge* e, Face* f)
{
	Edge* epL = e->GetPrevLeft();
	Edge* epR = e->GetPrevRight();
	Edge* enL = e->GetNextLeft();
	Edge* enR = e->GetNextRight(); 
		
	if (e->DirectionTo(epL) == EDGES_SEQUENTIAL)
		epL->UpdateNextLeft(enR);
	else if (e->DirectionTo(epL) == EDGES_NONSEQUENTIAL)
		epL->UpdateNextRight(enR);
	else
		throw "QuadEdge::DeleteEdgeAndRightFace_arrangeEdgeFaces()";

	if (e->DirectionTo(epR) == EDGES_SEQUENTIAL)
		epR->UpdateNextRight(enL);
	else if (e->DirectionTo(epR) == EDGES_NONSEQUENTIAL)
		epR->UpdateNextLeft(enL);
	else
		throw "QuadEdge::DeleteEdgeAndRightFace_arrangeEdgeFaces()";
			
	if (f->GetEdge() == e)
		f->SetEdge(e != enL ? enL : enR);
}
void QuadEdge::RemoveVertex(Vertex* v)
{
	if (!v)
		throw "QuadEdge::RemoveVertex()";
	Vertex* vn = v->next;
	Vertex* vp = v->prev;
	if (!vn || !vp)
		throw "QuadEdge::RemoveVertex()";
	if (this->v0 == v)
		this->v0 = vn;
	delete v;
	--this->Nv;
	if (!this->Nv)
		this->v0 = 0;
	else
	{
		vp->next = vn;
		vn->prev = vp;
	}
}
void QuadEdge::RemoveFace(Face* f)
{
	if (!f)
		throw "QuadEdge::RemoveFace()";
	Face* fn = f->next;
	Face* fp = f->prev;
	if (!fn || !fp)
		throw "QuadEdge::RemoveFace()";
	if (this->face0 == f)
		this->face0 = fn;//what if fn == f?
	delete f;
	--this->Nf;
	if (!this->Nf)
		this->face0 = 0;
	else
	{
		fp->next = fn;
		fn->prev = fp;
	}
}
void QuadEdge::RemoveEdge(EdgeContainer* E)
{
	if (!E)
		throw "QuadEdge::RemoveEdge()";
	Edge* e = E->edge;
	if (!e)
		throw "QuadEdge::RemoveEdge()";
	this->edges->RemoveContainer(E);
	delete e;	
}
void QuadEdge::RemoveEdge(Edge* e)
{
	EdgeContainer* E = this->edges->FindContainer(e);
	this->RemoveEdge(E);
}
void QuadEdge::SeparateEdgeAndVertex(Edge* e, Vertex* v)
{
	Edge *en, *ep;
	if (e->GetOrig() == v)
	{
		en = e->GetNextOrig();
		ep = e->GetPrevOrig();
	}
	else if (e->GetDest() == v)
	{
		en = e->GetNextDest();
		ep = e->GetPrevDest();
	}
	else
		throw "QuadEdge::SeparateEdgeAndVertex()";
	if (en->GetOrig() == v)
	{
		en->UpdatePrevOrig(ep);
	}
	else if (en->GetDest() == v)
	{
		en->UpdatePrevDest(ep);
	}
	else
		throw "QuadEdge::SeparateEdgeAndVertex()";
	if (v->GetEdge() == e)
	{
		if (e->GetOrig() == v)
			v->SetEdge(e->GetNextOrig());
		if (e->GetDest() == v)
			v->SetEdge(e->GetNextDest());
		if (v->GetEdge() == e)//relevant when vertex degree = 1
			v->SetEdge(0);
	}
	v->DecreaseDegree();
}
void QuadEdge::RemoveLoneVertices()
{
	Vertex* v = v0;
	for (int i = 0; i < Nv; ++i)
	{
		Vertex* vn = v->next;
		if (v && !v->GetEdge())
		{
			this->RemoveVertex(v);
			v = 0;
		}
		v = vn;
	}
}
QuadEdge::QuadEdge()
{
	this->v0 = 0;
	this->Nv = 0;
	this->face0 = 0;
	this->Nf = 0;
	this->edges = 0;
	this->facePrototype = new FaceBasic;
	this->edgePrototype = new EdgeBasic;
	this->vertexPrototype = new VertexBasic;
	this->nextVertexIndex = -1;
	this->isPlanar = false;
}
QuadEdge::QuadEdge(const QuadEdge& g)
{
	this->v0 = 0;
	this->Nv = 0;
	this->face0 = 0;
	this->Nf = 0;
	this->edges = 0;
	this->facePrototype = 0;
	this->edgePrototype = 0;
	this->vertexPrototype = 0;
	this->nextVertexIndex = -1;
	this->CopyBody(g);
}
QuadEdge::~QuadEdge()
{
	this->Reset();
	delete this->facePrototype;
	this->facePrototype = 0;
	delete this->edgePrototype;
	this->edgePrototype = 0;
	delete this->vertexPrototype;
	this->vertexPrototype = 0;
}
QuadEdge& QuadEdge::operator=(const QuadEdge& g)
{
	this->CopyBody(g);
	return *this;
};
void QuadEdge::SetPrototype(Face* prototype)
{
	delete this->facePrototype;
	this->facePrototype = prototype->Clone();
}
void QuadEdge::SetPrototype(Edge* prototype)
{
	delete this->edgePrototype;
	this->edgePrototype = prototype->Clone();
}
void QuadEdge::SetPrototype(Vertex* prototype)
{
	delete this->vertexPrototype;
	this->vertexPrototype = prototype->Clone();
}
void QuadEdge::Reset()
{
	if (v0 && v0->prev)
		v0->prev->next = 0;
	while (v0)
	{
		Vertex* v1 = v0->next;
		delete v0;
		v0 = v1;
	}
	this->Nv = 0;
	if (face0 && face0->prev)
		face0->prev->next = 0;
	while (face0)
	{
		Face* f1 = face0->next;
		delete face0;
		face0 = f1;
	}
	this->Nf = 0;
	if (this->edges)
	{
		this->edges->ResetAndDeleteEdges();
		delete this->edges;
	}
	this->edges = 0;
	if(!this->facePrototype)
		this->facePrototype = new FaceBasic;
	if (!this->edgePrototype)
		this->edgePrototype = new EdgeBasic;
	if (!this->vertexPrototype)
		this->vertexPrototype = new VertexBasic;
	this->nextVertexIndex = -1;
	this->isPlanar = false;
}
Vertex* QuadEdge::AddVertex(const Vector3D& v)
{
	++this->Nv;
	if (!v0)
	{
		v0 = this->vertexPrototype->Clone();
		v0->index = this->nextVertexIndex;
		*v0 << v;
		v0->next = v0;
		v0->prev = v0;
		return v0;
	}
	else
	{
		Vertex* v1 = v0->prev;
		v0->prev = this->vertexPrototype->Clone();
		v0->prev->index = this->nextVertexIndex;
		*(v0->prev) << v;
		v0->prev->next = v0;
		v0->prev->prev = v1;
		v1->next = v0->prev;
		return v0->prev;
	}
}
Vertex* QuadEdge::AddVertex(const Vertex* v)
{
	Vector3D u = v->GetPoint();
	return this->AddVertex(u);
}
Vertex* QuadEdge::GetVertex(int vi)
{
	if (vi < 0 || (vi > 0 && vi >= this->Nv))
		throw "QuadEdge::GetVertex()";
	Vertex* v = v0;
	for (int i = 0; i < vi; ++i)
	{
		if (!v)
			throw "QuadEdge::GetVertex()";
		v = v->next;
	}
	return v;
}
Vertex* QuadEdge::GetVertexOfIndex(int vi)
{
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	while (v)
	{
		if (v->index == vi)
			return v;
		v = itv.Next();
	}
	return 0;
}
Face* QuadEdge::AddFace()
{
	++this->Nf;
	if (!face0)
	{
		face0 = this->facePrototype->Clone();
		face0->next = face0;
		face0->prev = face0;
		return face0;
	}
	else
	{
		Face* f1 = face0->prev;
		face0->prev = this->facePrototype->Clone();
		face0->prev->next = face0;
		face0->prev->prev = f1;
		f1->next = face0->prev;
		return face0->prev;
	}
}
Edge* QuadEdge::AddEdge()
{
	Edge* e = this->edgePrototype->Clone();
	if (!this->edges)
		this->edges = new EdgeList;
	this->edges->AddContainer(e);
	return e;
}
Edge* QuadEdge::GetEdge(int ei)
{
	if (!this->edges)
		return 0;
	EdgeContainer* E = this->edges->GetContainer(ei);
	if (!E)
		return 0;
	return E->edge;
}
Face* QuadEdge::GetFace(int fi)
{
	if (fi < 0 || (fi > 0 && fi >= this->Nf))
		throw "QuadEdge::GetFace()";
	Face* f = face0;
	for (int i = 0; i < fi; ++i)
	{
		if (!f)
			throw "QuadEdge::GetFace()";
		f = f->next;
	}
	return f;
}
void QuadEdge::UpdatePrevs()
{
	EdgeIterator ite(this);
	Edge* e0 = ite.Next();
	while (e0)
	{
		Edge* enR = e0->GetNextRight();
		int dir_enR = e0->DirectionTo(enR);
		if (dir_enR == EDGES_SEQUENTIAL)
			enR->SetPrevRight(e0);
		else if (dir_enR == EDGES_NONSEQUENTIAL)
			enR->SetPrevLeft(e0);
		else
			throw "QuadEdge::UpdatePrevs()";

		Edge* enL = e0->GetNextLeft();
		int dir_enL = e0->DirectionTo(enL);
		if (dir_enL == EDGES_NONSEQUENTIAL)
			enL->SetPrevRight(e0);
		else if (dir_enL == EDGES_SEQUENTIAL)
			enL->SetPrevLeft(e0);
		else
			throw "QuadEdge::UpdatePrevs()";

		Edge* enO = e0->GetNextOrig();
		int dir_enO = e0->DirectionTo(enO);
		if (dir_enO == EDGES_NONSEQUENTIAL)
			enO->SetPrevOrig(e0);
		else if (dir_enO == EDGES_SEQUENTIAL)
			enO->SetPrevDest(e0);
		else
			throw "QuadEdge::UpdatePrevs()";

		Edge* enD = e0->GetNextDest();
		int dir_enD = e0->DirectionTo(enD);
		if (dir_enD == EDGES_SEQUENTIAL)
			enD->SetPrevOrig(e0);
		else if (dir_enD == EDGES_NONSEQUENTIAL)
			enD->SetPrevDest(e0);
		else
			throw "QuadEdge::UpdatePrevs()";

		e0 = ite.Next();
	} 
}
int QuadEdge::NumVertices() const
{
	return this->Nv;
}
int QuadEdge::NumEdges() const
{
	return this->edges ? this->edges->GetNumEdgeContainers() : 0;
}
int QuadEdge::NumFaces() const
{
	return this->Nf;
}
void QuadEdge::GetFaceEdges(Face* f, vector<Edge*>& edgeVec)
{
	edgeVec.resize(0);
	EdgeIterator ite(f);
	Edge* e = ite.Next();
	while (e)
	{
		edgeVec.push_back(e);
		e = ite.Next();
	}
}
void QuadEdge::GetVertexEdges(Vertex* v, vector<Edge*>& edgeVec)
{
	edgeVec.resize(0);
	EdgeIterator ite(v);
	Edge* e = ite.Next();
	while (e)
	{
		edgeVec.push_back(e);
		e = ite.Next();
	}
}
void QuadEdge::GetFaceEdgesReverse(Face* f, vector<Edge*>& edgeVec)
{
	edgeVec.resize(0);
	Edge* e0 = f->GetEdge();
	if (e0->GetLeft() == e0->GetRight())
	{
		if (e0->GetPrevLeft() != e0)
			e0 = e0->GetPrevLeft();
		else
			e0 = e0->GetPrevRight();
	}
	bool onLeft = (e0->GetLeft() == f);
	Edge* e = e0;
	do
	{
		edgeVec.push_back(e);
		Vertex* orig = e->GetOrig();
		Vertex* dest = e->GetDest();
		if (onLeft)
			e = e->GetPrevLeft();
		else
			e = e->GetPrevRight();
		if (e->GetOrig() == orig || e->GetDest() == dest)
			onLeft = !onLeft;
	} while (e != e0);
}
void QuadEdge::GetVertexEdgesReverse(Vertex* v, vector<Edge*>& edgeVec)
{
	edgeVec.resize(0);
	Edge* e0 = v->GetEdge();
	Edge* e = e0;
	do
	{
		edgeVec.push_back(e);
		if (e->GetOrig() == v)
			e = e->GetPrevOrig();
		else if (e->GetDest() == v)
			e = e->GetPrevDest();
		else
			throw "QuadEdge::GetVertexEdgesReverse()";
	} while (e != e0);
}
void QuadEdge::SetAsDoubleTriangle(const Vector3D& V0, const Vector3D& V1, const Vector3D& V2)
{
	this->Reset();
	Vertex* A = this->AddVertex(V0);
	Vertex* B = this->AddVertex(V1);
	Vertex* C = this->AddVertex(V2);
	Face* ABC = this->AddFace();
	Face* ACB = this->AddFace();
	Edge* eAB = this->AddEdge();
	Edge* eBC = this->AddEdge();
	Edge* eCA = this->AddEdge();
	this->AddEdgeToVertex(A, eAB, false);
	this->AddEdgeToVertex(A, eCA, true);
	this->AddEdgeToVertex(B, eBC, false);
	this->AddEdgeToVertex(B, eAB, true);
	this->AddEdgeToVertex(C, eCA, false);
	this->AddEdgeToVertex(C, eBC, true);
	this->AddEdgeToFace(ABC, eBC, true);
	this->AddEdgeToFace(ABC, eCA, true);
	this->AddEdgeToFace(ABC, eAB, true);
	this->AddEdgeToFace(ACB, eCA, false);
	this->AddEdgeToFace(ACB, eBC, false);
	this->AddEdgeToFace(ACB, eAB, false);
	this->UpdatePrevs();
}
void QuadEdge::SetAsDoubleTriangulatedSquare(const Vector3D& V0, const Vector3D& V1, const Vector3D& V2, const Vector3D& V3)
{
	this->Reset();
	Vector3D V01 = V1 - V0;
	Vector3D V02 = V2 - V0;
	Vector3D V03 = V3 - V0;
	Vertex* A, *B, *C, *D;
	if (((V02 && V01) || (V02 && V03)) < 0)
	{
		A = this->AddVertex(V0);
		B = this->AddVertex(V1);
		C = this->AddVertex(V2);
		D = this->AddVertex(V3);
	}
	else if (((V03 && V01) || (V03 && V02)) < 0)
	{
		A = this->AddVertex(V0);
		B = this->AddVertex(V1);
		C = this->AddVertex(V3);
		D = this->AddVertex(V2);
	}
	else
	{
		A = this->AddVertex(V0);
		B = this->AddVertex(V3);
		C = this->AddVertex(V1);
		D = this->AddVertex(V2);
	}
	Face* ABC = this->AddFace();
	Face* ACB = this->AddFace();
	Face* ACD = this->AddFace();
	Face* ADC = this->AddFace();
	Edge* eAB = this->AddEdge();
	Edge* eBC = this->AddEdge();
	Edge* eCA = this->AddEdge();
	Edge* eAC = this->AddEdge();
	Edge* eCD = this->AddEdge();
	Edge* eDA = this->AddEdge();
	this->AddEdgeToVertex(A, eAB, false);
	this->AddEdgeToVertex(A, eCA, true);
	this->AddEdgeToVertex(A, eDA, true);
	this->AddEdgeToVertex(A, eAC, false);
	this->AddEdgeToVertex(B, eBC, false);
	this->AddEdgeToVertex(B, eAB, true);
	this->AddEdgeToVertex(C, eCA, false);
	this->AddEdgeToVertex(C, eBC, true);
	this->AddEdgeToVertex(C, eCD, false);
	this->AddEdgeToVertex(C, eAC, true);
	this->AddEdgeToVertex(D, eCD, true);
	this->AddEdgeToVertex(D, eDA, false);
	this->AddEdgeToFace(ABC, eAB, true);
	this->AddEdgeToFace(ABC, eBC, true);
	this->AddEdgeToFace(ABC, eCA, true);
	this->AddEdgeToFace(ACB, eAB, false);
	this->AddEdgeToFace(ACB, eAC, true);
	this->AddEdgeToFace(ACB, eBC, false);
	this->AddEdgeToFace(ACD, eDA, true);
	this->AddEdgeToFace(ACD, eCA, false);
	this->AddEdgeToFace(ACD, eCD, true);
	this->AddEdgeToFace(ADC, eDA, false);
	this->AddEdgeToFace(ADC, eCD, false);
	this->AddEdgeToFace(ADC, eAC, false);
	this->UpdatePrevs();
}
Vertex* QuadEdge::AddPyramidAndDeleteFace(const Vector3D& v, Face* f0)
{
	if (!f0)
		throw "QuadEdge::AddPyramidAndDeleteFace()";
	Vertex* p = this->AddVertex(v);
	vector<Edge*> edgeVec;
	vector<Face*> newFaces;
	vector<Edge*> newEdges;
	this->GetFaceEdges(f0, edgeVec);
	p->SetDegree(edgeVec.size());
	this->MakeAllEdgesLeftTo(f0);
	for (int i = 0; i < edgeVec.size(); ++i)
	{
		Face* fNew = this->AddFace();
		newFaces.push_back(fNew);
		Edge* eNew = this->AddEdge();
		newEdges.push_back(eNew);
	}
	p->SetEdge(newEdges[edgeVec.size() - 1]);
	for (int i = 0; i < edgeVec.size(); ++i)
	{
		int ip1 = (i < edgeVec.size() - 1) ? i + 1 : 0;
		int im1 = (i != 0) ? i - 1 : edgeVec.size() - 1;
		Vertex* v = edgeVec[i]->GetOrig();
		newEdges[i]->SetDest(p);
		this->AddEdgeToVertex(v, newEdges[i], false);
	}
	for (int i = 0; i < edgeVec.size(); ++i)
	{
		int ip1 = (i < edgeVec.size() - 1) ? i + 1 : 0;
		int im1 = (i != 0) ? i - 1 : edgeVec.size() - 1;
		edgeVec[i]->SetLeft(newFaces[i]);
		newFaces[i]->SetEdge(edgeVec[i]);
		newEdges[i]->SetRight(newFaces[i]);
		newEdges[i]->SetLeft(newFaces[im1]);
		newEdges[i]->UpdateNextDest(newEdges[ip1]);
		newEdges[i]->UpdateNextRight(edgeVec[i]);
		newEdges[i]->UpdateNextLeft(newEdges[im1]);
		edgeVec[i]->UpdateNextLeft(newEdges[ip1]);
	}
	this->RemoveFace(f0);
	return p;
}
bool QuadEdge::SetAsTetrahedron(const Vector3D& V0, const Vector3D& V1, const Vector3D& V2, const Vector3D& V3)
{
	this->SetAsDoubleTriangle(V0, V1, V2);
	Edge* eAB = this->GetEdge(0);
	Face* ABC = eAB->GetLeft();//does not need to be the same as in SetAsDoubleTriangle(...)
	Face* ACB = eAB->GetRight();//does not need to be the same as in SetAsDoubleTriangle(...)
	Edge* eBC = eAB->GetNextLeft();
	if (eBC->GetLeft() != ABC && eBC->GetRight() != ABC)
		throw "QuadEdge::SetAsTetrahedron()";
	if (eBC->GetLeft() != ACB && eBC->GetRight() != ACB)
		throw "QuadEdge::SetAsTetrahedron()";
	if (eBC->GetLeft() != ABC)
		eBC->Reverse();
	Edge* eCA = eBC->GetNextLeft();
	if (eCA->GetLeft() != ABC && eCA->GetRight() != ABC)
		throw "QuadEdge::SetAsTetrahedron()";
	if (eCA->GetLeft() != ACB && eCA->GetRight() != ACB)
		throw "QuadEdge::SetAsTetrahedron()";
	if (eCA->GetLeft() != ABC)
		eCA->Reverse();
	Vector3D A = eAB->GetOrig()->GetPoint();
	Vector3D B = eBC->GetOrig()->GetPoint();
	Vector3D C = eCA->GetOrig()->GetPoint();
	Vector3D AB = B - A;
	Vector3D AC = C - A;
	Vector3D N = AB && AC;
	Vector3D AD = V3 - A;
	double NdotAD = N || AD;
	if (!NdotAD)
		return false;
	if (NdotAD < 0)
	{
		this->AddPyramidAndDeleteFace(V3, ACB);
	}
	else
	{
		this->AddPyramidAndDeleteFace(V3, ABC);
	}
	return true;
	
}
void QuadEdge::DeleteEdgeAndRightFace(EdgeContainer* E)
{//Right face will be absorbed in the left face. Updates Prev relations.
	Edge* e = E->edge;
	Face* fR = e->GetRight();
	Face* fL = e->GetLeft();
	if (!this->isRemovable(e))
		throw "QuadEdge::DeleteEdgeAndRightFace() : Euler's formula is broken.";
	vector<Edge*> edgeVec;
	this->GetFaceEdges(fR, edgeVec);
	this->SeparateEdgeAndVertex(e, e->GetOrig());
	this->SeparateEdgeAndVertex(e, e->GetDest());
	this->DeleteEdgeAndRightFace_arrangeEdgeFaces(e, fL);
	for (int i = 0; i < edgeVec.size(); ++i)
	{
		if (edgeVec[i]->GetRight() == fR)
			edgeVec[i]->SetRight(fL);
		else
			edgeVec[i]->SetLeft(fL);
	}
	if (fR != fL)
	{
		this->RemoveFace(fR); 
	}
	else
		this->RemoveLoneVertices();
	this->RemoveEdge(E);
}
void QuadEdge::DeleteEdgeAndRightFace(Edge* e)
{
	EdgeContainer* E = this->edges->FindContainer(e);
	this->DeleteEdgeAndRightFace(E);
}
Vector3D QuadEdge::GetNormal(Face* f) const
{
	VertexIterator itv(f);
	Vertex* v0 = itv.Next();
	Vertex* v1 = itv.Next();
	Vertex* v2 = itv.Next();
	if (!v0 || !v1 || !v2)
		throw "QuadEdge::GetNormal()";
	Vector3D O = v0->GetPoint();
	Vector3D n;
	while (n.abs() < 1.0e-10 && v2)
	{
		Vector3D u1 = v1->GetPoint() - v0->GetPoint();
		Vector3D u2 = v2->GetPoint() - v1->GetPoint();
		n = u1 && u2;
		v0 = v1;
		v1 = v2;
		v2 = itv.Next();
	}
	if (!n.abs())
		throw "QuadEdge::GetNormal()";
	n = n / n.abs();
	return n;

}
POINT_FACE_VISIBILITY QuadEdge::isVisible(Face* f, const Vector3D p) const
{
	Vertex* v0 = f->GetEdge()->GetOrig();
	if (!v0)
		throw "QuadEdge::isVisible()";
	Vector3D O = v0->GetPoint();
	Vector3D n = this->GetNormal(f);
	Vector3D OP = p - O;
	double nDotOP = OP || n;
	if (fabs(nDotOP) < 1.0e-15)
		nDotOP = 0;
	if (nDotOP > 0)
		return VISIBLE;
	if (nDotOP == 0)
		return COPLANAR;
	return INVISIBLE;
}
int QuadEdge::UpdateVertexDegree(Vertex* v)
{
	v->SetDegree(0);
	EdgeIterator ite(v);
	Edge* e = ite.Next();
	while (e)
	{
		v->IncreaseDegree();
		e = ite.Next();
	}
	return v->GetDegree();
}
void QuadEdge::UpdateVerticesDegree()
{
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	while (v)
	{
		this->UpdateVertexDegree(v);
		v = itv.Next();
	}
}
int QuadEdge::GetSumVertexDegrees() 
{
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	int SumDv = 0;
	while (v)
	{
		SumDv += v->GetDegree();
		v = itv.Next();
	}
	return SumDv;
}
int QuadEdge::GetMaxVertexDegree()
{
	this->UpdateVerticesDegree();
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	int MaxDv = 0;
	while (v)
	{
		if(v->GetDegree() > MaxDv)
			MaxDv = v->GetDegree();
		v = itv.Next();
	}
	return MaxDv;
}
int QuadEdge::GetMinVertexDegree()
{
	this->UpdateVerticesDegree();
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	int MinDv = 10000;
	while (v)
	{
		if (v->GetDegree()< MinDv)
			MinDv = v->GetDegree();
		v = itv.Next();
	}
	return MinDv;
}
int QuadEdge::UpdateFaceDegree(Face* f)
{
	f->SetDegree(0);
	EdgeIterator ite(f);
	Edge* e = ite.Next();
	while (e)
	{
		f->IncrementDegree();
		e = ite.Next();
	}
	return f->GetDegree();
}
void QuadEdge::UpdateFacesDegree()
{
	FaceIterator itf(this);
	Face* f = itf.Next();
	while (f)
	{
		this->UpdateFaceDegree(f);
		f = itf.Next();
	}
}
int QuadEdge::GetSumFaceDegrees()
{
	FaceIterator itf(this);
	Face* f = itf.Next();
	int SumDf = 0;
	while (f)
	{
		SumDf += f->GetDegree();
		f = itf.Next();
	}
	return SumDf;
}
bool QuadEdge::isRemovable(Edge* e) const
{
	int de = -1;
	int dv = (e->GetOrig()->GetDegree() == 1) ? -1 : 0;
	dv += (e->GetDest()->GetDegree() == 1) ? -1 : 0;
	int df = (e->GetRight() != e->GetLeft()) ? -1 : 0;
	int dN = dv - de + df;
	return (dN == 0);
}
EdgeList* QuadEdge::GetEdgeList()
{
	return this->edges;
}
void QuadEdge::SetNextVetexIndex(int index)
{
	this->nextVertexIndex = index;
}
void QuadEdge::IssueIndices()
{
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	int maxVertexIndex = -1;
	while (v)
	{
		if (v->index > maxVertexIndex)
			maxVertexIndex = v->index;
		v = itv.Next();
	}
	VertexIterator itv1(this);
	v = itv1.Next();
	while (v)
	{
		if (v->index == -1)
			v->index = ++maxVertexIndex;
		v = itv1.Next();
	}
	FaceIterator itf(this);
	Face* f = itf.Next();
	int maxFaceIndex = -1;
	while (f)
	{
		f->index = ++maxFaceIndex;
		f = itf.Next();
	}
	EdgeIterator ite(this);
	Edge* e = ite.Next();
	int maxEdgeIndex = -1;
	while (e)
	{
		e->index = ++maxEdgeIndex;
		ite.GetCurrentContainer()->index = maxEdgeIndex;
		e = ite.Next();
	}
}
bool QuadEdge::IsPlanar() const
{
	return this->isPlanar;
}
void QuadEdge::SetPlanar(bool v)
{
	this->isPlanar = v;
}
bool QuadEdge::operator==(const QuadEdge& rhs) const
{
	if (this->NumVertices() != rhs.NumVertices())
		return false;
	if (this->NumEdges() != rhs.NumEdges())
		return false;
	if (this->NumEdges() != rhs.NumEdges())
		return false;
	QuadEdge h0 = *this;
	QuadEdge h1 = rhs;
	h0.IssueIndices();
	h1.IssueIndices();
	VertexIterator itv0(&h0);
	VertexIterator itv1(&h1);
	Vertex* v0 = itv0.Next();
	Vertex* v1 = itv1.Next();
	while (v0 && v1)
	{
		if (v0->index != v1->index)
			return false;
		Vector3D p0 = v0->GetPoint();
		Vector3D p1 = v1->GetPoint();
		if (p0 != p1)
			return false;
		if (bool(v0->GetEdge()) !=  bool(v1->GetEdge()))
			return false;
		if (v0->GetEdge() && v1->GetEdge() && v0->GetEdge()->index != v1->GetEdge()->index)
			return false;
		v0 = itv0.Next();
		v1 = itv1.Next();
	}
	FaceIterator itf0(&h0);
	FaceIterator itf1(&h1);
	Face* f0 = itf0.Next();
	Face* f1 = itf1.Next();
	while (f0 && f1)
	{
		if (f0->index != f1->index)
			return false;
		if (bool(f0->GetEdge()) != bool(f1->GetEdge()))
			return false;
		if (f0->GetEdge() && f1->GetEdge() && f0->GetEdge()->index != f1->GetEdge()->index)
			return false;
		f0 = itf0.Next();
		f1 = itf1.Next();
	}
	EdgeIterator ite0(&h0);
	EdgeIterator ite1(&h1);
	Edge* e0 = ite0.Next();
	Edge* e1 = ite1.Next();
	while (e0 && e1)
	{
		if (e0->index != e1->index)
			return false;
		if (e0->GetOrig()->index != e1->GetOrig()->index)
			return false;
		if (e0->GetDest()->index != e1->GetDest()->index)
			return false;
		if (e0->GetRight()->index != e1->GetRight()->index)
			return false;
		if (e0->GetLeft()->index != e1->GetLeft()->index)
			return false;
		if (e0->GetNextOrig()->index != e1->GetNextOrig()->index)
			return false;
		if (e0->GetNextDest()->index != e1->GetNextDest()->index)
			return false;
		if (e0->GetNextRight()->index != e1->GetNextRight()->index)
			return false;
		if (e0->GetNextLeft()->index != e1->GetNextLeft()->index)
			return false;
		if (e0->GetPrevOrig()->index != e1->GetPrevOrig()->index)
			return false;
		if (e0->GetPrevDest()->index != e1->GetPrevDest()->index)
			return false;
		if (e0->GetPrevRight()->index != e1->GetPrevRight()->index)
			return false;
		if (e0->GetPrevLeft()->index != e1->GetPrevLeft()->index)
			return false;
		e0 = ite0.Next();
		e1 = ite1.Next();
	}
	return true;
}
bool QuadEdge::operator!=(const QuadEdge& rhs) const
{
	return !(this->operator==(rhs));
}
void QuadEdge::ReverseFace(Face* f)
{
	this->UpdateFaceDegree(f);
	int Df = f->GetDegree();
	vector<Edge*> edges;
	edges.resize(Df, 0);
	EdgeIterator ite(f);
	Edge* e = ite.Next();
	int i = 0;
	while (e){
		if (i >= Df)
			throw "QuadEdge::ReverseFace()";
		edges[i] = e;
		if (!e)
			throw "QuadEdge::ReverseFace()";
		++i;
		e = ite.Next();
	}
	for (int j = 0; j < Df; ++j){
		Edge* e = edges[j];
		if (f == e->GetLeft())
		{
			Edge* buffer = e->GetNextLeft();
			e->SetNextLeft(e->GetPrevLeft());
			e->SetPrevLeft(buffer);
		}
		if (f == e->GetRight())
		{
			Edge* buffer = e->GetNextRight();
			e->SetNextRight(e->GetPrevRight());
			e->SetPrevRight(buffer);
		}
	}
}
void QuadEdge::MakeAllEdgesLeftTo(Face* f)
{
	if (!f)
		throw "QuadEdge::MakeAllEdgesLeftTo()";
	this->UpdateFaceDegree(f);
	int Df = f->GetDegree();
	vector<Edge*> edVec;
	edVec.resize(Df, 0);
	EdgeIterator ite(f);
	Edge* e = ite.Next();
	int i = 0;
	while (e){
		if (i >= Df)
			throw "QuadEdge::MakeAllEdgesLeftTo()";
		edVec[i] = e;
		++i;
		e = ite.Next();
	}
	if (i != Df)
		throw "QuadEdge::MakeAllEdgesLeftTo()";
	for (int i = 0; i < Df; ++i){
		if (edVec[i]->GetLeft() != f && edVec[i]->GetRight() != f)
			throw "QuadEdge::MakeAllEdgesLeftTo()";
		if (edVec[i]->GetLeft() != f)
			edVec[i]->Reverse();
	}
}
Edge* QuadEdge::DivideFaceOnRight(Face* f, Vertex* v1, Vertex* v2)
{
	//The new face will be on the right of new v1-v2 edge. 
	//Recall: in this code of quad edge, all faces are CCW when looked at from outside of the solid
	this->MakeAllEdgesLeftTo(f);
	Edge* e1L = 0;//the edge left of v1
	Edge* e1R = 0;//the edge right of v1
	Edge* e2L = 0;//the edge left of v2
	Edge* e2R = 0;//the edge right of v2
	EdgeIterator ite(f);
	Edge* e = ite.Next();
	while (e){
		if (e->GetDest() == v1){
			if (e1L)
				throw "QuadEdge::DivideFaceOnRight()";
			e1L = e;
		}
		else if (e->GetOrig() == v1){
			if (e1R)
				throw "QuadEdge::DivideFaceOnRight()";
			e1R = e;
		}
		else if (e->GetOrig() == v2){
			if (e2L)
				throw "QuadEdge::DivideFaceOnRight()";
			e2L = e;
		}
		else if (e->GetDest() == v2){
			if (e2R)
				throw "QuadEdge::DivideFaceOnRight()";
			e2R = e;
		}
		e = ite.Next();
	}
	if (!e1L || !e1R || !e2L || !e2R)
		throw "QuadEdge::DivideFaceOnRight()";
	Edge* eNew = this->AddEdge();
	Face* fNew = this->AddFace();
	this->AddEdgeToVertex(v1, eNew, false);
	this->AddEdgeToVertex(v2, eNew, true);
	eNew->SetLeft(f);
	eNew->SetRight(fNew);
	fNew->SetEdge(eNew);
	f->SetEdge(eNew);
	e1L->SetNextLeft(eNew); eNew->SetPrevLeft(e1L);
	eNew->SetNextLeft(e2L); e2L->SetPrevLeft(eNew);
	e2R->SetNextLeft(eNew); eNew->SetPrevRight(e2R);
	eNew->SetNextRight(e1R); e1R->SetPrevLeft(eNew);
	this->UpdateFaceDegree(fNew);
	this->UpdateFaceDegree(f);
	this->UpdateVertexDegree(v1);
	this->UpdateVertexDegree(v2);
	e = e1R;
	while (e != eNew) {
		e->SetLeft(fNew);
		e = e->GetNextLeft();
	}
	return eNew;
}
bool QuadEdge::IsConnected(Vertex* v1, Vertex* v2) const
{
	//this and GetEdgeConnecting functions are copy of each other.
	if (!v1 || !v2)
		throw "QuadEdge::IsConnected()";
	EdgeIterator ite(v1);
	Edge* e = ite.Next();
	while (e)
	{
		if (e->GetDest() == v2 || e->GetOrig() == v2)
			return true;
		e = ite.Next();
	}
	return false;
}
Edge* QuadEdge::GetEdgeConnecting(Vertex* v1, Vertex* v2)
{
	//this and GetEdgeConnecting functions are copy of each other.
	if (!v1 || !v2)
		throw "QuadEdge::IsConnected()";
	EdgeIterator ite(v1);
	Edge* e = ite.Next();
	while (e)
	{
		if (e->GetDest() == v2 || e->GetOrig() == v2)
			return e;
		e = ite.Next();
	}
	return 0;
}
void QuadEdge::print(string address)
{
	ofstream file;
	file.open(address.c_str(), std::fstream::out);
	this->print(file);
	file.close();
}
void QuadEdge::print(ofstream& file)
{
	file << endl << endl << "-------ADDRESSES ----------------";
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D p = v->GetPoint();
		file << endl << "v" << v->index << "[" << v << "]";
		file << " #" << v->GetDegree();
		file << " (x,y,z) = (" << p(0) << " , " << p(1) << "," << p(2) << ")";
		file << " --> e" << v->GetEdge()->index << "[" << v->GetEdge() << "]";
		file << endl << "		Faces: ";
		FaceIterator itvf(v);
		Face* vf = itvf.Next();
		while (vf) {
			file << "f" << vf->index << "[" << vf << "]-->";
			vf = itvf.Next();
			if (!vf)
				file << "Return";
		}
		file << endl << "		Edges: ";
		EdgeIterator itev(v);
		Edge* ev = itev.Next();
		while (ev) {
			file << "e" << ev->index << "[" << ev << "]-->";
			ev = itev.Next();
			if (!ev)
				file << "Return";
		}
		v = itv.Next();
	}
	file << endl << "----------------------------------------";
	FaceIterator itf(this);
	Face* f = itf.Next();
	while (f)
	{
		file << endl << "f" << f->index << "[" << f << "]";
		file << " #" << f->GetDegree();
		file << " --> e" << f->GetEdge()->index << "[" << f->GetEdge() << "]";
		file << endl << "		Vertices: ";
		VertexIterator itvf(f);
		Vertex* vf = itvf.Next();
		while (vf) {
			file << "V" << vf->index << "[" << vf << "]-->";
			vf = itvf.Next();
			if (!vf)
				file << "Return";
		}
		file << endl << "		Edges: ";
		EdgeIterator itef(f);
		Edge* ef = itef.Next();
		while (ef) {
			file << "e" << ef->index << "[" << ef << "]-->";
			ef = itef.Next();
			if (!ef)
				file << "Return";
		}
		f = itf.Next();
	}
	file << endl << "----------------------------------------";
	EdgeIterator ite(this);
	Edge* e = ite.Next();
	while (e)
	{
		file << endl << "e" << e->index << "[" << e << "]";
		file << ", v" << e->GetOrig()->index << "[" << e->GetOrig() << "]";
		file << "-->v" << e->GetDest()->index << "[" << e->GetDest() << "]";
		file << ", fL" << e->GetLeft()->index << "[" << e->GetLeft() << "]";
		file << ", fR" << e->GetRight()->index << "[" << e->GetRight() << "]";
		e = ite.Next();
	}
	file << endl << "----------------------------------------";
	EdgeIterator ite1(this);
	e = ite1.Next();
	while (e)
	{
		file << endl << "e" << e->index << "[" << e << "]";
		file << ", nextOrig = e" << e->GetNextOrig()->index << "[" << e->GetNextOrig() << "]";
		file << ", nextDest = e" << e->GetNextDest()->index << "[" << e->GetNextDest() << "]";
		file << ", prevOrig = e" << e->GetPrevOrig()->index << "[" << e->GetPrevOrig() << "]";
		file << ", prevDest = e" << e->GetPrevDest()->index << "[" << e->GetPrevDest() << "]";
		e = ite1.Next();
	}
	file << endl << "----------------------------------------";
	EdgeIterator ite2(this);
	e = ite2.Next();
	while (e)
	{
		file << endl << "e" << e->index << "[" << e << "]";
		file << ", nextLeft = e" << e->GetNextLeft()->index << "[" << e->GetNextLeft() << "]";
		file << ", nextRight = e" << e->GetNextRight()->index << "[" << e->GetNextRight() << "]";
		file << ", prevLeft = e" << e->GetPrevLeft()->index << "[" << e->GetPrevLeft() << "]";
		file << ", prevRight = e" << e->GetPrevRight()->index << "[" << e->GetPrevRight() << "]";
		e = ite2.Next();
	}
}
void QuadEdge::PrintTestScript(ofstream& file)
{
	file << endl << "if(T.GetMesh2D()->NumVertices() != " << this->NumVertices();
	file << " || T.GetMesh2D()->NumEdges() != " << this->NumEdges();
	file << " || T.GetMesh2D()->NumFaces() != " << this->NumFaces() << ")";
	file << endl << "	return false;";
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	while (v)
	{
		file << endl << "if(v[" << v->index << "]->GetEdge() != e[" << v->GetEdge()->index << "]";
		file << " || v[" << v->index << "]->GetDegree() != " << v->GetDegree() << ")";
		file << endl << "	return false;";
		v = itv.Next();
	}
	FaceIterator itf(this);
	Face* f = itf.Next();
	while (f)
	{
		file << endl << "if(f[" << f->index << "]->GetEdge() != e[" << f->GetEdge()->index << "]";
		file << " || f[" << f->index << "]->GetDegree() != " << f->GetDegree() << ")";
		file << endl << "	return false;";
		f = itf.Next();
	}
	EdgeIterator ite(this);
	Edge* e = ite.Next();
	while (e)
	{
		file << endl << "if(e[" << e->index << "]->GetOrig() != v[" << e->GetOrig()->index << "]";
		file << " || e[" << e->index << "]->GetDest() != v[" << e->GetDest()->index << "]";
		file << " || e[" << e->index << "]->GetLeft() != f[" << e->GetLeft()->index << "]";
		file << " || e[" << e->index << "]->GetRight() != f[" << e->GetRight()->index << "])";
		file << endl << "	return false;";
		file << endl << "if(e[" << e->index << "]->GetNextOrig() != e[" << e->GetNextOrig()->index << "]";
		file << " || e[" << e->index << "]->GetNextDest() != e[" << e->GetNextDest()->index << "]";
		file << " || e[" << e->index << "]->GetNextLeft() != e[" << e->GetNextLeft()->index << "]";
		file << " || e[" << e->index << "]->GetNextRight() != e[" << e->GetNextRight()->index << "])";
		file << endl << "	return false;";
		file << endl << "if(e[" << e->index << "]->GetPrevOrig() != e[" << e->GetPrevOrig()->index << "]";
		file << " || e[" << e->index << "]->GetPrevDest() != e[" << e->GetPrevDest()->index << "]";
		file << " || e[" << e->index << "]->GetPrevLeft() != e[" << e->GetPrevLeft()->index << "]";
		file << " || e[" << e->index << "]->GetPrevRight() != e[" << e->GetPrevRight()->index << "])";
		file << endl << "	return false;";
		e = ite.Next();
	}
}
void QuadEdge::PrintTestScript(string address)
{
	ofstream file;
	file.open(address.c_str(), std::fstream::app);
	this->PrintTestScript(file);
	file.close();
}
bool QuadEdge::TestExclusivity()//O(n^2)
{
	int Nv = this->NumVertices();
	int Nf = this->NumFaces();
	int Ne = this->NumEdges();
	unordered_set<Vertex*> vSet;
	unordered_set<Face*> fSet;
	unordered_set<Edge*> eSet;
	vSet.rehash(Nv);
	fSet.rehash(Nf);
	eSet.rehash(Ne);
	VertexIterator itv(this);
	Vertex* vv = itv.Next();
	int Nv_prime = 0;
	while (vv)
	{
		++Nv_prime;
		vSet.insert(vv);
		vv = itv.Next();
	}
	if (Nv_prime != Nv)
		return false;
	FaceIterator itf(this);
	Face* ff = itf.Next();
	int Nf_prime = 0;
	while (ff)
	{
		++Nf_prime;
		fSet.insert(ff);
		ff = itf.Next();
	}
	if (Nf_prime != Nf)
		return false;
	EdgeIterator ite(this);
	Edge* ee = ite.Next();
	int Ne_prime = 0;
	while (ee)
	{
		++Ne_prime;
		eSet.insert(ee);
		ee = ite.Next();
	}
	if (Ne_prime != Ne)
		return false;
	VertexIterator itvp(this);
	vv = itvp.Next();
	while (vv)
	{
		Edge* ve = vv->GetEdge();
		unordered_set<Edge*>::iterator it = eSet.find(ve);
		if (it == eSet.end())
			return false;
		vv = itvp.Next();
	}
	FaceIterator itfp(this);
	ff = itfp.Next();
	while (ff)
	{
		Edge* vf = ff->GetEdge();
		unordered_set<Edge*>::iterator it = eSet.find(vf);
		if (it == eSet.end())
			return false;
		ff = itfp.Next();
	}
	EdgeIterator itep(this);
	ee = itep.Next();
	while (ee)
	{
		Vertex* vO = ee->GetOrig();
		unordered_set<Vertex*>::iterator itv1 = vSet.find(vO);
		if (itv1 == vSet.end())
			return false;
		Vertex* vD = ee->GetDest();
		unordered_set<Vertex*>::iterator itv2 = vSet.find(vD);
		if (itv2 == vSet.end())
			return false;
		Face* fR = ee->GetRight();
		unordered_set<Face*>::iterator itf1 = fSet.find(fR);
		if (itf1 == fSet.end())
			return false;
		Face* fL = ee->GetLeft();
		unordered_set<Face*>::iterator itf2 = fSet.find(fL);
		if (itf2 == fSet.end())
			return false;
		Edge* nO = ee->GetNextOrig();
		unordered_set<Edge*>::iterator ite1 = eSet.find(nO);
		if (ite1 == eSet.end())
			return false;
		Edge* pO = ee->GetPrevOrig();
		unordered_set<Edge*>::iterator ite2 = eSet.find(pO);
		if (ite2 == eSet.end())
			return false;
		Edge* nD = ee->GetNextDest();
		unordered_set<Edge*>::iterator ite3 = eSet.find(nD);
		if (ite3 == eSet.end())
			return false;
		Edge* pD = ee->GetPrevDest();
		unordered_set<Edge*>::iterator ite4 = eSet.find(pD);
		if (ite4 == eSet.end())
			return false;
		Edge* nR = ee->GetNextRight();
		unordered_set<Edge*>::iterator ite5 = eSet.find(nR);
		if (ite5 == eSet.end())
			return false;
		Edge* pR = ee->GetPrevRight();
		unordered_set<Edge*>::iterator ite6 = eSet.find(pR);
		if (ite6 == eSet.end())
			return false;
		Edge* nL = ee->GetNextLeft();
		unordered_set<Edge*>::iterator ite7 = eSet.find(nL);
		if (ite7 == eSet.end())
			return false;
		Edge* pL = ee->GetPrevLeft();
		unordered_set<Edge*>::iterator ite8 = eSet.find(pL);
		if (ite8 == eSet.end())
			return false;
		ee = itep.Next();
	}
	return true;
}
bool QuadEdge::TestIntegrity()
{
	int Nv = this->NumVertices();
	int Nf = this->NumFaces();
	int Ne = this->NumEdges();
	if (!Nv && !Ne && !Nf)
		return true;
	if (!this->TestExclusivity())
		return false;
	if (Nv - Ne + Nf != 2)
		return false;
	this->UpdateVerticesDegree();
	if (this->GetSumVertexDegrees() != 2 * Ne)
		return false;
	this->UpdateFacesDegree();
	if (this->GetSumFaceDegrees() != 2 * Ne)
		return false;
	return true;
}
Vertex* QuadEdge::InsertVertexAtMid(Edge* e)
{
	// vO---->----vNew---->----vD
	//		e			 eNew
	Vertex* vO = e->GetOrig();	Vertex* vD = e->GetDest();
	Face* fL = e->GetLeft();	Face* fR = e->GetRight();
	Edge* enO = e->GetNextOrig();	Edge* epO = e->GetPrevOrig();
	Edge* enD = e->GetNextDest();	Edge* epD = e->GetPrevDest();
	Edge* enL = e->GetNextLeft();	Edge* epL = e->GetPrevLeft();
	Edge* enR = e->GetNextRight();	Edge* epR = e->GetPrevRight();
	Vertex* vNew = this->AddVertex(e->GetPoint());
	Edge* eNew = this->AddEdge();
	vNew->SetEdge(eNew);
	eNew->CloneAdjacency(e);
	e->SetDest(vNew);			eNew->SetOrig(vNew);
	eNew->UpdateNextOrig(e);	//Does also e->SetNextDest(eNew);
	eNew->UpdatePrevOrig(e);	//Does also e->SetPrevDest(eNew);
	eNew->UpdateNextDest(enD);	eNew->UpdatePrevDest(epD);
	if (vD->GetEdge() == e)
		vD->SetEdge(eNew);
	if (enL->GetOrig() == vO || enL->GetDest() == vO) {
		e->SetPrevLeft(eNew);	//Already done in CloneAdjacency e->SetNextLeft(enL);
		eNew->UpdateNextLeft(e);	eNew->UpdatePrevLeft(epL);
	}
	else {
		e->SetNextLeft(eNew);	//Already done in CloneAdjacency e->SetPrevLeft(epL);
		eNew->UpdateNextLeft(enL);	eNew->UpdatePrevLeft(e);
	}
	if (enR->GetOrig() == vO || enR->GetDest() == vO) {
		e->SetPrevRight(eNew);	//Already done in CloneAdjacency e->SetNextRight(enR);
		eNew->UpdateNextRight(e);	eNew->UpdatePrevRight(epR);
	}
	else {
		e->SetNextRight(eNew);	//Already done in CloneAdjacency e->SetPrevRight(epR);
		eNew->UpdateNextRight(enR);	eNew->UpdatePrevRight(e);
	}
	vNew->SetDegree(2);
	fR->IncrementDegree();
	fL->IncrementDegree();
	return vNew;
}
double QuadEdge::GetMinEdgeLength()
{
	EdgeIterator ite(this);
	Edge* e = ite.Next();
	double minLength = -1;
	while (e)
	{
		double L = e->GetLength();
		if (minLength < 0 || L < minLength)
			minLength = L;
		e = ite.Next();
	}
	return minLength;
}
double QuadEdge::GetMaxEdgeLength()
{
	EdgeIterator ite(this);
	Edge* e = ite.Next();
	double maxLength = -1;
	while (e)
	{
		double L = e->GetLength();
		if (L > maxLength)
			maxLength = L;
		e = ite.Next();
	}
	return maxLength;
}
double QuadEdge::GetAvgEdgeLength()
{
	EdgeIterator ite(this);
	Edge* e = ite.Next();
	double avgLength = 0;
	int counter = 0;
	while (e)
	{
		++counter;
		double L = e->GetLength();
		avgLength += L;
		e = ite.Next();
	}
	return counter ? avgLength / double(counter) : 0 ;
}
double QuadEdge::GetWeight()
{
	EdgeIterator ite(this);
	Edge* e = ite.Next();
	double weight = 0;
	while (e)
	{
		weight += e->GetLength();
		e = ite.Next();
	}
	return weight;
}
void QuadEdge::Write(string fileName)
{
	ofstream file;
	file.open(fileName.c_str(), fstream::out);
	file << (*this);
	file.close();
}
void QuadEdge::Read(string fileName)
{
	ifstream file;
	file.open(fileName.c_str(), fstream::in); 
	file >> (*this);
	file.close();
}
void QuadEdge::GetPoints(vector<Vector3D>& points)
{
	points.resize(0);
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D P = v->GetPoint();
		points.push_back(P);
		v = itv.Next();
	}
}
ostream& operator<<(ostream& out, QuadEdge& qe)
{
	QuadEdgeIndex qui(&qe);
	int Nv = qe.NumVertices();
	int Nf = qe.NumFaces();
	int Ne = qe.NumEdges();
	out << Nv << "," << Nf << "," << Ne << ";";
	VertexIterator itv(&qe);
	Vertex* v = itv.Next();
	for (int i = 0; i < Nv; ++i)
	{
		Vertex* v = qui.GetVertex(i);
		if (v)
		{
			Vector3D P = v->GetPoint();
			out << endl << std::setprecision(10) << P(0) << "," << P(1) << "," << P(2) << "," << qui.GetIndex(v->GetEdge()) << ";";
		}
		else
			out << endl << "NULL";

	}
	for (int i = 0; i < Nf; ++i)
	{
		Face* f = qui.GetFace(i);
		if (f)
		{
			out << endl << qui.GetIndex(f->GetEdge()) << ";";
		}
		else
			out << endl << "NULL";
	}
	for (int i = 0; i < Ne; ++i)
	{
		Edge* e = qui.GetEdge(i);
		if (e)
		{
			out << endl << qui.GetIndex(e->GetOrig()) << "," << qui.GetIndex(e->GetDest()) << ",";
			out << qui.GetIndex(e->GetLeft()) << "," << qui.GetIndex(e->GetRight()) << ",";
			out << qui.GetIndex(e->GetNextOrig()) << "," << qui.GetIndex(e->GetNextDest()) << ",";
			out << qui.GetIndex(e->GetPrevOrig()) << "," << qui.GetIndex(e->GetPrevDest()) << ",";
			out << qui.GetIndex(e->GetNextLeft()) << "," << qui.GetIndex(e->GetNextRight()) << ",";
			out << qui.GetIndex(e->GetPrevLeft()) << "," << qui.GetIndex(e->GetPrevRight()) << ";";
		}
		else
			out << endl << "NULL";
	}
	return out;
}
istream& operator>>(istream& in, QuadEdge& qe)
{
	string str;
	std::getline(in, str, ',');	int Nv = stoi(str);
	std::getline(in, str, ',');	int Nf = stoi(str);
	std::getline(in, str, ';');	int Ne = stoi(str);
	vector<Vertex*> v;
	vector<Face*> f;
	vector<Edge*> e;
	v.resize(Nv, 0);
	e.resize(Ne, 0);
	f.resize(Nf, 0);
	for (int i = 0; i < Nf; ++i)
		f[i] = qe.AddFace();
	for (int i = 0; i < Ne; ++i)
		e[i] = qe.AddEdge();
	for (int i = 0; i < Nv; ++i)
	{
		std::getline(in, str, ',');		double x = stod(str);
		std::getline(in, str, ',');		double y = stod(str);
		std::getline(in, str, ',');		double z = stod(str);
		Vector3D P(x, y, z);
		v[i] = qe.AddVertex(P);
		std::getline(in, str, ';');		int j = stoi(str);		v[i]->SetEdge(e[j]);
	}
	for (int i = 0; i < Nf; ++i)
	{
		std::getline(in, str, ';');
		int j = stoi(str);
		f[i]->SetEdge(e[j]);
	}
	for (int i = 0; i < Ne; ++i)
	{
		std::getline(in, str, ',');		int j = stoi(str);		e[i]->SetOrig(v[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetDest(v[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetLeft(f[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetRight(f[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetNextOrig(e[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetNextDest(e[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetPrevOrig(e[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetPrevDest(e[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetNextLeft(e[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetNextRight(e[j]);
		std::getline(in, str, ',');		j = stoi(str);			e[i]->SetPrevLeft(e[j]);
		std::getline(in, str, ';');		j = stoi(str);			e[i]->SetPrevRight(e[j]);
	}
	return in;
}
bool tester_QuadEdge(int& NumTests)
{
	VertexBasic v0, v1(1), v2(0, 1), v3(0, 0, 1);
	QuadEdge Tet;
	Tet.AddVertex(&v0);
	Tet.AddVertex(&v1);
	Tet.AddVertex(&v2);
	Tet.AddVertex(&v3);
	QuadEdge Tet2(Tet);
	QuadEdge Tet3;
	Tet3 = Tet2;
	Vertex* v00 = Tet3.GetVertex(0);
	Vertex* v11 = Tet3.GetVertex(1);
	Vertex* v22 = Tet3.GetVertex(2);
	Vertex* v33 = Tet3.GetVertex(3);
	if (v00->GetPoint() != v0.GetPoint() || v11->GetPoint() != v1.GetPoint() 
		|| v22->GetPoint() != v2.GetPoint() || v33->GetPoint() != v3.GetPoint())
		return false;
	if (Tet != Tet2)
		return false; 
	if (Tet2 != Tet3)
		return false; 
	//-----------------------------------------------------------------------------
	QuadEdge doubleTriangle;
	doubleTriangle.SetAsDoubleTriangle(v0.GetPoint(), v1.GetPoint(), v2.GetPoint());
	if (!doubleTriangle.TestIntegrity())
		return false;
	if (doubleTriangle.NumEdges() != 3 || doubleTriangle.NumFaces() != 2 || doubleTriangle.NumVertices() != 3)
		return false;
	Edge* ed[3];
	ed[0] = doubleTriangle.GetEdge(0);
	ed[1] = doubleTriangle.GetEdge(1);
	ed[2] = doubleTriangle.GetEdge(2);
	if (ed[0] == ed[1] || ed[1] == ed[2] || ed[2] == ed[0])
		return false;
	for (int i = 0; i < 3; ++i)
	{
		if (ed[i] == ed[i]->GetNextRight() || ed[i] == ed[i]->GetNextLeft() || ed[i] == ed[i]->GetNextOrig() || ed[i] == ed[i]->GetNextDest())
			return false;
		if (ed[i] == ed[i]->GetPrevRight() || ed[i] == ed[i]->GetPrevLeft() || ed[i] == ed[i]->GetPrevOrig() || ed[i] == ed[i]->GetPrevDest())
			return false;
	}
	vector<Edge*> edgeVector0, edgeVector1;
	doubleTriangle.GetFaceEdges(doubleTriangle.GetFace(0), edgeVector0);
	doubleTriangle.GetFaceEdges(doubleTriangle.GetFace(1), edgeVector1);
	if (edgeVector0.size() != 3 || edgeVector1.size() != 3)
		return false;
	if (edgeVector0[0] != ed[0] || edgeVector0[1] != ed[1] || edgeVector0[2] != ed[2])
		return false;
	if (edgeVector1[0] != ed[0] || edgeVector1[1] != ed[2] || edgeVector1[2] != ed[1])
		return false;
	vector<Edge*> edgeVectorA, edgeVectorB, edgeVectorC;
	doubleTriangle.GetVertexEdges(doubleTriangle.GetVertex(0), edgeVectorA);
	doubleTriangle.GetVertexEdges(doubleTriangle.GetVertex(1), edgeVectorB);
	doubleTriangle.GetVertexEdges(doubleTriangle.GetVertex(2), edgeVectorC);
	if (!doubleTriangle.IsConnected(doubleTriangle.GetVertex(0), doubleTriangle.GetVertex(1)))
		return false;
	if (!doubleTriangle.IsConnected(doubleTriangle.GetVertex(1), doubleTriangle.GetVertex(0)))
		return false;
	if (!doubleTriangle.IsConnected(doubleTriangle.GetVertex(0), doubleTriangle.GetVertex(2)))
		return false;
	if (!doubleTriangle.IsConnected(doubleTriangle.GetVertex(2), doubleTriangle.GetVertex(0)))
		return false;
	if (!doubleTriangle.IsConnected(doubleTriangle.GetVertex(1), doubleTriangle.GetVertex(2)))
		return false;
	if (!doubleTriangle.IsConnected(doubleTriangle.GetVertex(2), doubleTriangle.GetVertex(1)))
		return false;
	for (int i = 0; i < 3; ++i)
	{
		if (doubleTriangle.GetVertex(i)->GetDegree() != 2)
			return false;
		if (doubleTriangle.UpdateVertexDegree(doubleTriangle.GetVertex(i)) != 2)
			return false;
	}
	if (edgeVectorA.size() != 2 || edgeVectorB.size() != 2 || edgeVectorC.size() != 2)
		return false;
	if (edgeVectorA[0] != ed[0] || edgeVectorA[1] != ed[2])
		return false;
	if (edgeVectorB[0] != ed[1] || edgeVectorB[1] != ed[0])
		return false;
	if (edgeVectorC[0] != ed[2] || edgeVectorC[1] != ed[1])
		return false;
	Vector3D P(11.5, 12, +5.5);
	Vector3D Q(11.5, 12, -5.5);
	Face* ABC = doubleTriangle.GetFace(0);
	Face* ACB = doubleTriangle.GetFace(1);
	if (doubleTriangle.isVisible(ABC, P) != VISIBLE || doubleTriangle.isVisible(ABC, Q) != INVISIBLE)
		return false;
	if (doubleTriangle.isVisible(ACB, P) != INVISIBLE || doubleTriangle.isVisible(ACB, Q) != VISIBLE)
		return false;
	if (!doubleTriangle.isRemovable(ed[0]) || !doubleTriangle.isRemovable(ed[1]) || !doubleTriangle.isRemovable(ed[2]))
		return false;
	QuadEdge QE_copy = doubleTriangle;
	if (QE_copy != doubleTriangle)
		return false;
	ed[0]->Reverse();
	if (QE_copy == doubleTriangle)
		return false;
	if (!doubleTriangle.TestIntegrity())
		return false;
	if (doubleTriangle.isVisible(ABC, P) != VISIBLE || doubleTriangle.isVisible(ABC, Q) != INVISIBLE)
		return false;
	if (doubleTriangle.isVisible(ACB, P) != INVISIBLE || doubleTriangle.isVisible(ACB, Q) != VISIBLE)
		return false;
	if (!doubleTriangle.isRemovable(ed[0]) || !doubleTriangle.isRemovable(ed[1]) || !doubleTriangle.isRemovable(ed[2]))
		return false;
	ed[1]->Reverse();
	if (!doubleTriangle.TestIntegrity())
		return false;
	if (doubleTriangle.isVisible(ABC, P) != VISIBLE || doubleTriangle.isVisible(ABC, Q) != INVISIBLE)
		return false;
	if (doubleTriangle.isVisible(ACB, P) != INVISIBLE || doubleTriangle.isVisible(ACB, Q) != VISIBLE)
		return false;
	if (!doubleTriangle.isRemovable(ed[0]) || !doubleTriangle.isRemovable(ed[1]) || !doubleTriangle.isRemovable(ed[2]))
		return false;
	ed[2]->Reverse();
	if (!doubleTriangle.TestIntegrity())
		return false;
	if (doubleTriangle.isVisible(ABC, P) != VISIBLE || doubleTriangle.isVisible(ABC, Q) != INVISIBLE)
		return false;
	if (doubleTriangle.isVisible(ACB, P) != INVISIBLE || doubleTriangle.isVisible(ACB, Q) != VISIBLE)
		return false;
	if (!doubleTriangle.isRemovable(ed[0]) || !doubleTriangle.isRemovable(ed[1]) || !doubleTriangle.isRemovable(ed[2]))
		return false;
	//-----------------------------------------------------------------------------
	//Reverse Test
	doubleTriangle.GetFaceEdgesReverse(doubleTriangle.GetFace(0), edgeVector0);
	doubleTriangle.GetFaceEdgesReverse(doubleTriangle.GetFace(1), edgeVector1);
	if (edgeVector0.size() != 3 || edgeVector1.size() != 3)
		return false;
	if (edgeVector0[0] != ed[0] || edgeVector0[1] != ed[2] || edgeVector0[2] != ed[1])
		return false;
	if (edgeVector1[0] != ed[0] || edgeVector1[1] != ed[1] || edgeVector1[2] != ed[2])
		return false;
	doubleTriangle.GetVertexEdgesReverse(doubleTriangle.GetVertex(0), edgeVectorA);
	doubleTriangle.GetVertexEdgesReverse(doubleTriangle.GetVertex(1), edgeVectorB);
	doubleTriangle.GetVertexEdgesReverse(doubleTriangle.GetVertex(2), edgeVectorC);
	if (edgeVectorA.size() != 2 || edgeVectorB.size() != 2 || edgeVectorC.size() != 2)
		return false;
	if (edgeVectorA[0] != ed[0] || edgeVectorA[1] != ed[2])
		return false;
	if (edgeVectorB[0] != ed[1] || edgeVectorB[1] != ed[0])
		return false;
	if (edgeVectorC[0] != ed[2] || edgeVectorC[1] != ed[1])
		return false;
	//-----------------------------------------------------------------------------
	//Tetrahedron
	QuadEdge Tetrahedron;
	Tetrahedron.SetAsTetrahedron(v0.GetPoint(), v1.GetPoint(), v2.GetPoint(), v3.GetPoint());//A,B,C,D
	if (!Tetrahedron.TestIntegrity())
		return false;
	Edge* e[6];//eAB, eBC, eCA
	Face* f[4];
	Vertex* v[4];
	vector<Vector3D> points;
	points.push_back(v0.GetPoint());
	points.push_back(v1.GetPoint());
	points.push_back(v2.GetPoint());
	points.push_back(v3.GetPoint());
	if (!tester_QuadEdge_Tetrahedron(Tetrahedron, v, f, e, points))
		return false;
	QuadEdge Tet_Copy = Tetrahedron;
	if (Tet_Copy != Tetrahedron)
		return false;
	//-----------------------------------------------------------------------------
	//Delete a face
	vector<Edge*> edgeVerts[4];
	vector<Edge*> edgeFaces[4];
	Tetrahedron.DeleteEdgeAndRightFace(e[4]);
	if (Tet_Copy == Tetrahedron)
		return false;
	if (!Tetrahedron.TestIntegrity())
		return false;
	if (Tetrahedron.NumVertices() != 4 || Tetrahedron.NumEdges() != 5 || Tetrahedron.NumFaces() != 3)
		return false;
	if (v[0]->GetPoint() != v0.GetPoint() || v[1]->GetPoint() != v1.GetPoint() 
		|| v[2]->GetPoint() != v2.GetPoint() || v[3]->GetPoint() != v3.GetPoint())
		return false;
	for (int i = 0; i < 4; ++i)
	{
		Tetrahedron.GetVertexEdges(v[i], edgeVerts[i]);
		if ((i == 1 || i == 3) && edgeVerts[i].size() != 2)
			return false;
		if ((i == 0 || i == 2) && edgeVerts[i].size() != 3)
			return false;
		if ((i == 1 || i == 3) && v[i]->GetDegree() != 2)
			return false;
		if ((i == 0 || i == 2) && v[i]->GetDegree() != 3)
			return false;
		if ((i == 1 || i == 3) && Tetrahedron.UpdateVertexDegree(v[i]) != 2)
			return false;
		if ((i == 0 || i == 2) && Tetrahedron.UpdateVertexDegree(v[i]) != 3)
			return false;
	}
	for (int i = 0; i < 3; ++i)
	{
		Tetrahedron.GetFaceEdges(Tetrahedron.GetFace(i), edgeFaces[i]);
		if (i == 1 && edgeFaces[i].size() != 4)
			return false;
		if (i != 1 && edgeFaces[i].size() != 3)
			return false;
	}
	for (int i = 0; i < 6; ++i)
		if (i != 4 && !Tetrahedron.isRemovable(e[i]))
			return false;
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2] || v[3]->GetEdge() != e[5])
		return false;
	if (f[0]->GetEdge() != e[0] || f[1]->GetEdge() != e[0] || f[3]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[1])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[3])
		return false;
	if (e[3]->GetOrig() != v[0] || e[3]->GetDest() != v[3] || e[3]->GetRight() != f[1] || e[3]->GetLeft() != f[3])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[1])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[3] || e[0]->GetPrevOrig() != e[2] || e[0]->GetNextDest() != e[1] || e[0]->GetPrevDest() != e[1])
		return false;
	if (e[1]->GetNextOrig() != e[0] || e[1]->GetPrevOrig() != e[0] || e[1]->GetNextDest() != e[2] || e[1]->GetPrevDest() != e[5])
		return false;
	if (e[2]->GetNextOrig() != e[5] || e[2]->GetPrevOrig() != e[1] || e[2]->GetNextDest() != e[0] || e[2]->GetPrevDest() != e[3])
		return false;
	if (e[3]->GetNextOrig() != e[2] || e[3]->GetPrevOrig() != e[0] || e[3]->GetNextDest() != e[5] || e[3]->GetPrevDest() != e[5])
		return false;
	if (e[5]->GetNextOrig() != e[1] || e[5]->GetPrevOrig() != e[2] || e[5]->GetNextDest() != e[3] || e[5]->GetPrevDest() != e[3])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[2] || e[0]->GetPrevRight() != e[1] || e[0]->GetNextLeft() != e[1] || e[0]->GetPrevLeft() != e[3])
		return false; 
	if (e[1]->GetNextRight() != e[0] || e[1]->GetPrevRight() != e[2] || e[1]->GetNextLeft() != e[5] || e[1]->GetPrevLeft() != e[0])
		return false;
	if (e[2]->GetNextRight() != e[1] || e[2]->GetPrevRight() != e[0] || e[2]->GetNextLeft() != e[3] || e[2]->GetPrevLeft() != e[5])
		return false;
	if (e[3]->GetNextRight() != e[0] || e[3]->GetPrevRight() != e[5] || e[3]->GetNextLeft() != e[5] || e[3]->GetPrevLeft() != e[2])
		return false;
	if (e[5]->GetNextRight() != e[2] || e[5]->GetPrevRight() != e[3] || e[5]->GetNextLeft() != e[3] || e[5]->GetPrevLeft() != e[1])
		return false;
	QuadEdge Tet_Copy2 = Tetrahedron;
	if (Tet_Copy2 != Tetrahedron)
		return false;
	//-----------------------------------------------------------------------------
	//Delete another face
	e[0]->Reverse();
	if (!Tetrahedron.TestIntegrity())
		return false;
	Tetrahedron.DeleteEdgeAndRightFace(e[0]);
	if (Tet_Copy2 == Tetrahedron)
		return false;
	if (!Tetrahedron.TestIntegrity())
		return false;
	if (Tetrahedron.NumVertices() != 4 || Tetrahedron.NumEdges() != 4 || Tetrahedron.NumFaces() != 2)
		return false;
	if (v[0]->GetPoint() != v0.GetPoint() || v[1]->GetPoint() != v1.GetPoint()
		|| v[2]->GetPoint() != v2.GetPoint() || v[3]->GetPoint() != v3.GetPoint())
		return false;
	for (int i = 0; i < 4; ++i)
	{
		Tetrahedron.GetVertexEdges(v[i], edgeVerts[i]);
	}
	if (edgeVerts[0].size() != 2 || v[0]->GetDegree() != 2 || Tetrahedron.UpdateVertexDegree(v[0]) != 2)
		return false;
	if (edgeVerts[1].size() != 1 || v[1]->GetDegree() != 1 || Tetrahedron.UpdateVertexDegree(v[1]) != 1)
		return false;
	if (edgeVerts[2].size() != 3 || v[2]->GetDegree() != 3 || Tetrahedron.UpdateVertexDegree(v[2]) != 3)
		return false;
	if (edgeVerts[3].size() != 2 || v[3]->GetDegree() != 2 || Tetrahedron.UpdateVertexDegree(v[3]) != 2)
		return false;
	for (int i = 0; i < 2; ++i)
	{
		Tetrahedron.GetFaceEdges(Tetrahedron.GetFace(i), edgeFaces[i]);
	}
	if (edgeFaces[0].size() != 5)
		return false;
	if (edgeFaces[1].size() != 3)
		return false;
	for (int i = 0; i < 6; ++i)
		if (i != 4 && i != 0 && !Tetrahedron.isRemovable(e[i]))
			return false;
	if (v[0]->GetEdge() != e[3] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2] || v[3]->GetEdge() != e[5])
		return false;
	if (f[0]->GetEdge() != e[2] || f[3]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[0])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[3])
		return false;
	if (e[3]->GetOrig() != v[0] || e[3]->GetDest() != v[3] || e[3]->GetRight() != f[0] || e[3]->GetLeft() != f[3])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[0])
		return false;
	//-----------------------------------------------------------------------------
	if (e[1]->GetNextOrig() != e[1] || e[1]->GetPrevOrig() != e[1] || e[1]->GetNextDest() != e[2] || e[1]->GetPrevDest() != e[5])
		return false;
	if (e[2]->GetNextOrig() != e[5] || e[2]->GetPrevOrig() != e[1] || e[2]->GetNextDest() != e[3] || e[2]->GetPrevDest() != e[3])
		return false;
	if (e[3]->GetNextOrig() != e[2] || e[3]->GetPrevOrig() != e[2] || e[3]->GetNextDest() != e[5] || e[3]->GetPrevDest() != e[5])
		return false;
	if (e[5]->GetNextOrig() != e[1] || e[5]->GetPrevOrig() != e[2] || e[5]->GetNextDest() != e[3] || e[5]->GetPrevDest() != e[3])
		return false;
	//-----------------------------------------------------------------------------
	if (e[1]->GetNextRight() != e[1] || e[1]->GetPrevRight() != e[2] || e[1]->GetNextLeft() != e[5] || e[1]->GetPrevLeft() != e[1])
		return false;
	if (e[2]->GetNextRight() != e[1] || e[2]->GetPrevRight() != e[3] || e[2]->GetNextLeft() != e[3] || e[2]->GetPrevLeft() != e[5])
		return false;
	if (e[3]->GetNextRight() != e[2] || e[3]->GetPrevRight() != e[5] || e[3]->GetNextLeft() != e[5] || e[3]->GetPrevLeft() != e[2])
		return false;
	if (e[5]->GetNextRight() != e[2] || e[5]->GetPrevRight() != e[3] || e[5]->GetNextLeft() != e[3] || e[5]->GetPrevLeft() != e[1])
		return false;
	QuadEdge Tet_Copy3 = Tetrahedron;
	if (Tet_Copy3 != Tetrahedron)
		return false;
	//-----------------------------------------------------------------------------
	//Delete another face
	Tetrahedron.DeleteEdgeAndRightFace(e[1]);
	if (Tet_Copy3 == Tetrahedron)
		return false;
	if (!Tetrahedron.TestIntegrity())
		return false;
	if (Tetrahedron.NumVertices() != 3 || Tetrahedron.NumEdges() != 3 || Tetrahedron.NumFaces() != 2)
		return false;
	if (v[0]->GetPoint() != v0.GetPoint()
		|| v[2]->GetPoint() != v2.GetPoint() || v[3]->GetPoint() != v3.GetPoint())
		return false;
	for (int i = 0; i < 3; ++i)
	{
		Tetrahedron.GetVertexEdges(Tetrahedron.GetVertex(i), edgeVerts[i]);
	}
	if (edgeVerts[0].size() != 2 || v[0]->GetDegree() != 2 || Tetrahedron.UpdateVertexDegree(v[0]) != 2)
		return false;
	if (edgeVerts[1].size() != 2 || v[2]->GetDegree() != 2 || Tetrahedron.UpdateVertexDegree(v[2]) != 2)
		return false;
	if (edgeVerts[2].size() != 2 || v[3]->GetDegree() != 2 || Tetrahedron.UpdateVertexDegree(v[3]) != 2)
		return false;
	for (int i = 0; i < 2; ++i)
	{
		Tetrahedron.GetFaceEdges(Tetrahedron.GetFace(i), edgeFaces[i]);
	}
	if (edgeFaces[0].size() != 3)
		return false;
	if (edgeFaces[1].size() != 3)
		return false;
	for (int i = 0; i < 6; ++i)
		if (i != 4 && i != 0 && i != 1 && !Tetrahedron.isRemovable(e[i]))
			return false;
	if (v[0]->GetEdge() != e[3] || v[2]->GetEdge() != e[2] || v[3]->GetEdge() != e[5])
		return false;
	if (f[0]->GetEdge() != e[2] || f[3]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[3])
		return false;
	if (e[3]->GetOrig() != v[0] || e[3]->GetDest() != v[3] || e[3]->GetRight() != f[0] || e[3]->GetLeft() != f[3])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[0])
		return false;
	//-----------------------------------------------------------------------------
	if (e[2]->GetNextOrig() != e[5] || e[2]->GetPrevOrig() != e[5] || e[2]->GetNextDest() != e[3] || e[2]->GetPrevDest() != e[3])
		return false;
	if (e[3]->GetNextOrig() != e[2] || e[3]->GetPrevOrig() != e[2] || e[3]->GetNextDest() != e[5] || e[3]->GetPrevDest() != e[5])
		return false;
	if (e[5]->GetNextOrig() != e[2] || e[5]->GetPrevOrig() != e[2] || e[5]->GetNextDest() != e[3] || e[5]->GetPrevDest() != e[3])
		return false;
	//-----------------------------------------------------------------------------
	if (e[2]->GetNextRight() != e[5] || e[2]->GetPrevRight() != e[3] || e[2]->GetNextLeft() != e[3] || e[2]->GetPrevLeft() != e[5])
		return false;
	if (e[3]->GetNextRight() != e[2] || e[3]->GetPrevRight() != e[5] || e[3]->GetNextLeft() != e[5] || e[3]->GetPrevLeft() != e[2])
		return false;
	if (e[5]->GetNextRight() != e[2] || e[5]->GetPrevRight() != e[3] || e[5]->GetNextLeft() != e[3] || e[5]->GetPrevLeft() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	//Copying test
	QuadEdge Tetra1;
	Tetra1.SetAsTetrahedron(v0.GetPoint(), v1.GetPoint(), v2.GetPoint(), v3.GetPoint());
	if (!Tetra1.TestIntegrity())
		return false;
	QuadEdge Tetra2(Tetra1), Tetra3;
	if (!Tetra2.TestIntegrity())
		return false;
	Tetra3 = Tetra2;
	if (!Tetra3.TestIntegrity())
		return false;
	Edge* eOrig[6];
	for (int i = 0; i < 6; ++i)
	{
		eOrig[i] = Tetra1.GetEdge(i);
		e[i] = Tetra3.GetEdge(i);
	}
	Face* fOrig[4];
	for (int i = 0; i < 4; ++i)
	{
		fOrig[i] = Tetra1.GetFace(i);
		f[i] = Tetra3.GetFace(i);
	}
	Vertex* vOrig[4];
	for (int i = 0; i < 4; ++i)
	{
		vOrig[i] = Tetra1.GetVertex(i);
		v[i] = Tetra3.GetVertex(i);
	}
	for (int i = 0; i < 6; ++i)
	{
		for (int j = 0; j < 6; ++j)
		{
			if (eOrig[i] == e[j])
				return false;
		}
	}
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (vOrig[i] == v[j])
				return false;
			if (fOrig[i] == f[j])
				return false;
		}
	}
	if (!tester_QuadEdge_Tetrahedron(Tetra3, v, f, e, points))
		return false;
	if (!tester_QuadEdge_1(NumTests))
		return false;
	if (!tester_QuadEdge_2(NumTests))
		return false;
	if (!tester_QuadEdge_Pyramid(NumTests))
		return false;
	if (!tester_QuadEdge_SetAsDoubleTriangulatedSquare(NumTests))
		return false;
	if (!tester_QuadEdge_3(NumTests))
		return false;
	if (!tester_QuadEdge_4(NumTests))
		return false;
	NumTests += 1;
	return true;
}
bool tester_QuadEdge_1(int& NumTests)
{
	//3 lines attached .--.--.--.
	Vector3D p[4];
	Vertex* v[4];
	QuadEdge g;
	for (int i = 0; i < 4; ++i)
	{
		p[i](0) = i;
		v[i] = g.AddVertex(p[i]);
	}
	Edge* e[3];
	Face* f = g.AddFace();
	for (int i = 0; i < 3; ++i)
	{
		e[i] = g.AddEdge();
		e[i]->SetOrig(v[i]);
		e[i]->SetDest(v[i + 1]);
		e[i]->SetRight(f);
		e[i]->SetLeft(f);
		v[i]->SetEdge(e[i]);
	}
	v[3]->SetEdge(e[2]);
	f->SetEdge(e[1]);
	for (int i = 0; i < 3; ++i)
	{
		e[i]->SetNextOrig(i == 0 ? e[i] : e[i - 1]);
		e[i]->SetNextDest(i == 2 ? e[i] : e[i + 1]);
		e[i]->SetNextRight(i == 2 ? e[i] : e[i + 1]);
		e[i]->SetNextLeft(i == 0 ? e[i] : e[i - 1]);
	}
	g.UpdatePrevs();
	if (!g.TestIntegrity())
		return false;
	for (int i = 0; i < 3; ++i)
	{
		if (e[i]->GetPrevOrig() != (i == 0 ? e[i] : e[i - 1]))
			return false;
		if (e[i]->GetPrevDest() != (i == 2 ? e[i] : e[i + 1]))
			return false;
		if (e[i]->GetPrevRight() != (i == 0 ? e[i] : e[i - 1]))
			return false;
		if (e[i]->GetPrevLeft() != (i == 2 ? e[i] : e[i + 1]))
			return false;
	}
	g.UpdateVerticesDegree();
	if (v[0]->GetDegree() != 1 || v[1]->GetDegree() != 2 || v[2]->GetDegree() != 2 || v[3]->GetDegree() != 1)
		return false;
	if (!g.isRemovable(e[0]) || g.isRemovable(e[1]) || !g.isRemovable(e[2]))
		return false;
	if (!g.IsConnected(v[0], v[1]) || !g.IsConnected(v[1], v[0]))
		return false;
	if (!g.IsConnected(v[1], v[2]) || !g.IsConnected(v[2], v[1]))
		return false;
	if (!g.IsConnected(v[2], v[3]) || !g.IsConnected(v[3], v[2]))
		return false;
	if (g.IsConnected(v[0], v[2]) || g.IsConnected(v[2], v[0]))
		return false;
	if (g.IsConnected(v[0], v[3]) || g.IsConnected(v[3], v[0]))
		return false;
	if (g.IsConnected(v[1], v[3]) || g.IsConnected(v[3], v[1]))
		return false;
	if (g.GetEdgeConnecting(v[0], v[1]) != e[0] || g.GetEdgeConnecting(v[1], v[0]) != e[0])
		return false;
	if (g.GetEdgeConnecting(v[1], v[2]) != e[1] || g.GetEdgeConnecting(v[2], v[1]) != e[1])
		return false;
	if (g.GetEdgeConnecting(v[2], v[3]) != e[2] || g.GetEdgeConnecting(v[3], v[2]) != e[2])
		return false;
	if (g.GetEdgeConnecting(v[0], v[2]) || g.GetEdgeConnecting(v[2], v[0]))
		return false;
	if (g.GetEdgeConnecting(v[0], v[3]) || g.GetEdgeConnecting(v[3], v[0]))
		return false;
	if (g.GetEdgeConnecting(v[1], v[3]) || g.GetEdgeConnecting(v[3], v[1]))
		return false;
	g.DeleteEdgeAndRightFace(e[2]);//e[1] nextRight is incorrect.
	if (!g.TestIntegrity())
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2])
		return false;
	if (e[0]->GetRight() != f || e[0]->GetLeft() != f || e[1]->GetRight() != f || e[1]->GetLeft() != f)
		return false;
	if (e[0]->GetNextOrig() != e[0] || e[0]->GetNextDest() != e[1] || e[1]->GetNextOrig() != e[0] || e[1]->GetNextDest() != e[1])
		return false;
	if (e[0]->GetPrevOrig() != e[0] || e[0]->GetPrevDest() != e[1] || e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[1])
		return false;
	if (e[0]->GetNextRight() != e[1] || e[0]->GetNextLeft() != e[0] || e[1]->GetNextRight() != e[1] /*|| e[1]->GetNextLeft() != e[0]*/)
		return false;
	if (!g.isRemovable(e[0]) || !g.isRemovable(e[1]))
		return false;
	if (!g.IsConnected(v[0], v[1]) || !g.IsConnected(v[1], v[0]))
		return false;
	if (!g.IsConnected(v[1], v[2]) || !g.IsConnected(v[2], v[1]))
		return false;
	if (g.IsConnected(v[0], v[2]) || g.IsConnected(v[2], v[0]))
		return false;
	if (g.GetEdgeConnecting(v[0], v[1]) != e[0] || g.GetEdgeConnecting(v[1], v[0]) != e[0])
		return false;
	if (g.GetEdgeConnecting(v[1], v[2]) != e[1] || g.GetEdgeConnecting(v[2], v[1]) != e[1])
		return false;
	if (g.GetEdgeConnecting(v[0], v[2]) || g.GetEdgeConnecting(v[2], v[0]))
		return false;
	g.DeleteEdgeAndRightFace(e[1]);
	if (!g.TestIntegrity())
		return false;
	if (g.isRemovable(e[0]))
		return false;
	NumTests += 1;
	return true;
}
bool tester_QuadEdge_Tetrahedron(QuadEdge& Tetrahedron, Vertex* v[4], Face* f[4], Edge* e[4], vector<Vector3D>& points)
{
	if (!Tetrahedron.TestIntegrity())
		return false;
	if (Tetrahedron.NumVertices() != 4 || Tetrahedron.NumEdges() != 6 || Tetrahedron.NumFaces() != 4)
		return false;	
	for (int i = 0; i < 6; ++i)
		e[i] = Tetrahedron.GetEdge(i);
	for (int i = 0; i < 4; ++i)
		f[i] = Tetrahedron.GetFace(i);
	for (int i = 0; i < 4; ++i)
		v[i] = Tetrahedron.GetVertex(i);
	for (int i = 0; i < 4; ++i)
	{
		if (v[i]->GetDegree() != 3)
			return false;
	}
	Tetrahedron.UpdateVerticesDegree();
	for (int i = 0; i < 4; ++i)
	{
		if (v[i]->GetDegree() != 3)
			return false;
	}
	if (v[0]->GetPoint() != points[0] || v[1]->GetPoint() != points[1] 
		|| v[2]->GetPoint() != points[2] || v[3]->GetPoint() != points[3])
		return false;
	vector<Edge*> edgeVerts[4];
	for (int i = 0; i < 4; ++i)
	{
		Tetrahedron.GetVertexEdges(v[i], edgeVerts[i]);
		if (edgeVerts[i].size() != 3)
			return false;
		if (v[i]->GetDegree() != 3)
			return false;
		if (Tetrahedron.UpdateVertexDegree(v[i]) != 3)
			return false;
	}
	vector<Edge*> edgeFaces[4];
	for (int i = 0; i < 4; ++i)
	{
		Tetrahedron.GetFaceEdges(Tetrahedron.GetFace(i), edgeFaces[i]);
		if (edgeFaces[i].size() != 3)
			return false;
	}
	for (int i = 0; i < 6; ++i)
		if (!Tetrahedron.isRemovable(e[i]))
			return false;
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2] || v[3]->GetEdge() != e[5])
		return false;
	if (f[0]->GetEdge() != e[0] || f[1]->GetEdge() != e[0] || f[2]->GetEdge() != e[1] || f[3]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[2])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[3])
		return false;
	if (e[3]->GetOrig() != v[0] || e[3]->GetDest() != v[3] || e[3]->GetRight() != f[1] || e[3]->GetLeft() != f[3])
		return false;
	if (e[4]->GetOrig() != v[1] || e[4]->GetDest() != v[3] || e[4]->GetRight() != f[2] || e[4]->GetLeft() != f[1])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[3] || e[0]->GetPrevOrig() != e[2] || e[0]->GetNextDest() != e[1] || e[0]->GetPrevDest() != e[4])
		return false;
	if (e[1]->GetNextOrig() != e[4] || e[1]->GetPrevOrig() != e[0] || e[1]->GetNextDest() != e[2] || e[1]->GetPrevDest() != e[5])
		return false;
	if (e[2]->GetNextOrig() != e[5] || e[2]->GetPrevOrig() != e[1] || e[2]->GetNextDest() != e[0] || e[2]->GetPrevDest() != e[3])
		return false;
	if (e[3]->GetNextOrig() != e[2] || e[3]->GetPrevOrig() != e[0] || e[3]->GetNextDest() != e[4] || e[3]->GetPrevDest() != e[5])
		return false;
	if (e[4]->GetNextOrig() != e[0] || e[4]->GetPrevOrig() != e[1] || e[4]->GetNextDest() != e[5] || e[4]->GetPrevDest() != e[3])
		return false;
	if (e[5]->GetNextOrig() != e[1] || e[5]->GetPrevOrig() != e[2] || e[5]->GetNextDest() != e[3] || e[5]->GetPrevDest() != e[4])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[2] || e[0]->GetPrevRight() != e[1] || e[0]->GetNextLeft() != e[4] || e[0]->GetPrevLeft() != e[3])
		return false;
	if (e[1]->GetNextRight() != e[0] || e[1]->GetPrevRight() != e[2] || e[1]->GetNextLeft() != e[5] || e[1]->GetPrevLeft() != e[4])
		return false;
	if (e[2]->GetNextRight() != e[1] || e[2]->GetPrevRight() != e[0] || e[2]->GetNextLeft() != e[3] || e[2]->GetPrevLeft() != e[5])
		return false;
	if (e[3]->GetNextRight() != e[0] || e[3]->GetPrevRight() != e[4] || e[3]->GetNextLeft() != e[5] || e[3]->GetPrevLeft() != e[2])
		return false;
	if (e[4]->GetNextRight() != e[1] || e[4]->GetPrevRight() != e[5] || e[4]->GetNextLeft() != e[3] || e[4]->GetPrevLeft() != e[0])
		return false;
	if (e[5]->GetNextRight() != e[2] || e[5]->GetPrevRight() != e[3] || e[5]->GetNextLeft() != e[4] || e[5]->GetPrevLeft() != e[1])
		return false;
	Vector3D Pv[4], Center;
	for (int i = 0; i < 4; ++i)
	{
		Pv[i] = v[i]->GetPoint();
		Center = Center + Pv[i] / 4.0;
	}
	if (Tetrahedron.isVisible(f[0], Pv[3]))
		return false;
	if (Tetrahedron.isVisible(f[1], Pv[2]))
		return false;
	if (Tetrahedron.isVisible(f[2], Pv[0]))
		return false;
	if (Tetrahedron.isVisible(f[3], Pv[1]))
		return false;
	if (Tetrahedron.isVisible(f[0], Center))
		return false;
	if (Tetrahedron.isVisible(f[1], Center))
		return false;
	if (Tetrahedron.isVisible(f[2], Center))
		return false;
	if (Tetrahedron.isVisible(f[3], Center))
		return false;
	Vector3D P4(0, 0, 2);
	if (Tetrahedron.isVisible(f[0], P4))
		return false;
	if (!Tetrahedron.isVisible(f[1], P4))
		return false;
	if (!Tetrahedron.isVisible(f[2], P4))
		return false;
	if (!Tetrahedron.isVisible(f[3], P4))
		return false;
	return true;
}
bool tester_QuadEdge_2(int& NumTests)
{
	vector<Vector3D> points;
	points.resize(5, 0);
	points[1](0) = 1;
	points[2](1) = 1;
	points[3](0) = points[3](1) = points[3](2) = 1.0 / 3.0;
	Vertex* v[4];
	Face* f[4];
	Edge* e[6];
	QuadEdge T;
	T.SetAsTetrahedron(points[0], points[1], points[2], points[3]);
	if (!T.TestIntegrity())
		return false;
	if (!tester_QuadEdge_Tetrahedron(T, v, f, e, points))
		return false;
	T.DeleteEdgeAndRightFace(e[3]);
	if (!T.TestIntegrity())
		return false;
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2] || v[3]->GetEdge() != e[5])
		return false;
	if (f[0]->GetEdge() != e[0] || f[2]->GetEdge() != e[1] || f[3]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[3])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[2])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[3])
		return false;
	if (e[4]->GetOrig() != v[1] || e[4]->GetDest() != v[3] || e[4]->GetRight() != f[2] || e[4]->GetLeft() != f[3])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[2] || e[0]->GetPrevOrig() != e[2] || e[0]->GetNextDest() != e[1] || e[0]->GetPrevDest() != e[4])
		return false;
	if (e[1]->GetNextOrig() != e[4] || e[1]->GetPrevOrig() != e[0] || e[1]->GetNextDest() != e[2] || e[1]->GetPrevDest() != e[5])
		return false;
	if (e[2]->GetNextOrig() != e[5] || e[2]->GetPrevOrig() != e[1] || e[2]->GetNextDest() != e[0] || e[2]->GetPrevDest() != e[0])
		return false;
	if (e[4]->GetNextOrig() != e[0] || e[4]->GetPrevOrig() != e[1] || e[4]->GetNextDest() != e[5] || e[4]->GetPrevDest() != e[5])
		return false;
	if (e[5]->GetNextOrig() != e[1] || e[5]->GetPrevOrig() != e[2] || e[5]->GetNextDest() != e[4] || e[5]->GetPrevDest() != e[4])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[2] || e[0]->GetPrevRight() != e[1] || e[0]->GetNextLeft() != e[4] || e[0]->GetPrevLeft() != e[2])
		return false;
	if (e[1]->GetNextRight() != e[0] || e[1]->GetPrevRight() != e[2] || e[1]->GetNextLeft() != e[5] || e[1]->GetPrevLeft() != e[4])
		return false;
	if (e[2]->GetNextRight() != e[1] || e[2]->GetPrevRight() != e[0] || e[2]->GetNextLeft() != e[0] || e[2]->GetPrevLeft() != e[5])
		return false;
	if (e[4]->GetNextRight() != e[1] || e[4]->GetPrevRight() != e[5] || e[4]->GetNextLeft() != e[5] || e[4]->GetPrevLeft() != e[0])
		return false;
	if (e[5]->GetNextRight() != e[2] || e[5]->GetPrevRight() != e[4] || e[5]->GetNextLeft() != e[4] || e[5]->GetPrevLeft() != e[1])
		return false;
	T.DeleteEdgeAndRightFace(e[4]);
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2] || v[3]->GetEdge() != e[5])
		return false;
	if (f[0]->GetEdge() != e[0] || f[3]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[3])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[3])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[3])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetRight() != f[3] || e[5]->GetLeft() != f[3])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextOrig() != e[2] || e[0]->GetPrevOrig() != e[2] || e[0]->GetNextDest() != e[1] || e[0]->GetPrevDest() != e[1])
		return false;
	if (e[1]->GetNextOrig() != e[0] || e[1]->GetPrevOrig() != e[0] || e[1]->GetNextDest() != e[2] || e[1]->GetPrevDest() != e[5])
		return false;
	if (e[2]->GetNextOrig() != e[5] || e[2]->GetPrevOrig() != e[1] || e[2]->GetNextDest() != e[0] || e[2]->GetPrevDest() != e[0])
		return false;
	if (e[5]->GetNextOrig() != e[1] || e[5]->GetPrevOrig() != e[2] || e[5]->GetNextDest() != e[5] || e[5]->GetPrevDest() != e[5])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetNextRight() != e[2] || e[0]->GetPrevRight() != e[1] || e[0]->GetNextLeft() != e[1] || e[0]->GetPrevLeft() != e[2])
		return false;
	if (e[1]->GetNextRight() != e[0] || e[1]->GetPrevRight() != e[2] || e[1]->GetNextLeft() != e[5] || e[1]->GetPrevLeft() != e[0])
		return false;
	if (e[2]->GetNextRight() != e[1] || e[2]->GetPrevRight() != e[0] || e[2]->GetNextLeft() != e[0] || e[2]->GetPrevLeft() != e[5])
		return false;
	if (e[5]->GetNextRight() != e[2] || e[5]->GetPrevRight() != e[5] || e[5]->GetNextLeft() != e[5] || e[5]->GetPrevLeft() != e[1])
		return false;
	T.DeleteEdgeAndRightFace(e[5]);
	if (!T.TestIntegrity())
		return false;
	if (v[0]->GetEdge() != e[0] || v[1]->GetEdge() != e[1] || v[2]->GetEdge() != e[2])
		return false;
	if (f[0]->GetEdge() != e[0] || f[3]->GetEdge() != e[2])
		return false;
	//-----------------------------------------------------------------------------
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetRight() != f[0] || e[0]->GetLeft() != f[3])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetRight() != f[0] || e[1]->GetLeft() != f[3])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetRight() != f[0] || e[2]->GetLeft() != f[3])
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
bool tester_QuadEdge_Pyramid(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(5);
	input[1](0) = 1;
	input[2](1) = 1;
	input[3](0) = input[3](1) = 1;
	input[4](0) = input[4](1) = input[4](2) = 0.5;
	QuadEdge hull;
	hull.SetAsTetrahedron(input[0], input[1], input[2], input[4]);
	if (!hull.TestIntegrity())
		return false;
	if (hull.NumVertices() != 4 || hull.NumEdges() != 6 || hull.NumFaces() != 4)
		return false;
	EdgeContainer* E1 = hull.GetEdgeList()->GetContainer(1);
	Face* f2 = hull.GetFace(2);
	hull.DeleteEdgeAndRightFace(E1);
	if (!hull.TestIntegrity())
		return false;
	hull.AddPyramidAndDeleteFace(input[3], f2);
	if (!hull.TestIntegrity())
		return false;
	if (hull.NumVertices() != 5 || hull.NumEdges() != 9 || hull.NumFaces() != 6)
		return false;
	Vertex* v[5];
	Face* f[6];
	Edge* e[9];
	EdgeContainer* E[9];
	for (int i = 0; i < 5; ++i)
	{
		v[i] = hull.GetVertex(i);
	}
	for (int i = 0; i < 6; ++i)
	{
		f[i] = hull.GetFace(i);
	}
	for (int i = 0; i < 9; ++i)
	{
		e[i] = hull.GetEdge(i);
		E[i] = hull.GetEdgeList()->GetContainer(i);
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
	NumTests += 1;
	return true;
}
bool tester_QuadEdge_SetAsDoubleTriangulatedSquare(int& NumTests)
{
	Vector3D pA, pB(0, 1), pC(1, 1), pD(1, 0);
	QuadEdge Q;
	Q.SetAsDoubleTriangulatedSquare(pA, pB, pC, pD);
	if (!Q.TestIntegrity())
		return false;
	if (Q.NumVertices() != 4 || Q.NumEdges() != 6 || Q.NumFaces() != 4)
		return false;
	Vertex* A = Q.GetVertex(0);
	Vertex* B = Q.GetVertex(1);
	Vertex* C = Q.GetVertex(2);
	Vertex* D = Q.GetVertex(3);
	Face* ABC = Q.GetFace(0);
	Face* ACB = Q.GetFace(1);
	Face* ACD = Q.GetFace(2);
	Face* ADC = Q.GetFace(3);
	Edge* eAB = Q.GetEdge(0);
	Edge* eBC = Q.GetEdge(1);
	Edge* eCA = Q.GetEdge(2);
	Edge* eAC = Q.GetEdge(3);
	Edge* eCD = Q.GetEdge(4);
	Edge* eDA = Q.GetEdge(5);
	if (A->GetDegree() != 4 || B->GetDegree() != 2 || C->GetDegree() != 4 || D->GetDegree() != 2)
		return false;
	EdgeIterator itA(A);
	if (itA.Next() != eAB || itA.Next() != eAC || itA.Next() != eDA || itA.Next() != eCA || itA.Next())
		return false;
	EdgeIterator itB(B);
	if (itB.Next() != eBC || itB.Next() != eAB || itB.Next())
		return false;
	EdgeIterator itC(C);
	if (itC.Next() != eCA || itC.Next() != eAC || itC.Next() != eCD || itC.Next() != eBC || itC.Next())
		return false;
	EdgeIterator itD(D);
	if (itD.Next() != eCD || itD.Next() != eDA || itD.Next())
		return false;
	VertexIterator itADC(ADC);
	if (itADC.Next() != C || itADC.Next() != A || itADC.Next() != D || itADC.Next())
		return false;
	VertexIterator itACD(ACD);
	if (itACD.Next() != C || itACD.Next() != D || itACD.Next() != A || itACD.Next())
		return false;
	VertexIterator itABC(ABC);
	if (itABC.Next() != C || itABC.Next() != A || itABC.Next() != B || itABC.Next())
		return false;
	VertexIterator itACB(ACB);
	if (itACB.Next() != C || itACB.Next() != B || itACB.Next() != A || itACB.Next())
		return false;
	EdgeIterator iteADC(ADC);
	if (iteADC.Next() != eAC || iteADC.Next() != eDA || iteADC.Next() != eCD || iteADC.Next())
		return false;
	EdgeIterator iteACD(ACD);
	if (iteACD.Next() != eCD || iteACD.Next() != eDA || iteACD.Next() != eCA || iteACD.Next())
		return false;
	EdgeIterator iteABC(ABC);
	if (iteABC.Next() != eCA || iteABC.Next() != eAB || iteABC.Next() != eBC || iteABC.Next())
		return false;
	EdgeIterator iteACB(ACB);
	if (iteACB.Next() != eBC || iteACB.Next() != eAB || iteACB.Next() != eAC || iteACB.Next())
		return false;
	NumTests += 1;
	return true;
}
bool tester_QuadEdge_3(int& NumTests)
{
	//Tetrahedron
	VertexBasic v0, v1(1), v2(0, 1), v3(0, 0, 1); 
	QuadEdge Tetrahedron;
	Tetrahedron.SetAsTetrahedron(v0.GetPoint(), v1.GetPoint(), v2.GetPoint(), v3.GetPoint());//A,B,C,D
	if (!Tetrahedron.TestIntegrity())
		return false;
	Edge* e[7];//eAB, eBC, eCA
	Face* f[4];
	Vertex* v[5];
	vector<Vector3D> points;
	points.push_back(v0.GetPoint());
	points.push_back(v1.GetPoint());
	points.push_back(v2.GetPoint());
	points.push_back(v3.GetPoint());
	if (!tester_QuadEdge_Tetrahedron(Tetrahedron, v, f, e, points))
		return false;
	v[4] = Tetrahedron.InsertVertexAtMid(e[4]);
	e[6] = v[4]->GetEdge();
	Tetrahedron.IssueIndices();
	Tetrahedron.print("..\\Data\\tester_QuadEdge_3.txt");
	Tetrahedron.PrintTestScript("..\\Data\\tester_QuadEdge_3.txt");
	if (Tetrahedron.NumVertices() != 5 || Tetrahedron.NumEdges() != 7 || Tetrahedron.NumFaces() != 4)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 3)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[2] || v[2]->GetDegree() != 3)
		return false;
	if (v[3]->GetEdge() != e[5] || v[3]->GetDegree() != 3)
		return false;
	if (v[4]->GetEdge() != e[6] || v[4]->GetDegree() != 2)
		return false;
	if (f[0]->GetEdge() != e[0] || f[0]->GetDegree() != 3)
		return false;
	if (f[1]->GetEdge() != e[0] || f[1]->GetDegree() != 4)
		return false;
	if (f[2]->GetEdge() != e[1] || f[2]->GetDegree() != 4)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[1] || e[0]->GetRight() != f[0])
		return false;
	if (e[0]->GetNextOrig() != e[3] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[4] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[4] || e[0]->GetPrevLeft() != e[3] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[2] || e[1]->GetRight() != f[0])
		return false;
	if (e[1]->GetNextOrig() != e[4] || e[1]->GetNextDest() != e[2] || e[1]->GetNextLeft() != e[5] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[5] || e[1]->GetPrevLeft() != e[4] || e[1]->GetPrevRight() != e[2])
		return false;
	if (e[2]->GetOrig() != v[2] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[0])
		return false;
	if (e[2]->GetNextOrig() != e[5] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[3] || e[2]->GetNextRight() != e[1])
		return false;
	if (e[2]->GetPrevOrig() != e[1] || e[2]->GetPrevDest() != e[3] || e[2]->GetPrevLeft() != e[5] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[3]->GetOrig() != v[0] || e[3]->GetDest() != v[3] || e[3]->GetLeft() != f[3] || e[3]->GetRight() != f[1])
		return false;
	if (e[3]->GetNextOrig() != e[2] || e[3]->GetNextDest() != e[6] || e[3]->GetNextLeft() != e[5] || e[3]->GetNextRight() != e[0])
		return false;
	if (e[3]->GetPrevOrig() != e[0] || e[3]->GetPrevDest() != e[5] || e[3]->GetPrevLeft() != e[2] || e[3]->GetPrevRight() != e[6])
		return false;
	if (e[4]->GetOrig() != v[1] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[1] || e[4]->GetRight() != f[2])
		return false;
	if (e[4]->GetNextOrig() != e[0] || e[4]->GetNextDest() != e[6] || e[4]->GetNextLeft() != e[6] || e[4]->GetNextRight() != e[1])
		return false;
	if (e[4]->GetPrevOrig() != e[1] || e[4]->GetPrevDest() != e[6] || e[4]->GetPrevLeft() != e[0] || e[4]->GetPrevRight() != e[6])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[2] || e[5]->GetRight() != f[3])
		return false;
	if (e[5]->GetNextOrig() != e[1] || e[5]->GetNextDest() != e[3] || e[5]->GetNextLeft() != e[6] || e[5]->GetNextRight() != e[2])
		return false;
	if (e[5]->GetPrevOrig() != e[2] || e[5]->GetPrevDest() != e[6] || e[5]->GetPrevLeft() != e[1] || e[5]->GetPrevRight() != e[3])
		return false;
	if (e[6]->GetOrig() != v[4] || e[6]->GetDest() != v[3] || e[6]->GetLeft() != f[1] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[4] || e[6]->GetNextDest() != e[5] || e[6]->GetNextLeft() != e[3] || e[6]->GetNextRight() != e[4])
		return false;
	if (e[6]->GetPrevOrig() != e[4] || e[6]->GetPrevDest() != e[3] || e[6]->GetPrevLeft() != e[4] || e[6]->GetPrevRight() != e[5])
		return false;
	NumTests += 1;
	return true;
}
bool tester_QuadEdge_4(int& NumTests)
{//test write to and read from a file
	Vector3D v0, v1(1.1), v2(0, 1.2), v3(0, 0, 1.3584); 
	QuadEdge Tetra1;
	Tetra1.SetAsTetrahedron(v0, v1, v2, v3);
	stringstream qe_stream;
	qe_stream << Tetra1;
	QuadEdge Tetra2;
	qe_stream >> Tetra2;
	Tetra2.IssueIndices();
	Tetra2.print("..\\Data\\Tetra2.txt");
	//---------------------------------------------------------------------------------------
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
	QuadEdge* qe = T.GetMesh2D();
	stringstream qe_stream2;
	qe_stream2 << (*qe);
	QuadEdge* qe2 = new QuadEdge;
	qe_stream2 >> (*qe2);
	Vertex* v[20];
	for (int i = 0; i < 20; ++i)
	{
		v[i] = qe2->GetVertex(i);
	}
	Face* f[28];
	for (int i = 0; i < 28; ++i)
	{
		f[i] = qe2->GetFace(i);
	}
	Edge* e[46];
	for (int i = 0; i < 46; ++i)
	{
		e[i] = qe2->GetEdge(i);
	}
	qe2->IssueIndices();
	qe2->UpdateVerticesDegree();
	qe2->UpdateFacesDegree();
	qe2->print("..\\Data\\RS2_10x10.txt");
	qe2->PrintTestScript("..\\Data\\RS2_10x10.testScript.txt");
	Triangulation T2;
	T2 << (*qe2);
	if (!T2.TestIntegrity())
		return false;
	if (!T2.TestDelaunay())
		return false;
	double alpha_max = T2.GetMaxAngle() / 3.14159 * 180.0;
	if (fabs(alpha_max - 120.9729) > 0.01)
		return false;
	int Nv2 = T2.NumVertices();
	if (Nv2 != 20)
		return false;
	if (T.GetMesh2D()->NumVertices() != 20 || T.GetMesh2D()->NumEdges() != 46 || T.GetMesh2D()->NumFaces() != 28)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 3)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 4)
		return false;
	if (v[2]->GetEdge() != e[21] || v[2]->GetDegree() != 3)
		return false;
	if (v[3]->GetEdge() != e[2] || v[3]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[3] || v[4]->GetDegree() != 4)
		return false;
	if (v[5]->GetEdge() != e[4] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[5] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[6] || v[7]->GetDegree() != 4)
		return false;
	if (v[8]->GetEdge() != e[7] || v[8]->GetDegree() != 4)
		return false;
	if (v[9]->GetEdge() != e[8] || v[9]->GetDegree() != 4)
		return false;
	if (v[10]->GetEdge() != e[10] || v[10]->GetDegree() != 3)
		return false;
	if (v[11]->GetEdge() != e[24] || v[11]->GetDegree() != 5)
		return false;
	if (v[12]->GetEdge() != e[13] || v[12]->GetDegree() != 4)
		return false;
	if (v[13]->GetEdge() != e[23] || v[13]->GetDegree() != 5)
		return false;
	if (v[14]->GetEdge() != e[37] || v[14]->GetDegree() != 6)
		return false;
	if (v[15]->GetEdge() != e[27] || v[15]->GetDegree() != 5)
		return false;
	if (v[16]->GetEdge() != e[27] || v[16]->GetDegree() != 7)
		return false;
	if (v[17]->GetEdge() != e[28] || v[17]->GetDegree() != 7)
		return false;
	if (v[18]->GetEdge() != e[40] || v[18]->GetDegree() != 7)
		return false;
	if (v[19]->GetEdge() != e[45] || v[19]->GetDegree() != 5)
		return false;
	if (f[0]->GetEdge() != e[10] || f[0]->GetDegree() != 11)
		return false;
	if (f[1]->GetEdge() != e[4] || f[1]->GetDegree() != 3)
		return false;
	if (f[2]->GetEdge() != e[7] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[14] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[6] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[2] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[1] || f[6]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[15] || f[7]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[12] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[3] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[20] || f[10]->GetDegree() != 3)
		return false;
	if (f[11]->GetEdge() != e[27] || f[11]->GetDegree() != 3)
		return false;
	if (f[12]->GetEdge() != e[22] || f[12]->GetDegree() != 3)
		return false;
	if (f[13]->GetEdge() != e[0] || f[13]->GetDegree() != 3)
		return false;
	if (f[14]->GetEdge() != e[9] || f[14]->GetDegree() != 3)
		return false;
	if (f[15]->GetEdge() != e[10] || f[15]->GetDegree() != 3)
		return false;
	if (f[16]->GetEdge() != e[8] || f[16]->GetDegree() != 3)
		return false;
	if (f[17]->GetEdge() != e[13] || f[17]->GetDegree() != 3)
		return false;
	if (f[18]->GetEdge() != e[17] || f[18]->GetDegree() != 3)
		return false;
	if (f[19]->GetEdge() != e[16] || f[19]->GetDegree() != 3)
		return false;
	if (f[20]->GetEdge() != e[23] || f[20]->GetDegree() != 3)
		return false;
	if (f[21]->GetEdge() != e[28] || f[21]->GetDegree() != 3)
		return false;
	if (f[22]->GetEdge() != e[33] || f[22]->GetDegree() != 3)
		return false;
	if (f[23]->GetEdge() != e[5] || f[23]->GetDegree() != 3)
		return false;
	if (f[24]->GetEdge() != e[11] || f[24]->GetDegree() != 3)
		return false;
	if (f[25]->GetEdge() != e[15] || f[25]->GetDegree() != 3)
		return false;
	if (f[26]->GetEdge() != e[16] || f[26]->GetDegree() != 3)
		return false;
	if (f[27]->GetEdge() != e[19] || f[27]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[1] || e[0]->GetDest() != v[0] || e[0]->GetLeft() != f[13] || e[0]->GetRight() != f[0])
		return false;
	if (e[0]->GetNextOrig() != e[1] || e[0]->GetNextDest() != e[31] || e[0]->GetNextLeft() != e[31] || e[0]->GetNextRight() != e[1])
		return false;
	if (e[0]->GetPrevOrig() != e[22] || e[0]->GetPrevDest() != e[9] || e[0]->GetPrevLeft() != e[30] || e[0]->GetPrevRight() != e[9])
		return false;
	if (e[1]->GetOrig() != v[2] || e[1]->GetDest() != v[1] || e[1]->GetLeft() != f[6] || e[1]->GetRight() != f[0])
		return false;
	if (e[1]->GetNextOrig() != e[21] || e[1]->GetNextDest() != e[30] || e[1]->GetNextLeft() != e[22] || e[1]->GetNextRight() != e[2])
		return false;
	if (e[1]->GetPrevOrig() != e[2] || e[1]->GetPrevDest() != e[0] || e[1]->GetPrevLeft() != e[21] || e[1]->GetPrevRight() != e[0])
		return false;
	if (e[2]->GetOrig() != v[3] || e[2]->GetDest() != v[2] || e[2]->GetLeft() != f[5] || e[2]->GetRight() != f[0])
		return false;
	if (e[2]->GetNextOrig() != e[26] || e[2]->GetNextDest() != e[1] || e[2]->GetNextLeft() != e[21] || e[2]->GetNextRight() != e[3])
		return false;
	if (e[2]->GetPrevOrig() != e[3] || e[2]->GetPrevDest() != e[21] || e[2]->GetPrevLeft() != e[20] || e[2]->GetPrevRight() != e[1])
		return false;
	if (e[3]->GetOrig() != v[4] || e[3]->GetDest() != v[3] || e[3]->GetLeft() != f[9] || e[3]->GetRight() != f[0])
		return false;
	if (e[3]->GetNextOrig() != e[25] || e[3]->GetNextDest() != e[2] || e[3]->GetNextLeft() != e[26] || e[3]->GetNextRight() != e[4])
		return false;
	if (e[3]->GetPrevOrig() != e[4] || e[3]->GetPrevDest() != e[20] || e[3]->GetPrevLeft() != e[25] || e[3]->GetPrevRight() != e[2])
		return false;
	if (e[4]->GetOrig() != v[5] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[1] || e[4]->GetRight() != f[0])
		return false;
	if (e[4]->GetNextOrig() != e[42] || e[4]->GetNextDest() != e[3] || e[4]->GetNextLeft() != e[12] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[5] || e[4]->GetPrevDest() != e[12] || e[4]->GetPrevLeft() != e[11] || e[4]->GetPrevRight() != e[3])
		return false;
	if (e[5]->GetOrig() != v[6] || e[5]->GetDest() != v[5] || e[5]->GetLeft() != f[23] || e[5]->GetRight() != f[0])
		return false;
	if (e[5]->GetNextOrig() != e[41] || e[5]->GetNextDest() != e[4] || e[5]->GetNextLeft() != e[42] || e[5]->GetNextRight() != e[6])
		return false;
	if (e[5]->GetPrevOrig() != e[6] || e[5]->GetPrevDest() != e[11] || e[5]->GetPrevLeft() != e[41] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[7] || e[6]->GetDest() != v[6] || e[6]->GetLeft() != f[4] || e[6]->GetRight() != f[0])
		return false;
	if (e[6]->GetNextOrig() != e[18] || e[6]->GetNextDest() != e[5] || e[6]->GetNextLeft() != e[19] || e[6]->GetNextRight() != e[7])
		return false;
	if (e[6]->GetPrevOrig() != e[7] || e[6]->GetPrevDest() != e[19] || e[6]->GetPrevLeft() != e[18] || e[6]->GetPrevRight() != e[5])
		return false;
	if (e[7]->GetOrig() != v[8] || e[7]->GetDest() != v[7] || e[7]->GetLeft() != f[2] || e[7]->GetRight() != f[0])
		return false;
	if (e[7]->GetNextOrig() != e[35] || e[7]->GetNextDest() != e[6] || e[7]->GetNextLeft() != e[14] || e[7]->GetNextRight() != e[8])
		return false;
	if (e[7]->GetPrevOrig() != e[8] || e[7]->GetPrevDest() != e[14] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[6])
		return false;
	if (e[8]->GetOrig() != v[9] || e[8]->GetDest() != v[8] || e[8]->GetLeft() != f[16] || e[8]->GetRight() != f[0])
		return false;
	if (e[8]->GetNextOrig() != e[34] || e[8]->GetNextDest() != e[7] || e[8]->GetNextLeft() != e[35] || e[8]->GetNextRight() != e[10])
		return false;
	if (e[8]->GetPrevOrig() != e[10] || e[8]->GetPrevDest() != e[13] || e[8]->GetPrevLeft() != e[34] || e[8]->GetPrevRight() != e[7])
		return false;
	if (e[9]->GetOrig() != v[0] || e[9]->GetDest() != v[10] || e[9]->GetLeft() != f[14] || e[9]->GetRight() != f[0])
		return false;
	if (e[9]->GetNextOrig() != e[0] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[32] || e[9]->GetNextRight() != e[0])
		return false;
	if (e[9]->GetPrevOrig() != e[31] || e[9]->GetPrevDest() != e[32] || e[9]->GetPrevLeft() != e[31] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[10] || e[10]->GetDest() != v[9] || e[10]->GetLeft() != f[15] || e[10]->GetRight() != f[0])
		return false;
	if (e[10]->GetNextOrig() != e[32] || e[10]->GetNextDest() != e[8] || e[10]->GetNextLeft() != e[33] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[33] || e[10]->GetPrevLeft() != e[32] || e[10]->GetPrevRight() != e[8])
		return false;
	if (e[11]->GetOrig() != v[5] || e[11]->GetDest() != v[11] || e[11]->GetLeft() != f[24] || e[11]->GetRight() != f[1])
		return false;
	if (e[11]->GetNextOrig() != e[5] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[43] || e[11]->GetNextRight() != e[4])
		return false;
	if (e[11]->GetPrevOrig() != e[42] || e[11]->GetPrevDest() != e[15] || e[11]->GetPrevLeft() != e[42] || e[11]->GetPrevRight() != e[12])
		return false;
	if (e[12]->GetOrig() != v[11] || e[12]->GetDest() != v[4] || e[12]->GetLeft() != f[8] || e[12]->GetRight() != f[1])
		return false;
	if (e[12]->GetNextOrig() != e[24] || e[12]->GetNextDest() != e[4] || e[12]->GetNextLeft() != e[25] || e[12]->GetNextRight() != e[11])
		return false;
	if (e[12]->GetPrevOrig() != e[11] || e[12]->GetPrevDest() != e[25] || e[12]->GetPrevLeft() != e[24] || e[12]->GetPrevRight() != e[4])
		return false;
	if (e[13]->GetOrig() != v[8] || e[13]->GetDest() != v[12] || e[13]->GetLeft() != f[17] || e[13]->GetRight() != f[2])
		return false;
	if (e[13]->GetNextOrig() != e[8] || e[13]->GetNextDest() != e[36] || e[13]->GetNextLeft() != e[36] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[35] || e[13]->GetPrevDest() != e[14] || e[13]->GetPrevLeft() != e[35] || e[13]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[12] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[3] || e[14]->GetRight() != f[2])
		return false;
	if (e[14]->GetNextOrig() != e[13] || e[14]->GetNextDest() != e[7] || e[14]->GetNextLeft() != e[18] || e[14]->GetNextRight() != e[13])
		return false;
	if (e[14]->GetPrevOrig() != e[17] || e[14]->GetPrevDest() != e[18] || e[14]->GetPrevLeft() != e[17] || e[14]->GetPrevRight() != e[7])
		return false;
	if (e[15]->GetOrig() != v[11] || e[15]->GetDest() != v[13] || e[15]->GetLeft() != f[25] || e[15]->GetRight() != f[7])
		return false;
	if (e[15]->GetNextOrig() != e[11] || e[15]->GetNextDest() != e[23] || e[15]->GetNextLeft() != e[44] || e[15]->GetNextRight() != e[24])
		return false;
	if (e[15]->GetPrevOrig() != e[43] || e[15]->GetPrevDest() != e[16] || e[15]->GetPrevLeft() != e[43] || e[15]->GetPrevRight() != e[23])
		return false;
	if (e[16]->GetOrig() != v[13] || e[16]->GetDest() != v[14] || e[16]->GetLeft() != f[26] || e[16]->GetRight() != f[19])
		return false;
	if (e[16]->GetNextOrig() != e[15] || e[16]->GetNextDest() != e[17] || e[16]->GetNextLeft() != e[45] || e[16]->GetNextRight() != e[38])
		return false;
	if (e[16]->GetPrevOrig() != e[38] || e[16]->GetPrevDest() != e[45] || e[16]->GetPrevLeft() != e[44] || e[16]->GetPrevRight() != e[37])
		return false;
	if (e[17]->GetOrig() != v[12] || e[17]->GetDest() != v[14] || e[17]->GetLeft() != f[18] || e[17]->GetRight() != f[3])
		return false;
	if (e[17]->GetNextOrig() != e[14] || e[17]->GetNextDest() != e[18] || e[17]->GetNextLeft() != e[37] || e[17]->GetNextRight() != e[14])
		return false;
	if (e[17]->GetPrevOrig() != e[36] || e[17]->GetPrevDest() != e[16] || e[17]->GetPrevLeft() != e[36] || e[17]->GetPrevRight() != e[18])
		return false;
	if (e[18]->GetOrig() != v[7] || e[18]->GetDest() != v[14] || e[18]->GetLeft() != f[3] || e[18]->GetRight() != f[4])
		return false;
	if (e[18]->GetNextOrig() != e[14] || e[18]->GetNextDest() != e[19] || e[18]->GetNextLeft() != e[17] || e[18]->GetNextRight() != e[6])
		return false;
	if (e[18]->GetPrevOrig() != e[6] || e[18]->GetPrevDest() != e[17] || e[18]->GetPrevLeft() != e[14] || e[18]->GetPrevRight() != e[19])
		return false;
	if (e[19]->GetOrig() != v[14] || e[19]->GetDest() != v[6] || e[19]->GetLeft() != f[27] || e[19]->GetRight() != f[4])
		return false;
	if (e[19]->GetNextOrig() != e[37] || e[19]->GetNextDest() != e[6] || e[19]->GetNextLeft() != e[41] || e[19]->GetNextRight() != e[18])
		return false;
	if (e[19]->GetPrevOrig() != e[18] || e[19]->GetPrevDest() != e[41] || e[19]->GetPrevLeft() != e[45] || e[19]->GetPrevRight() != e[6])
		return false;
	if (e[20]->GetOrig() != v[3] || e[20]->GetDest() != v[15] || e[20]->GetLeft() != f[10] || e[20]->GetRight() != f[5])
		return false;
	if (e[20]->GetNextOrig() != e[3] || e[20]->GetNextDest() != e[21] || e[20]->GetNextLeft() != e[27] || e[20]->GetNextRight() != e[2])
		return false;
	if (e[20]->GetPrevOrig() != e[26] || e[20]->GetPrevDest() != e[29] || e[20]->GetPrevLeft() != e[26] || e[20]->GetPrevRight() != e[21])
		return false;
	if (e[21]->GetOrig() != v[2] || e[21]->GetDest() != v[15] || e[21]->GetLeft() != f[5] || e[21]->GetRight() != f[6])
		return false;
	if (e[21]->GetNextOrig() != e[2] || e[21]->GetNextDest() != e[22] || e[21]->GetNextLeft() != e[20] || e[21]->GetNextRight() != e[1])
		return false;
	if (e[21]->GetPrevOrig() != e[1] || e[21]->GetPrevDest() != e[20] || e[21]->GetPrevLeft() != e[2] || e[21]->GetPrevRight() != e[22])
		return false;
	if (e[22]->GetOrig() != v[15] || e[22]->GetDest() != v[1] || e[22]->GetLeft() != f[12] || e[22]->GetRight() != f[6])
		return false;
	if (e[22]->GetNextOrig() != e[27] || e[22]->GetNextDest() != e[0] || e[22]->GetNextLeft() != e[30] || e[22]->GetNextRight() != e[21])
		return false;
	if (e[22]->GetPrevOrig() != e[21] || e[22]->GetPrevDest() != e[30] || e[22]->GetPrevLeft() != e[29] || e[22]->GetPrevRight() != e[1])
		return false;
	if (e[23]->GetOrig() != v[13] || e[23]->GetDest() != v[16] || e[23]->GetLeft() != f[20] || e[23]->GetRight() != f[7])
		return false;
	if (e[23]->GetNextOrig() != e[44] || e[23]->GetNextDest() != e[24] || e[23]->GetNextLeft() != e[39] || e[23]->GetNextRight() != e[15])
		return false;
	if (e[23]->GetPrevOrig() != e[15] || e[23]->GetPrevDest() != e[28] || e[23]->GetPrevLeft() != e[38] || e[23]->GetPrevRight() != e[24])
		return false;
	if (e[24]->GetOrig() != v[11] || e[24]->GetDest() != v[16] || e[24]->GetLeft() != f[7] || e[24]->GetRight() != f[8])
		return false;
	if (e[24]->GetNextOrig() != e[43] || e[24]->GetNextDest() != e[25] || e[24]->GetNextLeft() != e[23] || e[24]->GetNextRight() != e[12])
		return false;
	if (e[24]->GetPrevOrig() != e[12] || e[24]->GetPrevDest() != e[23] || e[24]->GetPrevLeft() != e[15] || e[24]->GetPrevRight() != e[25])
		return false;
	if (e[25]->GetOrig() != v[4] || e[25]->GetDest() != v[16] || e[25]->GetLeft() != f[8] || e[25]->GetRight() != f[9])
		return false;
	if (e[25]->GetNextOrig() != e[12] || e[25]->GetNextDest() != e[26] || e[25]->GetNextLeft() != e[24] || e[25]->GetNextRight() != e[3])
		return false;
	if (e[25]->GetPrevOrig() != e[3] || e[25]->GetPrevDest() != e[24] || e[25]->GetPrevLeft() != e[12] || e[25]->GetPrevRight() != e[26])
		return false;
	if (e[26]->GetOrig() != v[3] || e[26]->GetDest() != v[16] || e[26]->GetLeft() != f[9] || e[26]->GetRight() != f[10])
		return false;
	if (e[26]->GetNextOrig() != e[20] || e[26]->GetNextDest() != e[27] || e[26]->GetNextLeft() != e[25] || e[26]->GetNextRight() != e[20])
		return false;
	if (e[26]->GetPrevOrig() != e[2] || e[26]->GetPrevDest() != e[25] || e[26]->GetPrevLeft() != e[3] || e[26]->GetPrevRight() != e[27])
		return false;
	if (e[27]->GetOrig() != v[16] || e[27]->GetDest() != v[15] || e[27]->GetLeft() != f[11] || e[27]->GetRight() != f[10])
		return false;
	if (e[27]->GetNextOrig() != e[39] || e[27]->GetNextDest() != e[29] || e[27]->GetNextLeft() != e[29] || e[27]->GetNextRight() != e[26])
		return false;
	if (e[27]->GetPrevOrig() != e[26] || e[27]->GetPrevDest() != e[22] || e[27]->GetPrevLeft() != e[28] || e[27]->GetPrevRight() != e[20])
		return false;
	if (e[28]->GetOrig() != v[16] || e[28]->GetDest() != v[17] || e[28]->GetLeft() != f[21] || e[28]->GetRight() != f[11])
		return false;
	if (e[28]->GetNextOrig() != e[23] || e[28]->GetNextDest() != e[40] || e[28]->GetNextLeft() != e[40] || e[28]->GetNextRight() != e[27])
		return false;
	if (e[28]->GetPrevOrig() != e[39] || e[28]->GetPrevDest() != e[33] || e[28]->GetPrevLeft() != e[39] || e[28]->GetPrevRight() != e[29])
		return false;
	if (e[29]->GetOrig() != v[15] || e[29]->GetDest() != v[17] || e[29]->GetLeft() != f[11] || e[29]->GetRight() != f[12])
		return false;
	if (e[29]->GetNextOrig() != e[20] || e[29]->GetNextDest() != e[30] || e[29]->GetNextLeft() != e[28] || e[29]->GetNextRight() != e[22])
		return false;
	if (e[29]->GetPrevOrig() != e[27] || e[29]->GetPrevDest() != e[40] || e[29]->GetPrevLeft() != e[27] || e[29]->GetPrevRight() != e[30])
		return false;
	if (e[30]->GetOrig() != v[1] || e[30]->GetDest() != v[17] || e[30]->GetLeft() != f[12] || e[30]->GetRight() != f[13])
		return false;
	if (e[30]->GetNextOrig() != e[22] || e[30]->GetNextDest() != e[31] || e[30]->GetNextLeft() != e[29] || e[30]->GetNextRight() != e[0])
		return false;
	if (e[30]->GetPrevOrig() != e[1] || e[30]->GetPrevDest() != e[29] || e[30]->GetPrevLeft() != e[22] || e[30]->GetPrevRight() != e[31])
		return false;
	if (e[31]->GetOrig() != v[0] || e[31]->GetDest() != v[17] || e[31]->GetLeft() != f[13] || e[31]->GetRight() != f[14])
		return false;
	if (e[31]->GetNextOrig() != e[9] || e[31]->GetNextDest() != e[32] || e[31]->GetNextLeft() != e[30] || e[31]->GetNextRight() != e[9])
		return false;
	if (e[31]->GetPrevOrig() != e[0] || e[31]->GetPrevDest() != e[30] || e[31]->GetPrevLeft() != e[0] || e[31]->GetPrevRight() != e[32])
		return false;
	if (e[32]->GetOrig() != v[10] || e[32]->GetDest() != v[17] || e[32]->GetLeft() != f[14] || e[32]->GetRight() != f[15])
		return false;
	if (e[32]->GetNextOrig() != e[9] || e[32]->GetNextDest() != e[33] || e[32]->GetNextLeft() != e[31] || e[32]->GetNextRight() != e[10])
		return false;
	if (e[32]->GetPrevOrig() != e[10] || e[32]->GetPrevDest() != e[31] || e[32]->GetPrevLeft() != e[9] || e[32]->GetPrevRight() != e[33])
		return false;
	if (e[33]->GetOrig() != v[17] || e[33]->GetDest() != v[9] || e[33]->GetLeft() != f[22] || e[33]->GetRight() != f[15])
		return false;
	if (e[33]->GetNextOrig() != e[28] || e[33]->GetNextDest() != e[10] || e[33]->GetNextLeft() != e[34] || e[33]->GetNextRight() != e[32])
		return false;
	if (e[33]->GetPrevOrig() != e[32] || e[33]->GetPrevDest() != e[34] || e[33]->GetPrevLeft() != e[40] || e[33]->GetPrevRight() != e[10])
		return false;
	if (e[34]->GetOrig() != v[9] || e[34]->GetDest() != v[18] || e[34]->GetLeft() != f[22] || e[34]->GetRight() != f[16])
		return false;
	if (e[34]->GetNextOrig() != e[33] || e[34]->GetNextDest() != e[35] || e[34]->GetNextLeft() != e[40] || e[34]->GetNextRight() != e[8])
		return false;
	if (e[34]->GetPrevOrig() != e[8] || e[34]->GetPrevDest() != e[40] || e[34]->GetPrevLeft() != e[33] || e[34]->GetPrevRight() != e[35])
		return false;
	if (e[35]->GetOrig() != v[8] || e[35]->GetDest() != v[18] || e[35]->GetLeft() != f[16] || e[35]->GetRight() != f[17])
		return false;
	if (e[35]->GetNextOrig() != e[13] || e[35]->GetNextDest() != e[36] || e[35]->GetNextLeft() != e[34] || e[35]->GetNextRight() != e[13])
		return false;
	if (e[35]->GetPrevOrig() != e[7] || e[35]->GetPrevDest() != e[34] || e[35]->GetPrevLeft() != e[8] || e[35]->GetPrevRight() != e[36])
		return false;
	if (e[36]->GetOrig() != v[12] || e[36]->GetDest() != v[18] || e[36]->GetLeft() != f[17] || e[36]->GetRight() != f[18])
		return false;
	if (e[36]->GetNextOrig() != e[17] || e[36]->GetNextDest() != e[37] || e[36]->GetNextLeft() != e[35] || e[36]->GetNextRight() != e[17])
		return false;
	if (e[36]->GetPrevOrig() != e[13] || e[36]->GetPrevDest() != e[35] || e[36]->GetPrevLeft() != e[13] || e[36]->GetPrevRight() != e[37])
		return false;
	if (e[37]->GetOrig() != v[14] || e[37]->GetDest() != v[18] || e[37]->GetLeft() != f[18] || e[37]->GetRight() != f[19])
		return false;
	if (e[37]->GetNextOrig() != e[45] || e[37]->GetNextDest() != e[38] || e[37]->GetNextLeft() != e[36] || e[37]->GetNextRight() != e[16])
		return false;
	if (e[37]->GetPrevOrig() != e[19] || e[37]->GetPrevDest() != e[36] || e[37]->GetPrevLeft() != e[17] || e[37]->GetPrevRight() != e[38])
		return false;
	if (e[38]->GetOrig() != v[13] || e[38]->GetDest() != v[18] || e[38]->GetLeft() != f[19] || e[38]->GetRight() != f[20])
		return false;
	if (e[38]->GetNextOrig() != e[16] || e[38]->GetNextDest() != e[39] || e[38]->GetNextLeft() != e[37] || e[38]->GetNextRight() != e[23])
		return false;
	if (e[38]->GetPrevOrig() != e[44] || e[38]->GetPrevDest() != e[37] || e[38]->GetPrevLeft() != e[16] || e[38]->GetPrevRight() != e[39])
		return false;
	if (e[39]->GetOrig() != v[16] || e[39]->GetDest() != v[18] || e[39]->GetLeft() != f[20] || e[39]->GetRight() != f[21])
		return false;
	if (e[39]->GetNextOrig() != e[28] || e[39]->GetNextDest() != e[40] || e[39]->GetNextLeft() != e[38] || e[39]->GetNextRight() != e[28])
		return false;
	if (e[39]->GetPrevOrig() != e[27] || e[39]->GetPrevDest() != e[38] || e[39]->GetPrevLeft() != e[23] || e[39]->GetPrevRight() != e[40])
		return false;
	if (e[40]->GetOrig() != v[17] || e[40]->GetDest() != v[18] || e[40]->GetLeft() != f[21] || e[40]->GetRight() != f[22])
		return false;
	if (e[40]->GetNextOrig() != e[29] || e[40]->GetNextDest() != e[34] || e[40]->GetNextLeft() != e[39] || e[40]->GetNextRight() != e[33])
		return false;
	if (e[40]->GetPrevOrig() != e[28] || e[40]->GetPrevDest() != e[39] || e[40]->GetPrevLeft() != e[28] || e[40]->GetPrevRight() != e[34])
		return false;
	if (e[41]->GetOrig() != v[6] || e[41]->GetDest() != v[19] || e[41]->GetLeft() != f[27] || e[41]->GetRight() != f[23])
		return false;
	if (e[41]->GetNextOrig() != e[19] || e[41]->GetNextDest() != e[42] || e[41]->GetNextLeft() != e[45] || e[41]->GetNextRight() != e[5])
		return false;
	if (e[41]->GetPrevOrig() != e[5] || e[41]->GetPrevDest() != e[45] || e[41]->GetPrevLeft() != e[19] || e[41]->GetPrevRight() != e[42])
		return false;
	if (e[42]->GetOrig() != v[5] || e[42]->GetDest() != v[19] || e[42]->GetLeft() != f[23] || e[42]->GetRight() != f[24])
		return false;
	if (e[42]->GetNextOrig() != e[11] || e[42]->GetNextDest() != e[43] || e[42]->GetNextLeft() != e[41] || e[42]->GetNextRight() != e[11])
		return false;
	if (e[42]->GetPrevOrig() != e[4] || e[42]->GetPrevDest() != e[41] || e[42]->GetPrevLeft() != e[5] || e[42]->GetPrevRight() != e[43])
		return false;
	if (e[43]->GetOrig() != v[11] || e[43]->GetDest() != v[19] || e[43]->GetLeft() != f[24] || e[43]->GetRight() != f[25])
		return false;
	if (e[43]->GetNextOrig() != e[15] || e[43]->GetNextDest() != e[44] || e[43]->GetNextLeft() != e[42] || e[43]->GetNextRight() != e[15])
		return false;
	if (e[43]->GetPrevOrig() != e[24] || e[43]->GetPrevDest() != e[42] || e[43]->GetPrevLeft() != e[11] || e[43]->GetPrevRight() != e[44])
		return false;
	if (e[44]->GetOrig() != v[13] || e[44]->GetDest() != v[19] || e[44]->GetLeft() != f[25] || e[44]->GetRight() != f[26])
		return false;
	if (e[44]->GetNextOrig() != e[38] || e[44]->GetNextDest() != e[45] || e[44]->GetNextLeft() != e[43] || e[44]->GetNextRight() != e[16])
		return false;
	if (e[44]->GetPrevOrig() != e[23] || e[44]->GetPrevDest() != e[43] || e[44]->GetPrevLeft() != e[15] || e[44]->GetPrevRight() != e[45])
		return false;
	if (e[45]->GetOrig() != v[14] || e[45]->GetDest() != v[19] || e[45]->GetLeft() != f[26] || e[45]->GetRight() != f[27])
		return false;
	if (e[45]->GetNextOrig() != e[16] || e[45]->GetNextDest() != e[41] || e[45]->GetNextLeft() != e[44] || e[45]->GetNextRight() != e[19])
		return false;
	if (e[45]->GetPrevOrig() != e[37] || e[45]->GetPrevDest() != e[44] || e[45]->GetPrevLeft() != e[16] || e[45]->GetPrevRight() != e[41])
		return false;
	NumTests += 1;
	return true;
}
