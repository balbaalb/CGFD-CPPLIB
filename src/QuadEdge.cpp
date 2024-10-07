#include <unordered_set>
#include <iomanip>
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