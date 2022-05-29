#include "EdgeIterator.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "QuadEdge.h"
#include "EdgeList.h"
Edge* EdgeIterator::NextForEdgeList()
{
	Edge* e = 0;
	if (this->container)
	{
		e = this->container->edge;
		this->containerReported = this->container;
		this->container = this->container->next;
	}
	if (this->container == this->el->GetContainer(0))
		this->container = 0;
	return e;
}
Edge* EdgeIterator::NextForFace()
{
	Edge* e = this->edge;
	this->onLeftReported = this->onLeft;
	if (this->edge)
	{
		Edge* eNext = 0;
		if (this->onLeft)
			eNext = this->edge->GetNextLeft();
		else
			eNext = this->edge->GetNextRight();
		if (e->DirectionTo(eNext) == EDGES_NONSEQUENTIAL)
			this->onLeft = !this->onLeft;
		this->edge = eNext;
	}
	if (this->edge == this->edge0 && this->onLeft == this->onLeft0)
		this->edge = 0;
	return e;
}
Edge* EdgeIterator::NextForVertex()
{
	Edge* e = this->edge;
	this->onOrigReported = this->onOrig;
	if (this->edge)
	{
		Edge* eNext = 0;
		if (this->onOrig)
			eNext = this->edge->GetNextOrig();
		else
			eNext = this->edge->GetNextDest();
		if (e->DirectionTo(eNext) == EDGES_SEQUENTIAL)
			this->onOrig = !this->onOrig;
		this->edge = eNext;
	}
	if (this->edge == this->edge0 && this->onOrig == this->onOrig0)
		this->edge = 0;
	return e;
}
EdgeIterator::EdgeIterator(QuadEdge* QE) : EdgeIterator(QE->GetEdgeList())
{

}
EdgeIterator::EdgeIterator(EdgeList* EL)
{
	this->el = EL;
	this->container = this->containerReported = EL ? EL->GetContainer(0) : 0;
	this->face = 0;
	this->vertex = 0;
	this->edge = this->edge0 = this->edgeReported = 0;
}
EdgeIterator::EdgeIterator(Face* f)
{
	this->el = 0;
	this->container = this->containerReported = 0;
	this->face = f;
	this->vertex = 0;
	this->edge0 = this->edgeReported = this->edge = this->face->GetEdge();
	if (this->edge0)
		this->onLeft0 = this->onLeft = this->onLeftReported = (this->edge0->GetLeft() == f);
}
EdgeIterator::EdgeIterator(Vertex* v)
{
	this->el = 0;
	this->container = this->containerReported = 0;
	this->face = 0;
	this->vertex = v;
	this->edge = this->edge0 = this->edgeReported = this->vertex->GetEdge();
	if (this->edge0)
		this->onOrig0 = this->onOrig = this->onOrigReported = (this->edge0->GetOrig() == v);
}
Edge* EdgeIterator::Next()
{
	if (this->el)
		this->edgeReported = this->NextForEdgeList();
	else if (this->face)
		this->edgeReported =  this->NextForFace();
	else if (this->vertex)
		this->edgeReported =  this->NextForVertex();
	return this->edgeReported;
}
bool EdgeIterator::GetOnLeft() const
{
	return this->onLeftReported;
}
bool EdgeIterator::GetOnOrig() const
{
	return this->onOrigReported;
}
EdgeContainer* EdgeIterator::GetCurrentContainer()
{
	return this->containerReported;
}
Edge* EdgeIterator::GetNext()//only works for EdgeIterator(FACE*)
{
	if (this->el)
		return (this->containerReported ? this->containerReported->next->edge : 0);
	if (this->face){
		if (this->GetOnLeft())
			return this->edgeReported->GetNextLeft();
		else
			return this->edgeReported->GetNextRight();
	}
	//Vertex:
	if (this->GetOnOrig())
		return this->edgeReported->GetNextOrig();
	else
		return this->edgeReported->GetNextDest();

}
Edge* EdgeIterator::GetPrev()//only works for EdgeIterator(FACE*)
{
	if (this->el)
		return (this->containerReported ? this->containerReported->prev->edge : 0);
	if (this->face){
		if (this->GetOnLeft())
			return this->edgeReported->GetPrevLeft();
		else
			return this->edgeReported->GetPrevRight();
	}
	//Vertex:
	if (this->GetOnOrig())
		return this->edgeReported->GetPrevOrig();
	else
		return this->edgeReported->GetPrevDest();

}
bool tester_EdgeIterator(int& NumTests)
{
	Vector3D P0, P1(1), P2(0, 1), P3(0, 0, 1);
	QuadEdge T;
	T.SetAsTetrahedron(P0, P1, P2, P3);
	EdgeIterator ite(&T);
	vector<Edge*> e;
	Edge* ep = ite.Next();
	while (ep)
	{
		e.push_back(ep);
		ep = ite.Next();
	}
	if (e.size() != 6)
		return false;
	if (ite.Next() || ite.Next())
		return false; 
	for (int i = 0; i < 5; ++i)
	{
		for (int j = i + 1; j < 6; ++j)
		{
			if (e[i] == e[j])
				return false;
		}
	}
	//Testing GetNext() and GetPrev()
	EdgeIterator ite2(&T);
	if (ite2.Next() != e[0] || ite2.GetNext() != e[1] || ite2.GetPrev() != e[5])
		return false;
	if (ite2.Next() != e[1] || ite2.GetNext() != e[2] || ite2.GetPrev() != e[0])
		return false;
	if (ite2.Next() != e[2] || ite2.GetNext() != e[3] || ite2.GetPrev() != e[1])
		return false;
	if (ite2.Next() != e[3] || ite2.GetNext() != e[4] || ite2.GetPrev() != e[2])
		return false;
	if (ite2.Next() != e[4] || ite2.GetNext() != e[5] || ite2.GetPrev() != e[3])
		return false;
	if (ite2.Next() != e[5] || ite2.GetNext() != e[0] || ite2.GetPrev() != e[4])
		return false;
	if (ite2.Next())
		return false;
	//The next tests rely on the order that tetrahedron vertices, faces and edges are ordered
	Vertex* v[4];
	Face* f[4];
	Edge* ee[6];
	for (int i = 0; i < 4; ++i)
	{
		f[i] = T.GetFace(i);
		v[i] = T.GetVertex(i);
	}
	for (int i = 0; i < 6; ++i)
	{
		ee[i] = T.GetEdge(i);
		if (e[i] != ee[i])
			return false;
	}
	EdgeIterator itev0(v[0]);
	if (itev0.Next() != ee[0] || itev0.Next() != ee[3] || itev0.Next() != ee[2])
		return false;
	if (itev0.Next() || itev0.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itev00(v[0]);
	if (itev00.Next() != ee[0] || itev00.GetNext() != ee[3] || itev00.GetPrev() != ee[2])
		return false;
	if (itev00.Next() != ee[3] || itev00.GetNext() != ee[2] || itev00.GetPrev() != ee[0])
		return false;
	if (itev00.Next() != ee[2] || itev00.GetNext() != ee[0] || itev00.GetPrev() != ee[3])
		return false;
	EdgeIterator itev1(v[1]);
	if (itev1.Next() != ee[1] || itev1.Next() != ee[4] || itev1.Next() != ee[0])
		return false;
	if (itev1.Next() || itev1.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itev11(v[1]);
	if (itev11.Next() != ee[1] || itev11.GetNext() != ee[4] || itev11.GetPrev() != ee[0])
		return false;
	if (itev11.Next() != ee[4] || itev11.GetNext() != ee[0] || itev11.GetPrev() != ee[1])
		return false;
	if (itev11.Next() != ee[0] || itev11.GetNext() != ee[1] || itev11.GetPrev() != ee[4])
		return false;
	EdgeIterator itev2(v[2]);
	if (itev2.Next() != ee[2] || itev2.Next() != ee[5] || itev2.Next() != ee[1])
		return false;
	if (itev2.Next() || itev2.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itev22(v[2]);
	if (itev22.Next() != ee[2] || itev22.GetNext() != ee[5] || itev22.GetPrev() != ee[1])
		return false;
	if (itev22.Next() != ee[5] || itev22.GetNext() != ee[1] || itev22.GetPrev() != ee[2])
		return false;
	if (itev22.Next() != ee[1] || itev22.GetNext() != ee[2] || itev22.GetPrev() != ee[5])
		return false;
	EdgeIterator itev3(v[3]);
	if (itev3.Next() != ee[5] || itev3.Next() != ee[3] || itev3.Next() != ee[4])
		return false;
	if (itev3.Next() || itev3.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itev33(v[3]);
	if (itev33.Next() != ee[5] || itev33.GetNext() != ee[3] || itev33.GetPrev() != ee[4])
		return false;
	if (itev33.Next() != ee[3] || itev33.GetNext() != ee[4] || itev33.GetPrev() != ee[5])
		return false;
	if (itev33.Next() != ee[4] || itev33.GetNext() != ee[5] || itev33.GetPrev() != ee[3])
		return false;
	//------------------------------------------------------------------------------
	EdgeIterator itef0(f[0]);
	if (itef0.Next() != ee[0] || itef0.Next() != ee[2] || itef0.Next() != ee[1])
		return false;
	if (itef0.Next() || itef0.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itef00(f[0]);
	if (itef00.Next() != ee[0] || itef00.GetNext() != ee[2] || itef00.GetPrev() != ee[1])
		return false;
	if (itef00.Next() != ee[2] || itef00.GetNext() != ee[1] || itef00.GetPrev() != ee[0])
		return false;
	if (itef00.Next() != ee[1] || itef00.GetNext() != ee[0] || itef00.GetPrev() != ee[2])
		return false;
	EdgeIterator itef1(f[1]);
	if (itef1.Next() != ee[0] || itef1.Next() != ee[4] || itef1.Next() != ee[3])
		return false;
	if (itef1.Next() || itef1.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itef11(f[1]);
	if (itef11.Next() != ee[0] || itef11.GetNext() != ee[4] || itef11.GetPrev() != ee[3])
		return false;
	if (itef11.Next() != ee[4] || itef11.GetNext() != ee[3] || itef11.GetPrev() != ee[0])
		return false;
	if (itef11.Next() != ee[3] || itef11.GetNext() != ee[0] || itef11.GetPrev() != ee[4])
		return false;
	EdgeIterator itef2(f[2]);
	if (itef2.Next() != ee[1] || itef2.Next() != ee[5] || itef2.Next() != ee[4])
		return false;
	if (itef2.Next() || itef2.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itef22(f[2]);
	if (itef22.Next() != ee[1] || itef22.GetNext() != ee[5] || itef22.GetPrev() != ee[4])
		return false;
	if (itef22.Next() != ee[5] || itef22.GetNext() != ee[4] || itef22.GetPrev() != ee[1])
		return false;
	if (itef22.Next() != ee[4] || itef22.GetNext() != ee[1] || itef22.GetPrev() != ee[5])
		return false;
	EdgeIterator itef3(f[3]);
	if (itef3.Next() != ee[2] || itef3.Next() != ee[3] || itef3.Next() != ee[5])
		return false;
	if (itef3.Next() || itef3.Next())
		return false;
	//Testing GetNext and GetPrev
	EdgeIterator itef33(f[3]);
	if (itef33.Next() != ee[2] || itef33.GetNext() != ee[3] || itef33.GetPrev() != ee[5])
		return false;
	if (itef33.Next() != ee[3] || itef33.GetNext() != ee[5] || itef33.GetPrev() != ee[2])
		return false;
	if (itef33.Next() != ee[5] || itef33.GetNext() != ee[2] || itef33.GetPrev() != ee[3])
		return false;
	//-------------------------------------------------------------------------
	EdgeList eL;
	EdgeBasic e0, e1, *e2;
	e2 = new EdgeBasic;
	eL.AddContainer(&e0);
	eL.AddContainer(&e1);
	eL.AddContainer(e2);
	EdgeIterator itE(&eL);
	if (itE.Next() != &e0 || itE.Next() != &e1 || itE.Next() != e2)
		return false;
	if (itE.Next() || itE.Next())
		return false;
	eL.Reset();
	EdgeIterator itE2(&eL);
	if (itE2.Next() || itE2.Next())
		return false;
	//------------------------------------------------------
	eL.AddContainer(&e0);
	eL.AddContainer(&e1);
	eL.AddContainer(e2);
	EdgeBasic e3;
	eL.AddContainer(&e3);
	EdgeIterator* itE3 = new EdgeIterator(&eL);
	if (itE3->Next() != &e0 || itE3->Next() != &e1)
		return false;
	EdgeContainer* E1 = itE3->GetCurrentContainer();
	if (E1->edge != &e1)
		return false;
	eL.MoveToEnd(E1);
	if (itE3->Next() != e2 || itE3->Next() != &e3)
		return false;
	if (itE3->Next() != &e1 || itE3->Next())
		return false;
	if (itE3->Next() || itE3->Next())
		return false;
	delete itE3;
	itE3 = new EdgeIterator(&eL);
	if (itE3->Next() != &e0)
		return false;
	EdgeContainer* E0 = itE3->GetCurrentContainer();
	if (E0->edge != &e0)
		return false;
	eL.MoveToEnd(E0);
	if (itE3->Next() != e2 || itE3->Next() != &e3 || itE3->Next() != &e1)
		return false;
	if (itE3->Next() != &e0 || itE3->Next())
		return false;
	if (itE3->Next() || itE3->Next())
		return false;
	delete e2;
	NumTests += 1;
	return true;
}