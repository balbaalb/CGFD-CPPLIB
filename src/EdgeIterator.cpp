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