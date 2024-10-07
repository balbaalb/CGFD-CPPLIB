#include "VertexIterator.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "QuadEdge.h"
#include "EdgeIterator.h"
Vertex* VertexIterator::NextForQuadEdge()
{
	Vertex* v = this->vertex;
	if (this->vertex)
		this->vertex = this->vertex->next;
	if (this->vertex == this->vertex0)
		this->vertex = 0;
	return v;
}
Vertex* VertexIterator::NextForFace()
{
	Vertex* v = 0;
	this->edge = this->ite->Next();
	if (this->edge)
	{
		if (this->ite->GetOnLeft())
			v = this->edge->GetOrig();
		else
			v = this->edge->GetDest();
	}
	return v;
}
VertexIterator::VertexIterator(QuadEdge* QE)
{
	this->qe = QE;
	this->vertex0 = this->vertex = this->vertexReported =  this->qe->GetVertex(0);
	this->face = 0;
	this->edge = 0;
	this->ite = 0;
}
VertexIterator::VertexIterator(Face* f)
{
	this->qe = 0;
	this->vertex0 = this->vertex = this->vertexReported = 0;
	this->face = f;
	this->ite = new EdgeIterator(f);
	this->edge = 0;
}
VertexIterator::~VertexIterator()
{
	delete this->ite;
	this->ite = 0;
}
Vertex* VertexIterator::Next()
{
	if (this->qe)
		this->vertexReported =  this->NextForQuadEdge();
	else
		this->vertexReported = this->NextForFace();
	return this->vertexReported;
}
Vertex* VertexIterator::GetNext()
{
	if (this->qe)
		return this->vertexReported->next;
	else
	{
		Edge* e = this->ite->GetNext();
		bool onLeftNext = this->ite->GetOnLeft();
		if (this->edge->DirectionTo(e) == EDGES_NONSEQUENTIAL)
			onLeftNext = !onLeftNext;
		if (onLeftNext)
			return e->GetOrig();
		else
			return e->GetDest();
	}
}
Vertex* VertexIterator::GetPrev()
{
	if (this->qe)
		return this->vertexReported->prev;
	else
	{
		Edge* e = this->ite->GetPrev();
		bool onLeftPrev = this->ite->GetOnLeft();
		if (this->edge->DirectionTo(e) == EDGES_NONSEQUENTIAL)
			onLeftPrev = !onLeftPrev;
		if (onLeftPrev)
			return e->GetOrig();
		else
			return e->GetDest();
	}
}
Vertex* VertexIterator::GetReported()
{
	return this->vertexReported;
}