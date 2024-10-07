#include "FaceIterator.h"
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "QuadEdge.h"
#include "EdgeIterator.h"
Face* FaceIterator::NextForQuadEdge()
{
	Face* f = this->face;
	if (this->face)
		this->face = this->face->next;
	if (this->face == this->face0)
		this->face = 0;
	return f;
}
Face* FaceIterator::NextForVertex()
{
	Face* f = 0;
	if (this->edge)
	{
		if (this->ite->GetOnOrig())
			f = this->edge->GetLeft();
		else
			f = this->edge->GetRight();
		this->edge = this->ite->Next();
	}
	return f;
}
FaceIterator::FaceIterator(QuadEdge* QE)
{
	this->qe = QE;
	this->face0 = this->face = this->qe->GetFace(0);
	this->vertex = 0;
	this->edge = 0;
	this->ite = 0;
}
FaceIterator::FaceIterator(Vertex* v)
{
	this->qe = 0;
	this->face0 = this->face = 0;
	this->vertex = v;
	this->ite = new EdgeIterator(v);
	this->edge = this->ite->Next();
}
FaceIterator::~FaceIterator()
{
	delete this->ite;
	this->ite = 0;
}
Face* FaceIterator::Next()
{
	if (this->qe)
		return this->NextForQuadEdge();
	else
		return this->NextForVertex();
}