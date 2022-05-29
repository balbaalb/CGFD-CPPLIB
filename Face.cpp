#include "Face.h"
#include "Edge.h"
#include "Vertex.h"
#include "VertexIterator.h"
void Face::CopyBody(const Face& f)
{
	this->edge = 0;
	this->next = 0;
	this->prev = 0;
	this->type = 0;
	this->degree = 0;
	this->index = -1;
}
Face::Face()
{
	this->edge = 0;
	this->next = 0;
	this->prev = 0;
	this->type = 0;
	this->degree = 0;
	this->index = -1;
}
Face::Face(const Face& f)
{
	this->CopyBody(f);
}
void Face::operator=(const Face& f)
{
	this->CopyBody(f);
}
Edge* Face::GetEdge() const
{
	return this->edge;
}
void Face::SetEdge(Edge* e)
{
	this->edge = e;
}
int Face::GetDegree() const
{
	return this->degree;
}
int Face::IncrementDegree()
{
	return (++this->degree);
}
int Face::DecrementDegree()
{
	return (--this->degree);
}
int Face::SetDegree(int d)
{
	return (this->degree = d);
}
Vector3D Face::GetCenteriod()
{
	Vector3D C;
	int N = 0;
	VertexIterator itv(this);
	Vertex* v = itv.Next();
	while (v)
	{
		++N;
		C = C + v->GetPoint();
		v = itv.Next();
	}
	if (N > 0)
		C = C / double(N);
	return C;
}
Vector3D Face::GetPoint()
{
	return this->GetCenteriod();
}
FaceBasic::FaceBasic() : Face()
{
	this->type = 1;
}
Face* FaceBasic::Clone() const
{
	return (new FaceBasic);
}
bool tester_Face(int& NumTests)
{
	EdgeBasic e1;
	FaceBasic f1;
	FaceBasic f2(f1), f3;
	f3 = f2;
	f3.SetEdge(&e1);
	if (!f3.GetEdge() || f3.GetEdge() != &e1)
		return false;
	NumTests += 1;
	return true;
}