#include "Vertex.h"
#include "Edge.h"
#include "MathUtils.h"
void Vertex::CopyBody(const Vertex& v)
{
	this->edge = 0;
	this->next = 0;
	this->prev = 0;
	this->degree = 0;
	this->x = v.x;
	this->y = v.y;
	this->z = v.z;
	this->index = -1;
}
Vertex::Vertex(double X, double Y, double Z)
{
	this->edge = 0;
	this->next = 0;
	this->prev = 0;
	this->degree = 0;
	this->type = 0;
	this->x = X;
	this->y = Y;
	this->z = Z;
	this->index = -1;
}
Vertex::Vertex(const Vertex& v)
{
	this->CopyBody(v);
	this->type = 0;
}
Vertex::Vertex(const Vector3D& v)
{
	this->edge = 0;
	this->next = 0;
	this->prev = 0;
	this->type = 0;
	this->x = v(0);
	this->y = v(1);
	this->z = v(2);
	this->index = -1;
}
Vertex& Vertex::operator=(const Vertex& v)
{
	this->CopyBody(v);
	this->type = 0;
	return *this;
}
Edge* Vertex::GetEdge() const
{
	return this->edge;
}
void Vertex::SetEdge(Edge* e)
{
	this->edge = e;
}
Vector3D Vertex::GetPoint()//need this one because it's a mandatory part of class GeoGraphObject
{
	Vector3D u(this->x, this->y, this->z);
	return u;
}
Vector3D Vertex::GetPoint() const
{
	Vector3D u(this->x, this->y, this->z);
	return u;
}
int Vertex::GetDegree() const
{
	return this->degree;
}
int Vertex::IncreaseDegree()
{
	return (++this->degree);
}
int Vertex::DecreaseDegree()
{
	return (--this->degree);
}
int Vertex::SetDegree(int d)
{
	return (this->degree = d);
}
void Vertex::operator<<(const Vector3D& v)
{
	this->x = v(0);
	this->y = v(1);
	this->z = v(2);
}
void Vertex::SetX(double X)
{
	this->x = X;
}
void Vertex::SetY(double Y)
{
	this->y = Y;
}
void Vertex::SetZ(double Z)
{
	this->z = Z;
}
double Vertex::GetX() const
{
	return this->x;
}
double Vertex::GetY() const
{
	return this->y;
}
double Vertex::GetZ() const
{
	return this->z;
}
VertexBasic::VertexBasic(double x, double y, double z) : Vertex(x,y,z)
{
	this->type = 1;
}
VertexBasic::VertexBasic(const Vertex& v) : Vertex(v)
{
	this->type = 1;
}
VertexBasic::VertexBasic(const Vector3D& v) : Vertex(v)
{
	this->type = 1;
}
Vertex* VertexBasic::Clone() const
{
	return (new VertexBasic);
}