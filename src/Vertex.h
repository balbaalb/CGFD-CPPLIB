#ifndef VertexH
#define VertexH
#include "GeoGraphObject.h"
#include "Vector3D.h"
class Edge;
class Vertex : public GeoGraphObject
{
	Edge* edge;
	void CopyBody(const Vertex& v);
	int degree;
	double x, y, z;
public:
	int type;//For debugging only
	int index;//For debugging only
	Vertex* next;
	Vertex* prev;
	Vertex(double X = 0, double Y = 0, double Z = 0);
	Vertex(const Vertex& v);
	Vertex(const Vector3D& v);
	Vertex& operator=(const Vertex& v);
	Edge* GetEdge() const;
	void SetEdge(Edge* e);
	Vector3D GetPoint();//need this one because it's a mandatory part of class GeoGraphObject
	Vector3D GetPoint() const;
	int GetDegree() const;
	int IncreaseDegree();
	int DecreaseDegree();
	int SetDegree(int d);
	void operator<<(const Vector3D& v);
	void SetX(double X);
	void SetY(double Y);
	void SetZ(double Z);
	double GetX() const;
	double GetY() const;
	double GetZ() const;
	virtual Vertex* Clone() const = 0;
};
class VertexBasic : public Vertex
{
public:
	VertexBasic(double X = 0, double Y = 0, double Z = 0);
	VertexBasic(const Vertex& v);
	VertexBasic(const Vector3D& v);
	Vertex* Clone() const;
};
#endif