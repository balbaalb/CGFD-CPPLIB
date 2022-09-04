#ifndef EdgeH
#define EdgeH
#include "GeoGraphObject.h"
class Vertex;
class Face;
#include "Vector3D.h"
#include "LineSegment2D.h"
enum EDGE_RELATION
{
	EDGES_DISJOINT,
	EDGES_CONVERGENT,
	EDGES_DIVERGENT,
	EDGES_CONVERGENT_DIVERGENT,
	EDGES_SEQUENTIAL,
	EDGES_CYCLIC,
	EDGES_NONSEQUENTIAL
};
class Edge : public GeoGraphObject
{
	Vertex* Orig;
	Vertex* Dest;
	Face* Right;
	Face* Left;
	Edge* nextOrig;
	Edge* prevOrig;
	Edge* nextDest;
	Edge* prevDest;
	Edge* nextRight;
	Edge* prevRight;
	Edge* nextLeft;
	Edge* prevLeft;
	void Reset();
	void CopyBody(const Edge& e);
	void ReverseFaces();
	void ReverseVertices();
public:
	int type;//For debugging only
	int index;//For debugging only
	Edge();
	Edge(const Edge& e);
	Edge& operator=(const Edge& e);
	void SetOrig(Vertex* v);
	void SetDest(Vertex* v);
	void SetRight(Face* f);
	void SetLeft(Face* f);
	void SetNextOrig(Edge* e);
	void SetPrevOrig(Edge* e);
	void SetNextDest(Edge* e);
	void SetPrevDest(Edge* e);
	void SetNextRight(Edge* e);
	void SetPrevRight(Edge* e);
	void SetNextLeft(Edge* e);
	void SetPrevLeft(Edge* e);
	Vertex* GetOrig() const;
	Vertex* GetDest() const;
	Face* GetRight() const;
	Face* GetLeft() const;
	Edge* GetNextOrig() const;
	Edge* GetPrevOrig() const;
	Edge* GetNextDest() const;
	Edge* GetPrevDest() const;
	Edge* GetNextRight() const;
	Edge* GetPrevRight() const;
	Edge* GetNextLeft() const;
	Edge* GetPrevLeft() const;
	Edge* GetNext(Face* f) const;//needs testing
	Edge* GetNext(Vertex* v) const;//needs testing
	void Reverse();
	Vector3D GetVector() const;
	Vector3D GetNormal_2D() const;//Only works for 2D edges
	double GetLength() const;
	Vector3D GetMidPoint() const;
	Vector3D GetPoint();
	LineSegment2D GetLineSegment2D() const;//needs testing
	EDGE_RELATION RelationTo(Edge* e1) const;
	EDGE_RELATION DirectionTo(Edge* e1) const;
	void UpdateNextRight(Edge* e);
	void UpdateNextLeft(Edge* e);
	void UpdateNextOrig(Edge* e);
	void UpdateNextDest(Edge* e);
	void UpdatePrevRight(Edge* e);
	void UpdatePrevLeft(Edge* e);
	void UpdatePrevOrig(Edge* e);
	void UpdatePrevDest(Edge* e);
	virtual Edge* Clone() const = 0;
	Face* GetOtherFace(Face* f);
	void CloneAdjacency(Edge* e);
};
class EdgeBasic : public Edge
{
public:
	EdgeBasic();
	Edge* Clone() const;
};
bool tester_Edge(int& NumTests);
#endif