#ifndef bPolygonH
#define bPolygonH
#include <unordered_set>
#include <unordered_map>
#include "Vector3D.h"
#include "Vertex.h"
#include "Triangulation.h"
using namespace std;
class QuadEdge;
class Face;
class EdgeList;
class bPolygon
{
	enum POINT_POLYGON_RELATION
	{
		INSIDE_POLYGON,
		OUTSIDE_POLYGON,
		ON_POLYGON
	};
	/*
	Note1: The ear fining is not optimal because ear status of vertices is not stored and is calculated everytime its needed.
	Note2: The first ear that is found is selected and thus the triangulation geometry is not optimal.
	*/
	QuadEdge* qe;
	bool closed;
	Vertex* latestVertex;
	Vertex* firstVertex;
	Face* inFace;
	Face* outFace;
	unordered_map<Vertex*, Edge*> Arrows;//find the edge such that e->Dest = v;
	void CopyBody(const bPolygon& rhs);
	void DeleteBody();
	void UpdateArrows();
	bool isCCW() const;
	void Reverse();
	void ReArrangeVertices();
	bool isEar(Edge* e);//O(n)
	Edge* FindAnEar();//O(n^2) needs improvement to O(n)
	POINT_POLYGON_RELATION GetRelationTo(const Vector3D& p) const;//O(n)
public:
	bPolygon();
	bPolygon(const bPolygon& rhs);
	~bPolygon();
	void operator=(const bPolygon& rhs);
	QuadEdge* GetQuadEdgePtr() const;
	void AddVertex(const Vector3D& p);
	void AddVerticesAndClose(const vector<Vector3D>& points);
	void Close();
	Face* GetInsideFace();
	bool isEar(Vertex* v);
	Triangulation Triangulate() const;//O(n^3), to improve, change FindAnEar to O(n)
	bool isInside(const Vector3D& p) const;//O(n)
	bool isOn(const Vector3D& p) const;//O(n)
	bool isOnOrInside(const Vector3D& p) const;//O(n)
	bool TestIntegrity();
};
bool tester_bPolygon(int& NumTests);
bool tester_bPolygon_1(int& NumTests);
bool tester_bPolygon_2(int& NumTests);
bool tester_bPolygon_3(int& NumTests);
bool tester_bPolygon_4(int& NumTests);
bool tester_bPolygon_5(int& NumTests);
bool tester_bPolygon_6(int& NumTests);
bool tester_bPolygon_7(int& NumTests);
bool tester_bPolygon_8(int& NumTests);
#endif