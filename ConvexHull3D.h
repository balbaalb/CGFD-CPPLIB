#ifndef ConvexHull3DH
#define ConvexHull3DH
#include <vector>
using namespace std;
#include "Vector3D.h"
#include "Face.h"
#include "Edge.h"
class QuadEdge;
class EdgeContainer;
class EdgeList;
class FaceVisible : public Face
{
	bool visible;
public:
	FaceVisible();
	bool GetVisible() const;
	void SetVisible(bool v);
	Face* Clone() const;
};
enum EDGE_VISIBILITY
{
	EDGE_STATUS_UNKNOWN,
	BOTH_FACES_INVISIBLE,
	RIGHT_FACE_VISIBLE,
	LEFT_FACE_VISIBLE,
	BOTH_FACES_VISIBLE
};
class EdgeVisible : public Edge
{
	EDGE_VISIBILITY visibility;
public:
	EdgeVisible();
	int GetVisiblity() const;
	void SetVisiblity(EDGE_VISIBILITY v);
	Edge* Clone() const;
};
class ConvexHull3D
{
	QuadEdge* hull;
	vector<Vector3D> inputCopy;//before tetrahedron is built input points are stored here
	void CopyBody(const ConvexHull3D& rhs);
	void DeleteBody();
	bool RemoveEdge(EdgeContainer* E);
	bool UpdateFacesVisibility(const Vector3D& p);
	Edge* FindOneBoundaryEdge() const;
	void UpdateEdgesVisibility();
	void DeleteEdges();
	Vertex* AddPyramid(const Vector3D& p, Edge* eBoundary);
	bool PointAcceptable(const Vector3D& p);
	void BuildTetrahedron(const vector<Vector3D>& input, vector<int>& usedPoints);
	QuadEdge* GetHullPointer();
	friend bool tester_ConvexHull3D_1(int& NumTests);
	friend class DelaunayLifting;
public:
	ConvexHull3D();
	ConvexHull3D(const ConvexHull3D& rhs);
	~ConvexHull3D();
	void operator=(const ConvexHull3D& rhs);
	Vertex* AddPoint(const Vector3D& p, int index =-1);
	void Build(const vector<Vector3D>& input);
	QuadEdge* GetHull() const;
	bool GetFaceVisibiliy(int i) const;
	int GetEdgeVisibiliy(int i) const;
	bool TestConvexity() const;
};
bool tester_ConvexHull3D(int& NumTests);
bool tester_ConvexHull3D_1(int& NumTests);
bool tester_ConvexHull3D_2(int& NumTests);
bool tester_ConvexHull3D_PyramidToCube(int& NumTests);
bool tester_ConvexHull3D_3(int& NumTests);
bool tester_ConvexHull3D_Icosahedron(int& NumTests);
bool tester_ConvexHull3D_Dodecahedron(int& NumTests);
#endif