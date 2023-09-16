#ifndef QuadEdgeH
#define QuadEdgeH
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;
class Vertex;
class Edge;
class Face;
class EdgeContainer;
class EdgeList;
#include "Vector3D.h"
enum POINT_FACE_VISIBILITY
{
	INVISIBLE,
	COPLANAR,
	VISIBLE
};
class QuadEdge
{
	Vertex* v0;
	int Nv;
	EdgeList* edges;
	Edge* edgePrototype;
	Face* facePrototype;
	Vertex* vertexPrototype;
	Face* face0;
	int Nf;
	int nextVertexIndex;
	bool isPlanar;
	void CopyBody(const QuadEdge& g);//O(n^2) //Has bugs
	void AddEdgeToVertex(Vertex* v, Edge* e1, bool isDest);//Does NOT update Prev relations
	void AddEdgeToFace(Face* f, Edge* e1, bool isLeft);//Does NOT update Prev relations
	void SeparateEdgeAndVertex(Edge* e, Vertex* v);
	void DeleteEdgeAndRightFace_arrangeEdgeFaces(Edge* e, Face* f);
	void RemoveVertex(Vertex* v);
	void RemoveFace(Face* f);
	void RemoveEdge(EdgeContainer* E);
	void RemoveEdge(Edge* e);//Needs testing
	void RemoveLoneVertices();
	friend class bPolygon;
	friend class Triangulation;
	friend class bGrid;
public:
	QuadEdge();
	QuadEdge(const QuadEdge& g);//O(n^2)
	~QuadEdge();
	QuadEdge& operator=(const QuadEdge& g);//O(n^2)
	void SetPrototype(Face* prototype);
	void SetPrototype(Edge* prototype);
	void SetPrototype(Vertex* prototype);
	void Reset();
	Vertex* AddVertex(const Vector3D& v);
	Vertex* AddVertex(const Vertex* v);
	Vertex* GetVertex(int vi);
	Vertex* GetVertexOfIndex(int vi);
	Face* AddFace();
	Edge* AddEdge();
	Edge* GetEdge(int ei);
	Face* GetFace(int fi);
	void UpdatePrevs();
	int NumVertices() const;
	int NumEdges() const;
	int NumFaces() const;
	void GetFaceEdges(Face* f, vector<Edge*>& edgeVec);
	void GetVertexEdges(Vertex* v, vector<Edge*>& edgeVec);
	void GetFaceEdgesReverse(Face* f, vector<Edge*>& edgeVec);
	void GetVertexEdgesReverse(Vertex* v, vector<Edge*>& edgeVec);
	void SetAsDoubleTriangle(const Vector3D& V0, const Vector3D& V1, const Vector3D& V2);
	void SetAsDoubleTriangulatedSquare(const Vector3D& V0, const Vector3D& V1, const Vector3D& V2, const Vector3D& V3);
	Vertex* AddPyramidAndDeleteFace(const Vector3D& v, Face* f0);//return the new vertex;
	bool SetAsTetrahedron(const Vector3D& V0, const Vector3D& V1, const Vector3D& V2, const Vector3D& V3);
	void DeleteEdgeAndRightFace(EdgeContainer* E);//Right face will be absorbed in the left face. Updates Prev relations.
	void DeleteEdgeAndRightFace(Edge* e);
	Vector3D GetNormal(Face* f) const;
	POINT_FACE_VISIBILITY isVisible(Face* f, const Vector3D p) const;
	int UpdateVertexDegree(Vertex* v);//SHOULD BE STATIC !!!!!!
	void UpdateVerticesDegree();
	int GetSumVertexDegrees();
	int GetMaxVertexDegree();
	int GetMinVertexDegree();
	int UpdateFaceDegree(Face* f);//SHOULD BE STATIC !!!!!!
	void UpdateFacesDegree();
	int GetSumFaceDegrees();
	bool isRemovable(Edge* e) const;
	EdgeList* GetEdgeList();
	void SetNextVetexIndex(int index);
	void IssueIndices();
	bool IsPlanar() const;
	void SetPlanar(bool v);
	bool operator==(const QuadEdge& rhs) const;
	bool operator!=(const QuadEdge& rhs) const;
	void ReverseFace(Face* f);//Needs independednt testing
	void MakeAllEdgesLeftTo(Face* f);//Needs independednt testing
	Edge* DivideFaceOnRight(Face* f, Vertex* v1, Vertex* v2);//New face is on right of v1-v2. //Needs independednt testing
	bool IsConnected(Vertex* v1, Vertex* v2) const;
	Edge* GetEdgeConnecting(Vertex* v1, Vertex* v2);
	void print(string address);
	void print(ofstream& file);//clears
	void PrintTestScript(ofstream& file);
	void PrintTestScript(string address);//appends
	bool TestExclusivity();
	bool TestIntegrity();
	Vertex* InsertVertexAtMid(Edge* e);
	double GetMinEdgeLength();
	double GetMaxEdgeLength();
	double GetAvgEdgeLength();
	double GetWeight();
	void Write(string fileName);
	void Read(string fileName);
	void GetPoints(vector<Vector3D>& points);//Needs independednt testing
};
ostream& operator<<(ostream& out, QuadEdge& qe);
istream& operator>>(istream& in, QuadEdge& qe);
bool tester_QuadEdge(int& NumTests);
bool tester_QuadEdge_1(int& NumTests);
bool tester_QuadEdge_Tetrahedron(QuadEdge& Tetrahedron, Vertex* v[4], Face* f[4], Edge* e[4], vector<Vector3D>& points);
bool tester_QuadEdge_2(int& NumTests);
bool tester_QuadEdge_Pyramid(int& NumTests);
bool tester_QuadEdge_SetAsDoubleTriangulatedSquare(int& NumTests);
bool tester_QuadEdge_3(int& NumTests);
bool tester_QuadEdge_4(int& NumTests);//test write to and read from a file
#endif