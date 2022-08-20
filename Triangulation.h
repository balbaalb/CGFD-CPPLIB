#ifndef TriangulationH
#define TriangulationH
#include <string>
#include <memory>
#include "Circle.h"
#include "Triangle.h"
#include "MathUtils.h"
class QuadEdge;
class Face;
class Vertex;
class Edge;
class bPolygon;
class QuadEdgeIndex;
struct TriangulationConstraint
{
	int point1;
	int point2;
};
struct TriangulationInput
{
	vector<Vector3D> points;
	vector<int> boundary;
	vector<TriangulationConstraint> constraints;
};
using namespace std;
class Triangulation
{
protected:
	QuadEdge* Mesh2D;
	Face* boundary;
	double Area;
	void CopyBody(const Triangulation& rhs);
	void DeleteBody();
public:
	Triangulation();
	Triangulation(const Triangulation& rhs);
	~Triangulation();
	Triangulation& operator=(const Triangulation& rhs);
	void operator<<(const QuadEdge& QE);
	QuadEdge* GetMesh2D();
	Face* GetBoundary();
	void GetTriangles(vector<shared_ptr<Triangle> >& triangles) const;
	bool TestIntegrity();
	bool TestDelaunay();
	void PrintTriangulation(string folderName, string extension = "");
	void PrintTriangulation();
	Vertex* GetRightFaceOpposite(Edge* e);
	Vertex* GetLeftFaceOpposite(Edge* e);
	Edge* GetOppositeEdge(Face* f, Vertex* v);
	Face* GetOppositeFace(Face* f, Edge* e);
	Edge* CheckThenFlip(Edge* e);
	bool IsOnBoundary(Vertex* v);//Needs independednt testing
	bool IsFlippable(Edge* e);
	Edge* Flip(Edge* e);//The edge and right face will be replaced (will have new memory address)
	Edge* Connect(Vertex* v1, Vertex* v2);//Connects vertex v1 to v2 with an edge using flips 
	void ShrinkBoundary(Edge* e);
	void ImposeBoundary(const vector<Vertex*>& boundary);
	Triangle GetTriangle(Face* f);
	Circle GetCircumCircle(Face* f);
	double GetCircumCircleRadius(Face* f);
	double GetIncircleRadius(Face* f);
	double GetMinLength(Face* f);
	bool isInside(const Vector3D& p) const;//needs testing //O(n)
	bool isOnOrInside(const Vector3D& p) const;//needs testing //O(n)
	int NumVertices() const;
	int NumTriangles() const;
	int NumBoundaryPoints();
	double GetShapeFactor() const;
	double GetMaxAngle() const;
	double GetMinAngle() const;
	double GetPiAngles() const;
	double GetMaxAngleMeasure() const;
	double GetMinAngleMeasure() const;//Ref: Nguyen 2009
	int GetMaxVertexDegree();
	int GetMinVertexDegree();
	shared_ptr<bPolygon> GetBoundaryPolygon();//Needs testing, not being used yet
	shared_ptr<Triangle> GetTriangle(Face* f) const;
	double GetMinEdgeLength();
	double GetMaxEdgeLength();
	double GetAvgEdgeLength();
	double GetWeight();
	int GetEncrochness();
	double GetSkinnyness();
	double GetMinDistance(Vertex* v);//Min distance between v and other vertices
	double GetCovariance();//Ref: D2R33 Nguyen et al. 2009
	double GetMeshRatio();//Ref: D2R33 Nguyen et al. 2009
	void UpdateArea();
	double GetArea();
	double GetArea(Face* f);
	double GetMinX();
	double GetMaxX();
	double GetMinY();
	double GetMaxY();
	void GetPoints(vector<Vector3D>& points);
	void Write(string fileName);
	void Read(string fileName);
	void Draw(string fileName);
	static Triangulation OffDiagonalGrid(int Nx, int Ny, double Lx, double Ly);
};
bool tester_Triangulation(int& NumTests);
bool tester_Triangulation_1(int& NumTests);
bool tester_Triangulation_1_0(int& NumTests, Triangulation& T, vector<Vector3D>& input);
bool tester_Triangulation_1_1(int& NumTests, Triangulation& T, vector<Vector3D>& input);
bool tester_Triangulation_1_2(int& NumTests, Triangulation& T, vector<Vector3D>& input);
bool tester_Triangulation_1_3(int& NumTests, Triangulation& T, vector<Vector3D>& input);
bool tester_Triangulation_1_4(int& NumTests, Triangulation& T, vector<Vector3D>& input);
bool tester_Triangulation_1_5(int& NumTests, Triangulation& T, vector<Vector3D>& input);
bool tester_Triangulation_1_6(int& NumTests, Triangulation& T, vector<Vector3D>& input);
bool tester_Triangulation_2(int& NumTests);
bool tester_Triangulation_3(int& NumTests);
bool tester_Triangulation_4(int& NumTests);
bool tester_Triangulation_5(int& NumTests);
bool tester_Triangulation_6(int& NumTests);
#endif 