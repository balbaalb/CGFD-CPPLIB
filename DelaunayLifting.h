#ifndef DelaunayLiftingH
#define DelaunayLiftingH
#include <vector>
using namespace std;
#include "Vector3D.h"
#include "Triangulation.h"
class QuadEdge;
class ConvexHull3D;
class Vertex;
class Face;
class DelaunayLifting
{
	ConvexHull3D* lifted;
	Triangulation* DelTriangulation;
	vector<Vertex*> VertexHash;
	void CopyBody(const DelaunayLifting& rhs);
	void DeleteBody();
	static void Lift(const vector<Vector3D>& input, vector<Vector3D>& liftedPoints);
	void UpdateFaceVisibility();//Mark faces that are visible from the top.
	void BuildMesh2D();
	DelaunayLifting();
	DelaunayLifting(const DelaunayLifting& rhs);
	~DelaunayLifting();
	void operator=(const DelaunayLifting& rhs);
	void Build(const vector<Vector3D>& input);
	void Build(const vector<Vector3D>& input, const vector<int>& boundary);
	void Build(const vector<Vector3D>& input, const vector<int>& Boundary, 
		const vector<TriangulationConstraint>& Constraints);
	Triangulation GetTriangulation() const;
public:
	static Triangulation Triangulate(const vector<Vector3D>& inputPoints);
	static Triangulation Triangulate(const vector<Vector3D>& inputPoints, const vector<int>& Boundary);
	static Triangulation Triangulate(const vector<Vector3D>& inputPoints, const vector<int>& Boundary,
		const vector<TriangulationConstraint>& Constraints);
	static Triangulation Triangulate(const TriangulationInput& input);
};
bool tester_DelaunayLifting(int& NumTests);
bool tester_DelaunayLifting_1(int& NumTests);
bool tester_DelaunayLifting_2(int& NumTests);
bool tester_DelaunayLifting_3(int& NumTests);
bool tester_DelaunayLifting_4(int& NumTests);
bool tester_DelaunayLifting_5(int& NumTests);
bool tester_DelaunayLifting_6(int& NumTests);
bool tester_DelaunayLifting_7(int& NumTests);
bool tester_DelaunayLifting_8(int& NumTests);//constrained Delaunay Triangulation
bool tester_DelaunayLifting_9(int& NumTests);//bug revealling case
bool tester_DelaunayLifting_10(int& NumTests);//Shewchuk 2002 examole on page  748
#endif