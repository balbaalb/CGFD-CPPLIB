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
	DelaunayLifting& operator=(const DelaunayLifting& rhs);
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
#endif