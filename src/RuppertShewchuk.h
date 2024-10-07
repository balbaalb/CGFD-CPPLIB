#ifndef RuppertShewchukH
#define RuppertShewchukH
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include "Triangulation.h"
#include "Vector3D.h"
class Edge;
class Face;
class bPolygon;
class RuppertShewchuk
{
	double LengthScale;
	unique_ptr<Triangulation> triang;
	vector<Vector3D> points;//boundary points are in the [0-lastBP] range, the rest are internal points
	vector<int> segmentPoints;
	vector<TriangulationConstraint> constrainst;
	RuppertShewchuk();
	~RuppertShewchuk();
	void RemovePoint(const Vector3D& point);
	bool IsEncroched(Edge* e, Face* fb);
	void ReBuildInput(QuadEdge* QE, Face* fb, list<Edge*>& newConstraints);
	Vertex* findVertex(Vector3D* newPoint);
	bool IsAcceptable(Vector3D* P);
	int UnEncrochEdges(Vector3D* newPoint = 0);
	Vector3D* SplitOneSkinnyTriangle();
	void UpdateLengthScale();
	void Build(const vector<Vector3D>& BoundaryPoints, int maxIter);
public:
	static double B;
	static double minLengthToMinSegment;
	static Triangulation Tessellate(const vector<Vector3D>& BoundaryPoints, int maxIter);
	static Triangulation Tessellate(const TriangulationInput& input, int maxIter);
};
#endif