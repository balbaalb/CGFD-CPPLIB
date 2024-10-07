#ifndef QuadEdgeIndexH
#define QuadEdgeIndexH
#include <unordered_map>
#include <vector>
using namespace std;
class QuadEdge;
class Vertex;
class Edge;
class Face;
class QuadEdgeIndex
{
	QuadEdge* qe;
	vector<Vertex*> VertexArr;
	vector<Face*> FaceArr;
	vector<Edge*> EdgeArr;
	unordered_map<Vertex*, int> VertexMap;
	unordered_map<Face*, int> FaceMap;
	unordered_map<Edge*, int> EdgeMap;
public:
	QuadEdgeIndex(QuadEdge* qe_);
	Vertex* GetVertex(int vi);
	Face* GetFace(int fi);
	Edge* GetEdge(int ei);
	int GetIndex(Vertex* v);
	int GetIndex(Face* f);
	int GetIndex(Edge* e);
};
#endif
