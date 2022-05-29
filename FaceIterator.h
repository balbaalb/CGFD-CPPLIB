#ifndef FaceIteratorH
#define FaceIteratorH
class Vertex;
class Face;
class Edge;
class QuadEdge;
class EdgeIterator;
class FaceIterator
{
	QuadEdge* qe;
	Face* face;
	Face* face0;
	Vertex* vertex;
	Edge* edge;
	Face* NextForQuadEdge();
	Face* NextForVertex();
	EdgeIterator* ite;
public:
	FaceIterator(QuadEdge* QE);
	FaceIterator(Vertex* v);
	~FaceIterator();
	Face* Next();
};
bool tester_FaceIterator(int& NumTests);
#endif