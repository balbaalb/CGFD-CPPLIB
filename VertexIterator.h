#ifndef VertexIteratorH
#define VertexIteratorH
class Vertex;
class Face;
class Edge;
class QuadEdge;
class EdgeIterator;
class VertexIterator
{
	QuadEdge* qe;
	Vertex* vertex;
	Vertex* vertex0;
	Vertex* vertexReported;
	Face* face;
	Edge* edge;
	Vertex* NextForQuadEdge();
	Vertex* NextForFace();
	EdgeIterator* ite;
public:
	VertexIterator(QuadEdge* QE);
	VertexIterator(Face* f);
	~VertexIterator();
	Vertex* Next();
	Vertex* GetNext();//Get the next Vertex without changing in the iterator //NEEDS TESTING 
	Vertex* GetPrev();//Get the prev Vertex without changing in the iterator //NEEDS TESTING 
	Vertex* GetReported();
};
bool tester_VertexIterator(int& NumTests);
#endif