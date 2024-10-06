#ifndef EdgeIteratorH
#define EdgeIteratorH
class Vertex;
class Face;
class Edge;
class QuadEdge;
class EdgeList;
class EdgeContainer;
class EdgeIterator
{
	EdgeList* el;
	Face* face;
	Vertex* vertex;
	EdgeContainer* container;
	EdgeContainer* containerReported;
	Edge* edgeReported;
	Edge* edge;
	Edge* edge0;
	Edge* NextForEdgeList();
	Edge* NextForFace();
	Edge* NextForVertex();
	bool onLeft;
	bool onLeft0;
	bool onLeftReported;
	bool onOrig;
	bool onOrig0;
	bool onOrigReported;
public:
	EdgeIterator(QuadEdge* QE);
	EdgeIterator(EdgeList* EL);
	EdgeIterator(Face* f);
	EdgeIterator(Vertex* v);
	Edge* Next();
	bool GetOnLeft() const;
	bool GetOnOrig() const;
	EdgeContainer* GetCurrentContainer();
	Edge* GetNext();//Get the next Edge without changing in the iterator
	Edge* GetPrev();//Get the prev Edge without changing in the iterator
};
bool tester_EdgeIterator(int& NumTests);
#endif