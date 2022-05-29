#ifndef EdgeListH
#define EdgeListH
#include <unordered_map>
class Edge;
/*template <class T>*/ 
class EdgeContainer
{
public:
	Edge* edge;
	EdgeContainer* next;
	EdgeContainer* prev;
	int index;//For debugging only
	EdgeContainer(Edge* e);
	EdgeContainer(const EdgeContainer& E);
	void operator=(const EdgeContainer& E);
};
/*template <class T>*/ 
class EdgeList
{
	EdgeContainer* container0;
	std::unordered_map<Edge*, EdgeContainer*> edgeContainerHash;
	int NumEdgeContainers;
	void CopyBody(const EdgeList& rhs);
public:
	EdgeList(int Order = 1000);
	EdgeList(const EdgeList& rhs);
	~EdgeList();
	void operator=(const EdgeList& rhs);
	int GetNumEdgeContainers() const;
	void AddContainer(Edge* e);
	void RemoveContainer(EdgeContainer* E);
	void RemoveContainerAndDeleteEdge(EdgeContainer* E);
	void Reset();
	void ResetAndDeleteEdges();
	EdgeContainer* GetContainer(int ei) const;
	EdgeContainer* FindContainer(Edge* e);
	void MoveToEnd(EdgeContainer* E);
};
bool tester_EdgeList(int& NumTests);
#endif