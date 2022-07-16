#ifndef NodeCompositeH
#define NodeCompositeH
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>
using namespace std;
#include "Node.h"
#include "MathUtils.h"
#include "Vector.h"
class GeoGraphObject;
class Face;
class Edge;
class Vertex;
class QuadEdge;
class SquareMatrix;
class Vector;
enum NODE_COMPOSITE_TYPE
{
	CELLS, EDGES, VERTICES, CELLS_AND_BOUNDARY, EDGES_AND_BOUNDARY
};
class NodeComposite
{
	QuadEdge* qe;
	Face* fb;
	NODE_COMPOSITE_TYPE type;
	unordered_map<GeoGraphObject*, Node*> NodeHash;
	unordered_map<GeoGraphObject*, bool> isDependent;
	unordered_map<GeoGraphObject*, double> constValues;
	unordered_map<GeoGraphObject*, double> constNormalGrad;
	vector<Node*> nodes;
	vector<Node*> DependentNodes;
	int index, index2;
	SquareMatrix* K;
	Vector* C;
	Vector* X; //stores the last solution
	typedef unordered_map<Edge*, bool> EdgeConditionMap;
	EdgeConditionMap* EdgeCondition;
	bool initilized{ false };
	void CopyBody(const NodeComposite& rhs);
	void Reset();
	void AddNode(GeoGraphObject* gPtr);
	void AddDependentNode(GeoGraphObject* gPtr);
	void populate_basedOn_CELLS();
	void populate_basedOn_EDGES();
	void populate_basedOn_VERTICES();
	void populate_basedOn_CELLS_AND_BOUNDARY();
	void populate_basedOn_EDGES_AND_BOUNDARY();
	void ApplyBC_internal(string caller);
public:
	NodeComposite();
	NodeComposite(const NodeComposite& rhs);
	~NodeComposite();
	void operator=(const NodeComposite& rhs);
	void SetEdgeCondition(QuadEdge* QE, EdgeConditionMap* f);//use this before intilization
	void Initialize(QuadEdge* QE, Face* Fb, NODE_COMPOSITE_TYPE Type);
	void InitializeEquations();
	void AddToK(GeoGraphObject* row, GeoGraphObject* col, double value);
	void AddToC(GeoGraphObject* row, double value);
	void SetValueInEquations(GeoGraphObject* row, double value);
	void SetValue(GeoGraphObject* objPtr, double value);
	void SetValue(GeoGraphObject* objPtr, function<double(const Vector3D& P)> f);
	void AddToValue(GeoGraphObject* objPtr, double value);
	void SetConstantValue(GeoGraphObject* objPtr, double value);
	void SetConstantValue(GeoGraphObject* objPtr, function<double(const Vector3D& P)> f);
	void SetConstantNormalGradient(GeoGraphObject* objPtr, double value);
	void SetConstantNormalGradient(GeoGraphObject* objPtr, function<double(const Vector3D& P)> f);
	void SetValueAllNodes(function<double(const Vector3D& P)> func, bool basedOnVertices = true);
	void Solve();
	double SolveAndUpdate(double underRelaxation = 1);
	double SolveAndAdd(double underRelaxation = 1);
	void Populate();
	Node* GetNode(GeoGraphObject* objPtr);
	double GetValue(GeoGraphObject* objPtr);
	double GetValue(GeoGraphObject* objPtr, bool& isValid);
	double GetValue(const Vector3D& P);//O(n)
	int GetNodeNumber(GeoGraphObject* objPtr);
	void GetResults_Vertices(vector<Node*>& results);
	void GetResults_Edges(vector<Node*>& results);
	void GetResults_Cells(vector<Node*>& results);
	bool IsIndependent(GeoGraphObject* objPtr);
	bool IsIndependent(Node* n);
	double GetSolution(GeoGraphObject* objPtr);
	void GetK(SquareMatrix* Kr);
	void GetC(Vector* Cr);
	double GetK(GeoGraphObject* row, GeoGraphObject* col);
	double GetC(GeoGraphObject* row);
	void GetX(Vector* Xr);
	void PrintEquations(ofstream& f, string label);
	bool IsAtBoundary(Vertex* v);
	void ApplyBC();
	bool IsInitialized() const;
	bool IsConstValue(GeoGraphObject* objPtr);
};
bool tester_NodeComposite(int& NumTests);
#endif
