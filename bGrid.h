#ifndef bGridH
#define bGridH
#include <unordered_map>
using namespace std;
class QuadEdge;
class Face;
class Vertex;
class Edge;
enum GRID_DIRECTION
{
	EAST,
	WEST,
	NORTH,
	SOUTH,
	TOP,
	BOTTOM,
	ANY_BOUNDARY,
	NONE
};
class bGrid
{
	QuadEdge* Mesh2D;
	Face* boundary;
	unordered_map<Edge*, bool> isHorizontal;
	double cellArea;
	int Nx, Ny;
	double Lx, Ly;
public:
	bGrid(int Nx_, int Ny_, double Lx_, double Ly_);
	int NumVertices() const;
	int NumBoundaryEdges() const;
	int NumCells() const;
	QuadEdge* GetMesh2D();
	Face* GetBoundary();
	bool IsHorizontal(Edge* e);
	bool IsVertical(Edge* e);
	void GetIsHorizontal(unordered_map<Edge*, bool>& isHorizontalMap);
	void GetIsVertical(unordered_map<Edge*, bool>& isVerticalMap);
	double GetCellArea();
	double GetLx() const;
	double GetLy() const;
	int GetNx() const;
	int GetNy() const;
	bool IsBoundary(Edge* e, GRID_DIRECTION side = GRID_DIRECTION::ANY_BOUNDARY);
	GRID_DIRECTION WhichBoundary(Edge* e);
	GRID_DIRECTION WhichBoundary(Vertex* v);
};
bool tester_bGrid(int& NumTests);
#endif