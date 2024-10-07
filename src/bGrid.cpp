#include "bGrid.h"
#include "Vector3D.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include "EdgeIterator.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "QuadEdge.h"
bGrid::bGrid(int Nx_, int Ny_, double Lx_, double Ly_) : Nx(Nx_) , Ny(Ny_) , Lx(Lx_), Ly(Ly_)
{
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	this->cellArea = hx * hy;
	this->Mesh2D = new QuadEdge;
	Vertex*** vertices = new Vertex**[Nx + 1];
	for (int i = 0; i <= Nx; ++i){
		vertices[i] = new Vertex*[Ny + 1];
		for (int j = 0; j <= Ny; ++j){
			Vector3D P(double(i) * hx, double(j) * hy);
			vertices[i][j] = this->Mesh2D->AddVertex(P);
		}
	}
	Face*** Faces = new Face**[Nx];
	for (int i = 0; i < Nx; ++i){
		Faces[i] = new Face*[Ny];
		for (int j = 0; j < Ny; ++j){
			Faces[i][j] = this->Mesh2D->AddFace();
		}
	}
	Edge*** horizontalEdges = new Edge**[Nx];
	for (int i = 0; i < Nx; ++i){
		horizontalEdges[i] = new Edge*[Ny + 1];
		for (int j = 0; j <= Ny; ++j){
			horizontalEdges[i][j] = this->Mesh2D->AddEdge();
			this->isHorizontal[horizontalEdges[i][j]] = true;
		}
	}
	Edge*** verticalEdges = new Edge**[Nx + 1];
	for (int i = 0; i <= Nx; ++i){
		verticalEdges[i] = new Edge*[Ny];
		for (int j = 0; j < Ny; ++j){
			verticalEdges[i][j] = this->Mesh2D->AddEdge();
			this->isHorizontal[verticalEdges[i][j]] = false;
		}
	}
	for (int i = 0; i <= Nx; ++i){
		for (int j = 0; j <= Ny; ++j){
			if (i != 0) this->Mesh2D->AddEdgeToVertex(vertices[i][j], horizontalEdges[i - 1][j], true);
			if (i != Nx) this->Mesh2D->AddEdgeToVertex(vertices[i][j], horizontalEdges[i][j], false);
			if (j != 0) this->Mesh2D->AddEdgeToVertex(vertices[i][j], verticalEdges[i][j - 1], true);
			if (j != Ny) this->Mesh2D->AddEdgeToVertex(vertices[i][j], verticalEdges[i][j], false);
		}
	}
	for (int i = 0; i < Nx; ++i){
		for (int j = 0; j < Ny; ++j){
			this->Mesh2D->AddEdgeToFace(Faces[i][j], horizontalEdges[i][j + 1], false);
			this->Mesh2D->AddEdgeToFace(Faces[i][j], verticalEdges[i][j], false);
			this->Mesh2D->AddEdgeToFace(Faces[i][j], horizontalEdges[i][j], true);
			this->Mesh2D->AddEdgeToFace(Faces[i][j], verticalEdges[i + 1][j], true);
		}
	}
	Face* fb = this->Mesh2D->AddFace();
	for (int j = 0; j < Ny; ++j){
		this->Mesh2D->AddEdgeToFace(fb, verticalEdges[0][j], true);
	}
	for (int i = 0; i < Nx; ++i){
		this->Mesh2D->AddEdgeToFace(fb, horizontalEdges[i][Ny], true);
	}
	for (int j = Ny - 1; j >= 0; --j){
		this->Mesh2D->AddEdgeToFace(fb, verticalEdges[Nx][j], false);
	}
	for (int i = Nx - 1; i >= 0; --i){
		this->Mesh2D->AddEdgeToFace(fb, horizontalEdges[i][0], false);
	}
	this->boundary = fb;
	this->Mesh2D->UpdatePrevs();
	for (int i = 0; i <= Nx; ++i){
		delete[] vertices[i];
	}
	delete[] vertices;
	for (int i = 0; i < Nx; ++i){
		delete[] Faces[i];
	}
	delete[] Faces;
	for (int i = 0; i < Nx; ++i){
		delete[] horizontalEdges[i];
	}
	delete[] horizontalEdges;
	for (int i = 0; i <= Nx; ++i){
		delete[] verticalEdges[i];
	}
	delete[] verticalEdges;
	this->Mesh2D->IssueIndices();
}
int bGrid::NumVertices() const
{
	return this->Mesh2D ? this->Mesh2D->NumVertices() : 0;
}
int bGrid::NumCells() const
{
	return this->Mesh2D ? this->Mesh2D->NumFaces() - 1 : 0;
}
int bGrid::NumBoundaryEdges() const
{
	return this->boundary ? this->boundary->GetDegree() : 0;
}
QuadEdge* bGrid::GetMesh2D()
{
	return this->Mesh2D;
}
Face* bGrid::GetBoundary()
{
	return this->boundary;
}
bool bGrid::IsHorizontal(Edge* e)
{
	return this->isHorizontal[e];
}
bool bGrid::IsVertical(Edge* e)
{
	return !this->isHorizontal[e];
}
void bGrid::GetIsHorizontal(unordered_map<Edge*, bool>& isHorizontalMap)
{
	isHorizontalMap.clear();
	isHorizontalMap.rehash(this->Mesh2D->NumEdges());
	EdgeIterator ite(this->Mesh2D);
	Edge* e = ite.Next();
	while (e)
	{
		isHorizontalMap[e] = this->isHorizontal[e];
		e = ite.Next();
	}
}
void bGrid::GetIsVertical(unordered_map<Edge*, bool>& isVerticalMap)
{
	isVerticalMap.clear();
	isVerticalMap.rehash(this->Mesh2D->NumEdges());
	EdgeIterator ite(this->Mesh2D);
	Edge* e = ite.Next();
	while (e)
	{
		isVerticalMap[e] = !this->isHorizontal[e];
		e = ite.Next();
	}
}
double bGrid::GetCellArea()
{
	//assumes the cell is a rectangle
	return this->cellArea;
}
double bGrid::GetLx() const
{
	return Lx;
}
double bGrid::GetLy() const
{
	return Ly;
}
int bGrid::GetNx() const
{
	return Nx;
}
int bGrid::GetNy() const
{
	return Ny;
}
bool bGrid::IsBoundary(Edge* e, GRID_DIRECTION side)
{
	Face* fb = this->GetBoundary();
	if (e->GetRight() == fb || e->GetLeft() == fb) {
		if (side == GRID_DIRECTION::ANY_BOUNDARY)
			return true;
		bool H = this->IsHorizontal(e);
		if (side == GRID_DIRECTION::EAST)
			return !H && e->GetRight() == fb;
		if (side == GRID_DIRECTION::WEST)
			return !H && e->GetLeft() == fb;
		if (side == GRID_DIRECTION::NORTH)
			return H && e->GetLeft() == fb;
		if (side == GRID_DIRECTION::SOUTH)
			return H && e->GetRight() == fb;
	}
	return false;
}
GRID_DIRECTION bGrid::WhichBoundary(Edge* e)
{
	Face* fb = this->GetBoundary();
	if (e->GetRight() == fb || e->GetLeft() == fb) {
		bool H = this->IsHorizontal(e);
		if (!H && e->GetRight() == fb)
			return GRID_DIRECTION::EAST;
		if (!H && e->GetLeft() == fb)
			return GRID_DIRECTION::WEST;
		if (H && e->GetLeft() == fb)
			return GRID_DIRECTION::NORTH;
		if (H && e->GetRight() == fb)
			return GRID_DIRECTION::SOUTH;
	}
	return GRID_DIRECTION::NONE;
}
GRID_DIRECTION bGrid::WhichBoundary(Vertex* v)
{
	EdgeIterator ite(v);
	Edge* e = ite.Next();
	while (e)
	{
		GRID_DIRECTION side = this->WhichBoundary(e);
		if (side != GRID_DIRECTION::NONE)
			return side;
		e = ite.Next();
	}
	return GRID_DIRECTION::NONE;
}
Edge* bGrid::GetOppositeEdge(Edge* e, Face* f)
{
	Edge* eOpp = 0;
	bool found_e = false;
	if (f && e)
	{
		EdgeIterator ite(f);
		Edge* ee = ite.Next();
		while (ee)
		{
			if (!eOpp && ee != e && ee->GetNextLeft() != e && ee->GetNextRight() != e && ee->GetPrevLeft() != e && ee->GetPrevRight() != e)
				eOpp = ee;
			if (ee == e)
				found_e = true;
			ee = ite.Next();
		}
	}
	return found_e ? eOpp : 0;
}