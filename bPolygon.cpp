#include "bPolygon.h"
#include "QuadEdge.h"
#include "EdgeIterator.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "Vertex.h"
#include "Edge.h"
#include "LineSegment2D.h"
#include "MathUtils.h"
void bPolygon::CopyBody(const bPolygon& rhs)
{
	this->qe = 0;
	this->inFace = 0;
	this->outFace = 0;
	if (rhs.qe){
		this->qe = new QuadEdge(*rhs.qe);
		if (rhs.inFace == rhs.qe->GetFace(0))
		{
			this->inFace = this->qe->GetFace(0);
			this->outFace = this->qe->GetFace(1);
		}
		else
		{
			this->inFace = this->qe->GetFace(1);
			this->outFace = this->qe->GetFace(0);
		}
	}
	this->closed = rhs.closed;
	this->latestVertex = 0;
	this->firstVertex = 0;
	if (this->closed)
	{
		this->UpdateArrows();
	}
}
void bPolygon::DeleteBody()
{
	delete this->qe;
	this->qe = 0;
	this->closed = false;
	this->latestVertex = 0;
	this->firstVertex = 0;
	this->inFace = 0;
	this->outFace = 0;
	this->Arrows.empty();
}
void bPolygon::UpdateArrows()
{
	if (this->closed){
		int Df = this->qe->UpdateFaceDegree(this->inFace);
		this->Arrows.empty();
		this->Arrows.rehash(Df);
		EdgeIterator ite(this->inFace);
		Edge* e = ite.Next();
		while (e){
			Vertex* v = e->GetDest();
			this->Arrows[v] = e;
			e = ite.Next();
		}
	}
}
bool bPolygon::isCCW() const
{
	double sumTheta = 0;
	Vertex* a = 0;
	Vertex* b = 0;
	Vertex* c = 0;
	Vertex* v0 = 0;
	Vertex* v1 = 0;
	VertexIterator itv(this->qe);
	Vertex* v = itv.Next();
	while (v){
		a = b;
		b = c;
		c = v;
		if (a && b && c)
			sumTheta += Vector3D::RotationAngle_z(a->GetPoint(), b->GetPoint(), c->GetPoint());
		if (!v0)
			v0 = v;
		else if (!v1)
			v1 = v;
		v = itv.Next();
	}
	a = b;
	b = c;
	c = v0;
	sumTheta += Vector3D::RotationAngle_z(a->GetPoint(), b->GetPoint(), c->GetPoint());
	a = b;
	b = c;
	c = v1;
	sumTheta += Vector3D::RotationAngle_z(a->GetPoint(), b->GetPoint(), c->GetPoint());
	return (sumTheta > 0);
}
void bPolygon::Reverse()
{
	if (this->closed)
	{
		EdgeIterator ite(this->qe);
		Edge* e = ite.Next();
		while (e){
			e->Reverse();
			e = ite.Next();
		}
		Vertex* v = this->firstVertex;
		this->firstVertex = this->latestVertex;
		this->latestVertex = v;
		this->qe->ReverseFace(this->inFace);
		this->qe->ReverseFace(this->outFace);
	}
}
void bPolygon::ReArrangeVertices()
{
	int Nv = this->qe->NumVertices();
	vector<Vertex*> vArray;
	vArray.resize(Nv, 0);
	VertexIterator itv(this->inFace);
	Vertex* v = itv.Next();
	int i = 0;
	while (v){
		if (i >= Nv)
			throw "bPolygon::ReArrangeVertices()";
		vArray[i] = v;
		++i;
		v = itv.Next();
	}
	if (i != Nv)
		throw "bPolygon::ReArrangeVertices()";
	for (int i = 0; i < Nv; ++i){
		int ip1 = (i != Nv - 1) ? i + 1 : 0;
		int im1 = i ? i - 1 : Nv - 1;
		vArray[i]->next = vArray[ip1];
		vArray[i]->prev = vArray[im1];
	}
}
bool bPolygon::isEar(Edge* e)//O(n)
{
	Vertex* v = e->GetDest();
	Vertex* vp = e->GetOrig();
	Vertex* vn = e->GetNextLeft()->GetDest();
	Vector3D A = vp->GetPoint();
	Vector3D B = v->GetPoint();
	Vector3D C = vn->GetPoint();
	double theta = Vector3D::RotationAngle_z(A, B, C);
	if (theta <= 0)
		return false;
	LineSegment2D AC(A, C);
	Vertex* vpp = e->GetPrevLeft()->GetOrig();
	Vector3D D = vpp->GetPoint();
	if (AC.GetRelationTo(D) != ON_LEFT)
		return false;
	EdgeIterator ite(this->inFace);
	Edge* ee = ite.Next();
	while (ee){
		Vertex* vO = ee->GetOrig();
		Vertex* vD = ee->GetDest();
		Vector3D pO = vO->GetPoint();
		Vector3D pD = vD->GetPoint();
		LineSegment2D LE(pO, pD);
		if (AC.isIntersecting(LE))
			return false;
		ee = ite.Next();
	}
	return true;
}
Edge* bPolygon::FindAnEar()//O(n^2)
{
	EdgeIterator ite(this->inFace);
	Edge* e = ite.Next();
	while (e){
		bool is_Ear = this->isEar(e);//O(n)
		if (is_Ear)
			return e;
		e = ite.Next();
	}
	return 0;
}
bPolygon::POINT_POLYGON_RELATION bPolygon::GetRelationTo(const Vector3D& p) const//O(n)
{
	EdgeIterator ite(this->inFace);
	Edge* e = ite.Next();
	int numIntersecting = 0;
	while (e)
	{
		Vertex* vO = e->GetOrig();
		Vertex* vD = e->GetDest();
		Vector3D pO = vO->GetPoint();
		Vector3D pD = vD->GetPoint();
		LineSegment2D LE(pO, pD);
		if (LE.isOn(p))
			return POINT_POLYGON_RELATION::ON_POLYGON;
		if (LE.isIntersectingRayFrom(p,1.356844))
		{
			numIntersecting += 1;
			if (fabs(p(1) - pO(1)) < 1.0e-15 || fabs(p(1) - pD(1)) < 1.0e-15)
			{
				Vector3D mOD = (pO + pD) / 2.0;
				if(mOD(1) > p(1) - 1.0e-15)
					numIntersecting -= 1;
			}
		}
		e = ite.Next();
	}
	if (numIntersecting % 2 == 1)
		return POINT_POLYGON_RELATION::INSIDE_POLYGON;
	return POINT_POLYGON_RELATION::OUTSIDE_POLYGON;
}
bPolygon::bPolygon()
{
	this->qe = 0;
	this->closed = false;
	this->latestVertex = 0;
	this->firstVertex = 0;
	this->inFace = 0;
	this->outFace = 0;
}
bPolygon::bPolygon(const bPolygon& rhs)
{
	this->CopyBody(rhs);
}
bPolygon::~bPolygon()
{
	this->DeleteBody();
}
void bPolygon::operator=(const bPolygon& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
}
QuadEdge* bPolygon::GetQuadEdgePtr() const
{
	QuadEdge* qeCopy = 0;
	if(this->qe)
		qeCopy = new QuadEdge(*this->qe);
	return qeCopy;
}
void bPolygon::AddVertex(const Vector3D& p)
{
	if (this->closed)
		throw "bPolygon::AddVertex(): Addding vertex to a closed Polygon.";
	if (p(2))
		throw "bPolygon::AddVertex(): 3D points are not accepted.";
	if (!this->qe)
		this->qe = new QuadEdge;
	Vertex* v = this->qe->AddVertex(p);
	if (!this->firstVertex)
		this->firstVertex = v;
	if (this->latestVertex)
	{
		Edge* e = this->qe->AddEdge();
		this->qe->AddEdgeToVertex(v, e, true);
		this->qe->AddEdgeToVertex(this->latestVertex, e, false);
	}
	this->latestVertex = v;
}
void bPolygon::AddVerticesAndClose(const vector<Vector3D>& points)
{
	for (int i = 0; i < points.size(); ++i)
	{
		this->AddVertex(points[i]);
	}
	this->Close();
}
void bPolygon::Close()
{
	if (!this->qe || this->qe->NumVertices() < 3)
		throw "Polygon::Close(): Not Enough Vertices";
	if (!this->closed)
	{
		Edge* eLast = this->qe->AddEdge(); 
		this->qe->AddEdgeToVertex(this->latestVertex, eLast, false);
		this->qe->AddEdgeToVertex(this->firstVertex, eLast, true);
		this->closed = true;
		this->inFace = this->qe->AddFace();
		this->outFace = this->qe->AddFace();
		bool ccw = this->isCCW();
		EdgeIterator ite(this->qe);
		Edge* e = ite.Next();
		while (e){
			this->qe->AddEdgeToFace(this->inFace, e, ccw);//assuming ccw
			this->qe->AddEdgeToFace(this->outFace, e, !ccw);//assuming ccw
			e = ite.Next();
		}
		if (!ccw)
			this->Reverse();
		this->ReArrangeVertices();
		this->UpdateArrows();
	}
}
Face* bPolygon::GetInsideFace()
{
	return this->inFace;
}
bool bPolygon::isEar(Vertex* v)
{
	Edge* e = this->Arrows[v];
	return this->isEar(e);
}
Triangulation bPolygon::Triangulate() const//O(n^3)
{
	Triangulation T;
	if (this->closed){
		int Nv = this->qe->NumVertices();
		bPolygon P(*this);
		while (Nv > 3)
		{
			Edge* e = P.FindAnEar();
			if (!e)
				throw "bPolygon::Triangulate() : CANNOT FIND AN EAR!!";
			Vertex* v = e->GetDest();
			Vertex* vp = e->GetOrig();
			Vertex* vn = e->GetNextLeft()->GetDest();
			Edge* eNew = P.qe->DivideFaceOnRight(P.inFace, vp, vn);
			--Nv;
		}
		T << *P.qe;
	}
	return T;
}
bool bPolygon::isInside(const Vector3D& p) const//O(n)
{
	return (this->GetRelationTo(p) == POINT_POLYGON_RELATION::INSIDE_POLYGON);
}
bool bPolygon::isOn(const Vector3D& p) const//O(n)
{
	return (this->GetRelationTo(p) == POINT_POLYGON_RELATION::ON_POLYGON);
}
bool bPolygon::isOnOrInside(const Vector3D& p) const//O(n)
{
	POINT_POLYGON_RELATION relationToP = this->GetRelationTo(p);
	return (relationToP == POINT_POLYGON_RELATION::ON_POLYGON 
		|| relationToP == POINT_POLYGON_RELATION::INSIDE_POLYGON);
}
bool bPolygon::TestIntegrity()
{
	if (this->closed)
	{
		EdgeIterator ite(this->inFace);
		Edge* e = ite.Next();
		while (e)
		{
			LineSegment2D L(e->GetOrig()->GetPoint(), e->GetDest()->GetPoint());
			EdgeIterator ite2(this->inFace);
			Edge* e2 = ite2.Next();
			while (e2)
			{
				LineSegment2D L2(e2->GetOrig()->GetPoint(), e2->GetDest()->GetPoint()); 
				if (e2 != e && L.isIntersecting(L2))
					return false;
				e2 = ite2.Next();
			}
			e = ite.Next();
		}
		if (this->qe->NumVertices() > 3)
		{
			Edge* ear = this->FindAnEar();
			if (!ear)
				return false;
		}
	}
	return true;
}
bool tester_bPolygon(int& NumTests)
{
	if (!tester_bPolygon_1(NumTests))
		return false;
	if (!tester_bPolygon_2(NumTests))
		return false;
	if (!tester_bPolygon_3(NumTests))
		return false;
	if (!tester_bPolygon_4(NumTests))
		return false;
	if (!tester_bPolygon_5(NumTests))
		return false;
	if (!tester_bPolygon_6(NumTests))
		return false;
	if (!tester_bPolygon_7(NumTests))
		return false;
	if (!tester_bPolygon_8(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_1(int& NumTests)
{
	Vector3D A, B(1), C(2, 2), D(0, 1);
	bPolygon PPP;
	PPP.AddVertex(A);
	PPP.AddVertex(B);
	PPP.AddVertex(C);
	PPP.AddVertex(D);
	PPP.Close();
	if (!PPP.TestIntegrity())
		return false;
	bPolygon PP(PPP), P;
	P = PP;
	Face* inface = P.GetInsideFace();
	if (!inface)
		return false;
	VertexIterator itv(inface);
	EdgeIterator ite(inface);
	const int N = 4;
	Vertex* v[N];
	Edge* e[4];
	for (int i = 0; i < N; ++i){
		v[i] = 0;
		v[i] = itv.Next();
		e[i] = 0;
		e[i] = ite.Next();
		if (e[i]->GetLeft() != inface)
			return false;
	}
	if (itv.Next() || !v[N - 1])
		return false;
	if (ite.Next() || !e[N - 1])
		return false;
	if (v[0]->GetPoint() != D || v[1]->GetPoint() != A || v[2]->GetPoint() != B || v[3]->GetPoint() != C)
		return false;
	if (!P.isEar(v[0]) || !P.isEar(v[1]) || !P.isEar(v[2]) || !P.isEar(v[3])) 
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	Vector3D E(0.5, 0.5), F(3, 3), G(3, 0), H(0.5, 0);
	if (!P.isInside(E) || !P.isOnOrInside(E))
		return false;
	if (P.isInside(F) || P.isOnOrInside(F))
		return false;
	if (P.isInside(G) || P.isOnOrInside(G))
		return false;
	if (P.isInside(H) || !P.isOnOrInside(H))
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_2(int& NumTests)
{
	Vector3D A, B(1), C(2, 2), D(0, 1);
	bPolygon PPP;
	PPP.AddVertex(A);
	PPP.AddVertex(D);
	PPP.AddVertex(C);
	PPP.AddVertex(B);
	PPP.Close();
	if (!PPP.TestIntegrity())
		return false;
	bPolygon PP(PPP), P;
	P = PP;
	Face* inface = P.GetInsideFace();
	if (!inface)
		return false;
	VertexIterator itv(inface);
	EdgeIterator ite(inface);
	const int N = 4;
	Vertex* v[N];
	Edge* e[4];
	for (int i = 0; i < N; ++i){
		v[i] = 0;
		v[i] = itv.Next();
		e[i] = 0;
		e[i] = ite.Next();
		if (e[i]->GetLeft() != inface)
			return false;
	}
	if (itv.Next() || !v[N - 1])
		return false;
	if (ite.Next() || !e[N - 1])
		return false;
	if (v[0]->GetPoint() != A || v[1]->GetPoint() != B || v[2]->GetPoint() != C || v[3]->GetPoint() != D)
		return false;
	if (!P.isEar(v[0]) || !P.isEar(v[1]) || !P.isEar(v[2]) || !P.isEar(v[3]))
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_3(int& NumTests)
{
	Vector3D A, B(3,3), C(0, 1), D(-3, 3);
	bPolygon PPP;
	PPP.AddVertex(A);
	PPP.AddVertex(D);
	PPP.AddVertex(C);
	PPP.AddVertex(B);
	PPP.Close();
	if (!PPP.TestIntegrity())
		return false;
	bPolygon PP(PPP), P;
	P = PP;
	Face* inface = P.GetInsideFace();
	if (!inface)
		return false;
	VertexIterator itv(inface);
	EdgeIterator ite(inface);
	const int N = 4;
	Vertex* v[N];
	Edge* e[4];
	for (int i = 0; i < N; ++i){
		v[i] = 0;
		v[i] = itv.Next();
		e[i] = 0;
		e[i] = ite.Next();
		if (e[i]->GetLeft() != inface)
			return false;
	}
	if (itv.Next() || !v[N - 1])
		return false;
	if (ite.Next() || !e[N - 1])
		return false;
	if (v[0]->GetPoint() != A || v[1]->GetPoint() != B || v[2]->GetPoint() != C || v[3]->GetPoint() != D)
		return false;
	if (P.isEar(v[0]) || !P.isEar(v[1]/*B*/) || P.isEar(v[2]/*C*/) || !P.isEar(v[3]/*D*/))
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_4(int& NumTests)
{
	vector<Vector3D> pentagonPoints;
	pentagonPoints.resize(5);
	bPolygon Pentagon;
	for (int i = 0; i < 5; ++i){
		double theta = double(i) * 72.0 / 180.0 * pi + pi / 2.0;
		double r = 1.0;
		pentagonPoints[i](0) = r * cos(theta);
		pentagonPoints[i](1) = r * sin(theta);
		Pentagon.AddVertex(pentagonPoints[i]);
	}
	Pentagon.Close();
	if (!Pentagon.TestIntegrity())
		return false;
	//bPolygon PP(PPP), Pentagon;
	//Pentagon = PP;
	Triangulation T = Pentagon.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	int N = 5;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_5(int& NumTests)
{
	vector<Vector3D> hexagonPoints;
	int N = 6;
	hexagonPoints.resize(N);
	bPolygon Hexagon;
	for (int i = 0; i < N; ++i){
		double theta = double(i) * (2 * pi) / double(N);
		double r = 1.0;
		hexagonPoints[i](0) = r * cos(theta);
		hexagonPoints[i](1) = r * sin(theta);
		Hexagon.AddVertex(hexagonPoints[i]);
	}
	Hexagon.Close();
	if (!Hexagon.TestIntegrity())
		return false;
	//bPolygon PP(PPP), Pentagon;
	//Pentagon = PP;
	Triangulation T = Hexagon.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_6(int& NumTests)
{
	vector<Vector3D> points;
	points.resize(11);
	points[0](0) = 0; points[0](1) = 0;
	points[1](0) = 1; points[1](1) = 0;
	points[2](0) = 2; points[2](1) = 1;
	points[3](0) = 2; points[3](1) = 2;
	points[4](0) = 3; points[4](1) = 2;
	points[5](0) = 3; points[5](1) = 3;
	points[6](0) = 0.5; points[6](1) = 3;
	points[7](0) = 0; points[7](1) = 1.5;
	points[8](0) = -1; points[8](1) = 1.5;
	points[9](0) = -1; points[9](1) = 0.75;
	points[10](0) = 0; points[10](1) = 0.75;
	bPolygon P;
	for (int i = 0; i < 11; ++i){
		P.AddVertex(points[i]);
	}
	P.Close();
	if (!P.TestIntegrity())
		return false;
	//bPolygon PP(PPP), Pentagon;
	//Pentagon = PP;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	int N = 11;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_7(int& NumTests)
{
	vector<Vector3D> points;
	int Nx = 2;
	int Ny = 2;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	points.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		points[i](0) = hx * i; points[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		points[n](0) = Lx; points[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		points[n](0) = Lx - i * hx; points[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		points[n](0) = 0; points[n](1) = Ly - i * hy;
	}
	int N = points.size();
	bPolygon P;
	for (int i = 0; i < N; ++i){
		P.AddVertex(points[i]);
	}
	P.Close();
	if (!P.TestIntegrity())
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}
bool tester_bPolygon_8(int& NumTests)
{
	double h = sqrt(3) / 2.0;
	vector<Vector3D> points;
	points.resize(27);
	points[0](0) = 6;		points[0](1) = 5;
	points[1](0) = 6;		points[1](1) = 6;
	points[2](0) = 5;		points[2](1) = 6;
	points[3](0) = 4;		points[3](1) = 6;
	points[4](0) = 3;		points[4](1) = 6;
	points[5](0) = 2;		points[5](1) = 6;
	points[6](0) = 1;		points[6](1) = 6;
	points[7](0) = 0;		points[7](1) = 6;
	points[8](0) = 0;		points[8](1) = 5;
	points[9](0) = 0;		points[9](1) = 4;
	points[10](0) = 0;		points[10](1) = 3;
	points[11](0) = 0;		points[11](1) = 2;
	points[12](0) = 0;		points[12](1) = 1;
	points[13](0) = 1;		points[13](1) = 0;
	points[14](0) = 1.5;	points[14](1) = h;
	points[15](0) = 2;		points[15](1) = 0;
	points[16](0) = 2.5;	points[16](1) = h;
	points[17](0) = 3;		points[17](1) = 0;
	points[18](0) = 3.5;	points[18](1) = h;
	points[19](0) = 4;		points[19](1) = 0;
	points[20](0) = 4.5;	points[20](1) = h;
	points[21](0) = 6;		points[21](1) = 2;
	points[22](0) = 6-h;	points[22](1) = 2.5;
	points[23](0) = 6;		points[23](1) = 3;
	points[24](0) = 6-h;	points[24](1) = 3.5;
	points[25](0) = 6;		points[25](1) = 4;
	points[26](0) = 6-h;	points[26](1) = 4.5;
	int N = points.size();
	bPolygon P;
	for (int i = 0; i < N; ++i){
		P.AddVertex(points[i]);
	}
	P.Close();
	if (!P.TestIntegrity())
		return false;
	Triangulation T = P.Triangulate();
	//T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (T.GetMesh2D()->NumVertices() != N || T.GetMesh2D()->NumFaces() != N - 1 || T.GetMesh2D()->NumEdges() != 2 * N - 3)
		return false;
	++NumTests;
	return true;
}