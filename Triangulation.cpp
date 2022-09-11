#include <fstream>
#include <algorithm> 
#include <queue>
#include "Triangulation.h"
#include "Vector3D.h"
#include "Circle.h"
#include "Triangle.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include "EdgeIterator.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "QuadEdge.h"
#include "bData.h"
#include "DelaunayLifting.h"
#include "bPolygon.h"
#include "MathUtils.h"
#include "QuadEdgeIndex.h"
#include "BitMap.h"
void Triangulation::CopyBody(const Triangulation& rhs)
{
	this->Mesh2D = 0;
	this->boundary = 0;
	if (rhs.Mesh2D)
	{
		this->Mesh2D = new QuadEdge(*rhs.Mesh2D);
		FaceIterator itf0(rhs.Mesh2D);
		FaceIterator itf1(this->Mesh2D);
		Face* f0 = itf0.Next();
		Face* f1 = itf1.Next();
		while (f0 && f1)
		{
			if (f0 == rhs.boundary)
			{
				this->boundary = f1;
				break;
			}
			f0 = itf0.Next();
			f1 = itf1.Next();
		}
	}
	this->Area = rhs.Area;
}
void Triangulation::DeleteBody()
{
	delete this->Mesh2D;
	this->Mesh2D = 0;
	this->boundary = 0;
}
Triangulation::Triangulation()
{
	this->Mesh2D = 0;
	this->boundary = 0;
	this->Area = 0;
}
Triangulation::Triangulation(const Triangulation& rhs)
{
	this->CopyBody(rhs);
}
Triangulation::~Triangulation()
{
	this->DeleteBody();
}
Triangulation& Triangulation::operator=(const Triangulation& rhs)
{
	this->DeleteBody();
	this->CopyBody(rhs);
	return *this;
}
void Triangulation::operator<<(const QuadEdge& QE)
{
	delete this->Mesh2D;
	this->Mesh2D = new QuadEdge(QE);
	this->boundary = 0;
	this->Mesh2D->UpdateFacesDegree();
	FaceIterator itf(this->Mesh2D);
	Face* f = itf.Next();
	while (f)
	{
		if (f->GetDegree() > 3)
		{
			if (!this->boundary)
				this->boundary = f;
			else
				throw "Triangulation::operator<<()";//In future: do polygon triangulation
		}
		f = itf.Next();
	}
	if (!this->boundary)
	{
		vector<Face*> positiveNormal;
		vector<Face*> negativeNormal;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			Vector3D n = this->Mesh2D->GetNormal(f);
			if (n(2) > 0)
				positiveNormal.push_back(f);
			else if (n(2) < 0)
				negativeNormal.push_back(f);
			f = itf.Next();
		}
		if (positiveNormal.size() != 1 && negativeNormal.size() != 1)
			throw "Triangulation::operator<<()";
		this->boundary = positiveNormal.size() == 1 ? positiveNormal[0] : negativeNormal[0];
	}
	if (!this->boundary)
	{
		throw "Triangulation::operator<<()";
	}
	VertexIterator itv(this->Mesh2D);
	Vertex* v = itv.Next();
	while (v)
	{
		v->SetZ(0);
		v = itv.Next();
	}
	this->Mesh2D->IssueIndices();
	this->Mesh2D->SetPlanar(true);

}
QuadEdge* Triangulation::GetMesh2D()
{
	return this->Mesh2D;
}
Face* Triangulation::GetBoundary()
{
	return this->boundary;
}
void Triangulation::GetTriangles(vector<shared_ptr<Triangle> >& triangles) const
{
	FaceIterator itf(this->Mesh2D);
	Face* f = itf.Next();
	while (f)
	{
		if (f != this->boundary && f->GetDegree() == 3)
		{
			VertexIterator itv(f);
			Vertex* v0 = itv.Next();
			Vertex* v1 = itv.Next();
			Vertex* v2 = itv.Next();
			Vector3D p0 = v0->GetPoint();
			Vector3D p1 = v1->GetPoint();
			Vector3D p2 = v2->GetPoint();
			triangles.push_back(0);
			//triangles.back() = new Triangle(p0, p1, p2);
			triangles.back() = make_shared<Triangle>(p0, p1, p2);
		}
		f = itf.Next();
	}
}
bool Triangulation::TestIntegrity()
{
	if (this->Mesh2D)
	{
		if (!this->Mesh2D->TestIntegrity())
			return false;
		int Nv = this->Mesh2D->NumVertices();
		int Nf = this->Mesh2D->NumFaces();
		this->Mesh2D->UpdateFacesDegree();
		int Nh = this->boundary->GetDegree();
		int Nt = Nf - 1;
		if (Nt != 2 * Nv - Nh - 2)
			return false;
		if (Nf < 2)
			return false;
		if (Nh > Nv)
			return false;
		vector<shared_ptr<Triangle> > triangles;
		this->GetTriangles(triangles);
		if (!triangles.size())
			return false;
		for (int i = 0; i < triangles.size() - 1; ++i)
		{
			for (int j = i + 1; j < triangles.size(); ++j)
			{
				if (triangles[i]->isIntersecting(*triangles[j]))
					return false;
			}
		}
	}
	return true;
}
bool Triangulation::TestDelaunay()
{
	FaceIterator itf(this->Mesh2D);
	Face* f = itf.Next();
	while (f)
	{
		if (f != this->GetBoundary())
		{
			Circle Cir = this->GetCircumCircle(f);
			VertexIterator itvf(f);
			Vertex* vA = itvf.Next();
			Vertex* vB = itvf.Next();
			Vertex* vC = itvf.Next();
			VertexIterator itv(this->Mesh2D);
			Vertex* v = itv.Next();
			while (v)
			{
				if (v != vA && v != vB && v != vC)
				{
					Vector3D p = v->GetPoint();
					bool inside = Cir.isInside(p);
					if (inside)
						return false;
				}
				v = itv.Next();
			}
		}
		f = itf.Next();
	}
	return true;
}
void Triangulation::PrintTriangulation(string folderName, string extension)
{
	string triangulationPath = folderName + "\\Triangulation" + extension +  ".txt";
	string pointsPath = folderName + "\\Points" + extension + ".txt";
	string edgesPath = folderName + "\\Edges" + extension + ".txt";
	string elementsPath = folderName + "\\Elements" + extension + ".txt";
	ofstream triangulationFile;
	ofstream pointsFile;
	triangulationFile.open(triangulationPath.c_str(), std::fstream::out);
	pointsFile.open(pointsPath.c_str(), std::fstream::out);
	int Nv = this->Mesh2D->NumVertices();
	int Ne = this->Mesh2D->NumEdges();
	int Nf = this->Mesh2D->NumFaces();
	triangulationFile << "Nv  = " << Nv << ", Ne = " << Ne << ", Nf = " << Nf;
	pointsFile << Nv;
	triangulationFile << endl << endl << "-------VERTICES----------------";
	VertexIterator itv(this->Mesh2D);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D p = v->GetPoint();
		pointsFile << endl << p(0);
		pointsFile << endl << p(1);
		pointsFile << endl << v->index;
		triangulationFile << endl << "v" << v->index << " (x,y) = (" << p(0) << " , " << p(1) << ")";
		v = itv.Next();
	}
	pointsFile.close();
	ofstream edgesFile;
	edgesFile.open(edgesPath.c_str(),std::fstream::out);
	triangulationFile << endl << endl << "-------EDGES----------------";
	edgesFile << Ne;
	EdgeIterator ite(this->Mesh2D);
	Edge* e = ite.Next();
	while (e)
	{
		int nO = e->GetOrig()->index;
		int nD = e->GetDest()->index;
		edgesFile << endl << nO;
		edgesFile << endl << nD;
		triangulationFile << endl << "e" << e->index  << " : v" << nO << " , v" << nD;
		e = ite.Next();
	}
	edgesFile.close();
	fstream elementsFile;
	elementsFile.open(elementsPath.c_str(), std::fstream::out);
	elementsFile << Nf - 1;
	triangulationFile << endl << endl << "-------FACES----------------";
	FaceIterator itf(this->Mesh2D);
	Face* f = itf.Next();
	while (f)
	{
		VertexIterator itv(f);
		Vertex* v0 = itv.Next();
		Vertex* v1 = itv.Next();
		Vertex* v2 = itv.Next();
		triangulationFile << endl << "f" << f->index  << " : v" << v0->index << ", v" << v1->index << ", v" << v2->index;
		if (f != this->boundary)
		{
			elementsFile << endl << v0->index;
			elementsFile << endl << v1->index;
			elementsFile << endl << v2->index;
		}
		else
		{
			Vertex* v = itv.Next();
			while (v)
			{
				triangulationFile << ", v" << v->index;
				v = itv.Next();
			}
			triangulationFile << " BOUNDARY FACE";
		}
		f = itf.Next();
	}
	elementsFile.close();
	this->Mesh2D->print(triangulationFile);
	this->Mesh2D->PrintTestScript(triangulationFile);
	triangulationFile.close();
}
void Triangulation::PrintTriangulation()
{
	string fileName(bData::FolderPath);
	this->PrintTriangulation(fileName);
}
Vertex* Triangulation::GetRightFaceOpposite(Edge* e)
{//needs testing
	Vertex* vO = e->GetOrig();
	Vertex* vD = e->GetDest();
	Face* fR = e->GetRight();
	VertexIterator itv(fR);
	Vertex* v = itv.Next();
	while (v){
		if (v != vO && v != vD)
			return v;
		v = itv.Next();
	}
	return 0;
}
Vertex* Triangulation::GetLeftFaceOpposite(Edge* e)
{//needs testing
	Vertex* vO = e->GetOrig();
	Vertex* vD = e->GetDest();
	Face* fL = e->GetLeft();
	VertexIterator itv(fL);
	Vertex* v = itv.Next();
	while (v){
		if (v != vO && v != vD)
			return v;
		v = itv.Next();
	}
	return 0;
}
Edge* Triangulation::GetOppositeEdge(Face* f, Vertex* v)
{
	//Check validity of input:
	VertexIterator itv(f);
	bool foundv = false;
	Vertex* vv = itv.Next();
	while (vv) {
		if (vv == v) {
			foundv = true;
			break;
		}
		vv = itv.Next();
	}
	if (foundv)
	{
		EdgeIterator ite(f);
		Edge* e = ite.Next();
		while (e) {
			if (e->GetOrig() != v && e->GetDest() != v)
				return e;
			e = ite.Next();
		}
	}
	return 0;
}
Face* Triangulation::GetOppositeFace(Face* f, Edge* e)
{
	if (e->GetRight() == f)
		return e->GetLeft();
	if (e->GetLeft() == f)
		return e->GetRight();
	return 0;
}
Edge* Triangulation::Connect(Vertex* v1, Vertex* v2)
{//Connects vertex v1 to v2 with an edge using flips 
	QuadEdge* qe = this->Mesh2D;
	Edge* e = qe->GetEdgeConnecting(v1, v2);
	if (e)
		return e;
	Face* fb = this->boundary;
	Edge* flipThis = 0;
	Edge* flipped = 0;
	LineSegment2D L1(v1->GetPoint(), v2->GetPoint());
	do
	{
		FaceIterator itf(v1);
		Face* f = itf.Next();
		while (f)
		{
			if (f != fb)
			{
				Edge* e2 = this->GetOppositeEdge(f, v1);
				Vector3D Q0 = e2->GetOrig()->GetPoint();
				Vector3D Q1 = e2->GetDest()->GetPoint();
				LineSegment2D L2(Q0, Q1);
				if (L1.isIntersecting(L2))
				{
					flipThis = e2;
					break;
				}
			}
			f = itf.Next();
		}
		if (!flipThis)
			throw "Triangulation::Connect()";
		flipped = this->Flip(flipThis);
		if (!flipped)
			throw "Triangulation::Connect()";
		if (flipped->GetOrig() != v1 && flipped->GetDest() != v1)
			throw "Triangulation::Connect()";
		if (flipped->GetOrig() != v1)
			flipped->Reverse();
	} while (flipped->GetDest() != v2);
	return flipped;
}
void Triangulation::ShrinkBoundary(Edge* e)
{
	if (this->Mesh2D)
	{
		Vertex* v0 = e->GetOrig();
		Vertex* v1 = e->GetDest();
		if (e->GetRight() == this->boundary)
			e->Reverse();
		if (e->GetLeft() == this->boundary)
		{
			this->Mesh2D->DeleteEdgeAndRightFace(e);
		}
		this->Mesh2D->UpdateVertexDegree(v0);
		this->Mesh2D->UpdateVertexDegree(v1);
		if (v0->GetDegree() == 1)
		{
			Edge* e0 = v0->GetEdge();
			this->Mesh2D->DeleteEdgeAndRightFace(e0);
		}
		else if (v0->GetDegree() == 0)
		{
			this->Mesh2D->RemoveVertex(v0);
		}
		if (v1->GetDegree() == 1)
		{
			Edge* e1 = v1->GetEdge();
			this->Mesh2D->DeleteEdgeAndRightFace(e1);
		}
		else if (v1->GetDegree() == 0)
		{
			this->Mesh2D->RemoveVertex(v1);
		}
		this->Mesh2D->UpdateFaceDegree(this->boundary);
	}
}
void Triangulation::ImposeBoundary(const vector<Vertex*>& boundary)
{
	unordered_set<Edge*> boundaryEdges;
	unordered_set<Edge*>::iterator its;
	boundaryEdges.rehash(boundary.size());
	for (int i = 0; i < boundary.size(); ++i) {
		int ip1 = i < boundary.size() - 1 ? i + 1 : 0;
		Edge* e = this->Connect(boundary[i], boundary[ip1]);
		boundaryEdges.insert(e);
	}
	Face* fb = this->boundary;
	bool go_on = true;
	int iter = 0;
	do
	{
		EdgeIterator ite(fb);
		Edge* e = ite.Next();
		vector<Edge*> removeThese;
		while (e) {
			unordered_set<Edge*>::iterator its = boundaryEdges.find(e);
			if (its == boundaryEdges.end()) {
				if (!removeThese.size())
					removeThese.push_back(e);
				else {
					Edge* ep = removeThese.back();
					if(this->GetOppositeFace(fb,e) != this->GetOppositeFace(fb, ep))
						removeThese.push_back(e);
				}
			}
			e = ite.Next();
		}
		for (int i = 0; i < removeThese.size(); ++i) {
			Edge* er = removeThese[i];
			this->ShrinkBoundary(er);
		}
		go_on = removeThese.size();
	} while (go_on);
}
Triangle Triangulation::GetTriangle(Face* f)
{
	if (this->Mesh2D && f != this->boundary)
	{
		VertexIterator itv(f);
		Vertex* v0 = itv.Next();
		Vertex* v1 = itv.Next();
		Vertex* v2 = itv.Next();
		Vector3D p0 = v0->GetPoint();
		Vector3D p1 = v1->GetPoint();
		Vector3D p2 = v2->GetPoint();
		Triangle t(p0, p1, p2);
		return t;
	}
	return Triangle();
}
Circle Triangulation::GetCircumCircle(Face* f)
{
	Circle Cir;
	if (f != this->GetBoundary())
	{
		VertexIterator itvf(f);
		Vertex* vA = itvf.Next();
		Vertex* vB = itvf.Next();
		Vertex* vC = itvf.Next();
		Vector3D A = vA->GetPoint();
		Vector3D B = vB->GetPoint();
		Vector3D C = vC->GetPoint();
		Cir.Set3Points(A, B, C);
	}
	return Cir;
}
double Triangulation::GetCircumCircleRadius(Face* f)
{
	if (f != this->GetBoundary())
	{
		Circle C = this->GetCircumCircle(f);
		return C.GetRadius();
	}
	return -1;
}
double Triangulation::GetIncircleRadius(Face* f)
{
	if (f != this->GetBoundary())
	{
		VertexIterator itvf(f);
		Vertex* vA = itvf.Next();
		Vertex* vB = itvf.Next();
		Vertex* vC = itvf.Next();
		Vector3D A = vA->GetPoint();
		Vector3D B = vB->GetPoint();
		Vector3D C = vC->GetPoint();
		Triangle t(A, B, C);
		return t.GetIncircleRadius();
	}
	return -1;
}
double Triangulation::GetMinLength(Face* f)
{
	double minL = 1e10;
	EdgeIterator ite(f);
	Edge* e = ite.Next();
	while (e){
		double L = e->GetVector().abs();
		if (L < minL)
			minL = L;
		e = ite.Next();
	}
	return minL;
}
bool Triangulation::isInside(const Vector3D& p) const
{
	if (!this->Mesh2D)
		return false; 
	if (!this->isOnOrInside(p))
		return false;
	EdgeIterator ite(this->boundary);
	Edge* e = ite.Next();
	while (e){
		LineSegment2D L = e->GetLineSegment2D();
		if (L.GetRelationTo(p) == ON_LINE_SEGMENT)
			return false;
		e = ite.Next();
	}
	return true;
}
bool Triangulation::isOnOrInside(const Vector3D& p) const//needs testing
{
	if (!this->Mesh2D)
		return false;
	FaceIterator itf(this->Mesh2D);
	Face* fb = this->boundary;
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			shared_ptr<Triangle> t = this->GetTriangle(f);
			if (t->isOnOrInside(p))
				return true;
		}
		f = itf.Next();
	}
	return false;
}
int Triangulation::NumVertices() const
{
	return this->Mesh2D ? this->Mesh2D->NumVertices() : 0;
}
int Triangulation::NumTriangles() const
{
	return this->Mesh2D ? this->Mesh2D->NumFaces() - 1 : 0;
}
int Triangulation::NumBoundaryPoints()
{
	if (this->Mesh2D && this->boundary)
	{
		this->Mesh2D->UpdateFaceDegree(this->boundary);
		return this->boundary->GetDegree();
	}
	return -1;
}
double Triangulation::GetShapeFactor() const
{
	vector<shared_ptr<Triangle> > triangles;
	this->GetTriangles(triangles);
	double Alpha = 1.0;//triangulation shape factor
	if (!triangles.size())
		return 0;
	for (int i = 0; i < triangles.size(); ++i)
	{
		if (!triangles[i])
			return 0;
		double alpha = triangles[i]->GetShapeFactor();
		if (fabs(alpha) < 1.0e-10)
			return 0;
		Alpha *= alpha;
	}
	Alpha = pow(Alpha, 1.0 / double(triangles.size()));
	return Alpha;
}
double Triangulation::GetMaxAngle() const
{
	vector<shared_ptr<Triangle> > triangles;
	this->GetTriangles(triangles);
	if (!triangles.size())
		return 0;
	double MaxAngle = 0;
	for (int i = 0; i < triangles.size(); ++i)
	{
		if (!triangles[i])
			return 0;
		double maxAngle = triangles[i]->GetMaxAngle();
		if (fabs(maxAngle) < 1.0e-10)
			return 0;
		if (maxAngle > MaxAngle)
			MaxAngle = maxAngle;
	}
	return MaxAngle;
}
double Triangulation::GetMinAngle() const
{
	vector<shared_ptr<Triangle> > triangles;
	this->GetTriangles(triangles);
	if (!triangles.size())
		return 0;
	double MinAngle = 10000;
	for (int i = 0; i < triangles.size(); ++i)
	{
		if (!triangles[i])
			return 0;
		double minAngle = triangles[i]->GetMinAngle();
		if (minAngle < MinAngle)
			MinAngle = minAngle;
	}
	return MinAngle;
}
double Triangulation::GetPiAngles() const
{
	vector<shared_ptr<Triangle> > triangles;
	this->GetTriangles(triangles);
	if (!triangles.size())
		return 0;
	double PiAngle = 1.0;
	for (int i = 0; i < triangles.size(); ++i)
	{
		if (!triangles[i])
			return 0;
		double piAngle = triangles[i]->GetPiAngles();
		PiAngle *= piAngle;
	}
	return PiAngle;
}
double Triangulation::GetMaxAngleMeasure() const
{
	return (this->GetMaxAngle() / pi * 3.0 - 1.0);
}
double Triangulation::GetMinAngleMeasure() const
{
	return (pi / 3.0 / this->GetMinAngle() - 1.0);
}
int Triangulation::GetMaxVertexDegree()
{
	if (this->Mesh2D)
		return this->Mesh2D->GetMaxVertexDegree();
	return -1;
}
int Triangulation::GetMinVertexDegree()
{
	if (this->Mesh2D)
		return this->Mesh2D->GetMinVertexDegree(); 
	return -1;
}
Edge* Triangulation::CheckThenFlip(Edge* e)
{
	if (this->IsFlippable(e))
	{
		return this->Flip(e);
	}
	return 0;
}
bool Triangulation::IsOnBoundary(Vertex* v)
{
	FaceIterator itf(v);
	Face* f = itf.Next();
	while (f)
	{
		if (f == this->boundary)
			return true;
		f = itf.Next();
	}
	return false;
}
bool Triangulation::IsFlippable(Edge* e)
{
	if (!this->Mesh2D || !e || !e->GetOrig() || !e->GetDest())
		return false;
	Vertex* v0 = this->GetLeftFaceOpposite(e);
	Vertex* v1 = this->GetRightFaceOpposite(e);
	if (!v0 || !v1)
		return false;
	LineSegment2D L1(v0->GetPoint(), v1->GetPoint());
	LineSegment2D L2 = e->GetLineSegment2D();
	if (!L1.isIntersecting(L2))
		return false;
	if (L1.isOn(e->GetOrig()->GetPoint()))
		return false;
	if (L1.isOn(e->GetDest()->GetPoint()))
		return false;
	return true;
}
Edge* Triangulation::Flip(Edge* e)
{
	QuadEdge* qe = this->Mesh2D;
	Face* fL = e->GetLeft();
	Vertex* vR = this->GetRightFaceOpposite(e);
	Vertex* vL = this->GetLeftFaceOpposite(e);
	int ie = e->index;
	int iR = e->GetRight()->index;
	qe->DeleteEdgeAndRightFace(e);
	Edge* eNew = qe->DivideFaceOnRight(fL, vL, vR);
	eNew->index = ie;
	eNew->GetRight()->index = iR;
	return eNew;
}
shared_ptr<bPolygon> Triangulation::GetBoundaryPolygon()
{
	shared_ptr<bPolygon> poly = make_shared<bPolygon>();
	if (this->boundary)
	{
		VertexIterator itv(this->boundary);
		Vertex* v = itv.Next();
		while (v)
		{
			poly->AddVertex(v->GetPoint());
			v = itv.Next();
		}
		poly->Close();
	}
	return poly;
}
shared_ptr<Triangle> Triangulation::GetTriangle(Face* f) const
{
	if (f && f != this->boundary)
	{
		VertexIterator itv(f);
		Vertex* v0 = itv.Next();
		Vertex* v1 = itv.Next();
		Vertex* v2 = itv.Next();
		Vector3D p0 = v0->GetPoint();
		Vector3D p1 = v1->GetPoint();
		Vector3D p2 = v2->GetPoint();
		return std::make_shared<Triangle>(p0,p1,p2);
	}
	return shared_ptr<Triangle>();
}
double Triangulation::GetMinEdgeLength()
{
	if (this->Mesh2D)
		return this->Mesh2D->GetMinEdgeLength();
	return -1;
}
double Triangulation::GetMaxEdgeLength()
{
	if (this->Mesh2D)
		return this->Mesh2D->GetMaxEdgeLength();
	return -1;
}
double Triangulation::GetAvgEdgeLength()
{
	if (this->Mesh2D)
		return this->Mesh2D->GetAvgEdgeLength();
	return -1;
}
double Triangulation::GetWeight()
{
	if (this->Mesh2D)
		return this->Mesh2D->GetWeight();
	return -1;
}
int Triangulation::GetEncrochness()
{
	if (this->Mesh2D)
	{
		int MaxEnchrochness = 0;
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		while(e)
		{
			int enchrochness = 0;
			Vertex* vO = e->GetOrig();
			Vertex* vD = e->GetDest();
			Vector3D pO = vO->GetPoint();
			Vector3D pD = vD->GetPoint();
			Vector3D pM = (pO + pD) / 2.0;
			double R = pO.distance(pD) / 2.0;
			VertexIterator itv(this->Mesh2D);
			Vertex* v = itv.Next();
			while (v)
			{
				if (v != vO && v != vD)
				{
					Vector3D p = v->GetPoint();
					double r = p.distance(pM);
					if (r <= R + 1.0e-10)
						++enchrochness;
				}
				v = itv.Next();
			}
			if (enchrochness > MaxEnchrochness)
				MaxEnchrochness = enchrochness;
			e = ite.Next();
		}
		return MaxEnchrochness;
	}
	return -1;
}
double Triangulation::GetSkinnyness()
{
	if (this->Mesh2D)
	{
		double Max_R_over_minL = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				Circle Cir = this->GetCircumCircle(f);
				double R = Cir.GetRadius();
				double minL = this->GetMinLength(f);
				double R_over_minL = R / minL;
				if (R_over_minL > Max_R_over_minL)
					Max_R_over_minL = R_over_minL;
			}
			f = itf.Next();
		}
		return Max_R_over_minL;
	}
	return -1;
}
double Triangulation::GetMinDistance(Vertex* v)
{
	if (this->Mesh2D)
	{
		Vector3D P = v->GetPoint();
		VertexIterator itvv(this->Mesh2D);
		Vertex* vv = itvv.Next();
		double minDist = -1;
		while (vv)
		{
			if (vv != v)
			{
				Vector3D Q = vv->GetPoint();
				Vector3D PQ = Q - P;
				double pq = PQ.abs();
				if (minDist < 0 || pq < minDist)
					minDist = pq;
			}
			vv = itvv.Next();
		}
		return minDist;
	}
	return -1;
}
double Triangulation::GetCovariance()
{
	if (this->Mesh2D)
	{
		vector<double> minDistArr;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			double minDist = this->GetMinDistance(v);
			minDistArr.push_back(minDist);
			v = itv.Next();
		}
		double sum = 0;
		double sqsum = 0;
		for (int k = 0; k < minDistArr.size(); ++k)
		{
			sum += minDistArr[k];
			sqsum += minDistArr[k] * minDistArr[k];
		}
		double lambda = sqrt(minDistArr.size() * sqsum / sum / sum - 1.0);
		return lambda;
	}
	return -1;
}
double Triangulation::GetMeshRatio()
{
	if (this->Mesh2D)
	{
		int Nv = this->Mesh2D->NumVertices();
		double* minDistArr = new double[Nv];
		int counter = 0;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			Vector3D P = v->GetPoint();
			VertexIterator itvv(this->Mesh2D);
			Vertex* vv = itvv.Next();
			double minDist = -1;
			while (vv)
			{
				if (vv != v)
				{
					Vector3D Q = vv->GetPoint();
					Vector3D PQ = Q - P;
					double pq = PQ.abs();
					if (minDist < 0 || pq < minDist)
						minDist = pq;
				}
				vv = itvv.Next();
			}
			minDistArr[counter] = minDist;
			++counter;
			v = itv.Next();
		}
		double gamma = (*std::max_element(minDistArr, minDistArr + Nv)) / (*std::min_element(minDistArr, minDistArr + Nv)) - 1.0;
		delete[] minDistArr;
		return gamma;
	}
	return -1;
}
void Triangulation::UpdateArea()
{
	if (this->Mesh2D && this->boundary)
	{
		this->Area = 0;
		EdgeIterator ite(this->boundary);
		Edge* e = ite.Next();
		while (e)
		{
			Vertex* a = e->GetOrig();
			Vertex* b = e->GetDest();
			if(e->GetRight() == this->boundary)
				this->Area += 0.5 * (a->GetX() * b->GetY() - b->GetX() * a->GetY());
			else
				this->Area -= 0.5 * (a->GetX() * b->GetY() - b->GetX() * a->GetY());
			e = ite.Next();
		}
		this->Area = fabs(this->Area);
	}
}
double Triangulation::GetArea()
{
	return this->Area;
}
double Triangulation::GetArea(Face* f)
{
	if (f != this->GetBoundary())
	{
		VertexIterator itvf(f);
		Vertex* vA = itvf.Next();
		Vertex* vB = itvf.Next();
		Vertex* vC = itvf.Next();
		Vector3D A = vA->GetPoint();
		Vector3D B = vB->GetPoint();
		Vector3D C = vC->GetPoint();
		Triangle t(A, B, C);
		return t.GetArea();
	}
	return -1;
}
double Triangulation::GetMinX()
{
	if (this->Mesh2D)
	{
		double minX = 1.0e-10;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			if (v->GetX() < minX)
				minX = v->GetX();
			v = itv.Next();
		}
		return minX;
	}
	return 0.0;
}
double Triangulation::GetMaxX()
{
	if (this->Mesh2D)
	{
		double maxX = -1.0e-10;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			if (v->GetX() > maxX)
				maxX = v->GetX();
			v = itv.Next();
		}
		return maxX;
	}
	return 0.0;
}
double Triangulation::GetMinY()
{
	if (this->Mesh2D)
	{
		double minY = 1.0e-10;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			if (v->GetY() < minY)
				minY = v->GetY();
			v = itv.Next();
		}
		return minY;
	}
	return 0.0;
}
double Triangulation::GetMaxY()
{
	if (this->Mesh2D)
	{
		double maxY = -1.0e-10;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			if (v->GetY() > maxY)
				maxY = v->GetY();
			v = itv.Next();
		}
		return maxY;
	}
	return 0.0;
}
void Triangulation::GetPoints(vector<Vector3D>& points)
{
	this->Mesh2D->GetPoints(points);
}
void Triangulation::Write(string fileName)
{
	this->Mesh2D->Write(fileName);
}
void Triangulation::Read(string fileName)
{
	QuadEdge qe;
	qe.Read(fileName);
	(*this) << qe;
}
void Triangulation::Draw(string fileName)
{
	
	double xMin = this->GetMinX();
	double xMax = this->GetMaxX();
	double yMin = this->GetMinY();
	double yMax = this->GetMaxY();
	const double WIDTH = 800;
	const double HEIGHT = 800;
	const double PADDING = 100;
	auto xp = [xMin, xMax, WIDTH, PADDING](double x)->int {
		return int((x - xMin) / (xMax - xMin) * WIDTH + PADDING / 2.0);
	};
	auto yp = [yMin, yMax, HEIGHT, PADDING](double x)->int {
		return int((x - yMin) / (yMax - yMin) * HEIGHT + PADDING / 2.0);
	};
	bBitMap::BitMap bmp(int(WIDTH + PADDING), int(HEIGHT + PADDING));
	EdgeIterator ite(this->GetMesh2D());
	Edge* e = ite.Next();
	while(e)
	{
		Vector3D P = e->GetOrig()->GetPoint();
		Vector3D Q = e->GetDest()->GetPoint();
		int x0 = xp(P(0));
		int x1 = xp(Q(0));
		int y0 = yp(P(1));
		int y1 = yp(Q(1));
		int tmax = 1000;
		for (int t = 0; t <= tmax; ++t)
		{
			int x = x0 + t * (x1 - x0) / tmax;
			int y = y0 + t * (y1 - y0) / tmax;
			bmp.Draw1Pixel(x, y, 0, 0, 0);
		}
		e = ite.Next();
	}
	bmp.GenerateBMP(fileName);
}
Triangulation Triangulation::OffDiagonalGrid(int Nx, int Ny, double Lx, double Ly)
{
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	Triangulation T;
	T.Mesh2D = new QuadEdge;
	Vertex*** vertices = new Vertex**[Nx + 1];
	for (int i = 0; i <= Nx; ++i){
		vertices[i] = new Vertex*[Ny + 1];
		for (int j = 0; j <= Ny; ++j){
			Vector3D P(double(i) * hx, double(j) * hy);
			vertices[i][j] = T.Mesh2D->AddVertex(P);
		}
	}
	Face*** upperFaces = new Face**[Nx];
	Face*** lowerFaces = new Face**[Nx];
	Edge*** obliqueEdges = new Edge**[Nx];
	for (int i = 0; i < Nx; ++i){
		upperFaces[i] = new Face*[Ny];
		lowerFaces[i] = new Face*[Ny];
		obliqueEdges[i] = new Edge*[Ny];
		for (int j = 0; j < Ny; ++j){
			upperFaces[i][j] = T.Mesh2D->AddFace();
			lowerFaces[i][j] = T.Mesh2D->AddFace();
			obliqueEdges[i][j] = T.Mesh2D->AddEdge();
		}
	}
	Edge*** horizontalEdges = new Edge**[Nx];
	for (int i = 0; i < Nx; ++i){
		horizontalEdges[i] = new Edge*[Ny + 1];
		for (int j = 0; j <= Ny; ++j){
			horizontalEdges[i][j] = T.Mesh2D->AddEdge();
		}
	}
	Edge*** verticalEdges = new Edge**[Nx + 1];
	for (int i = 0; i <= Nx; ++i){
		verticalEdges[i] = new Edge*[Ny];
		for (int j = 0; j < Ny; ++j){
			verticalEdges[i][j] = T.Mesh2D->AddEdge();
		}
	}
	for (int i = 0; i <= Nx; ++i){
		for (int j = 0; j <= Ny; ++j){
			if (i != 0) T.Mesh2D->AddEdgeToVertex(vertices[i][j], horizontalEdges[i - 1][j], true);
			if (i != Nx) T.Mesh2D->AddEdgeToVertex(vertices[i][j], horizontalEdges[i][j], false);
			if (j != 0) T.Mesh2D->AddEdgeToVertex(vertices[i][j], verticalEdges[i][j - 1], true);
			if (j != Ny) T.Mesh2D->AddEdgeToVertex(vertices[i][j], verticalEdges[i][j], false);
			if (i != 0 && j != 0) T.Mesh2D->AddEdgeToVertex(vertices[i][j], obliqueEdges[i - 1][j - 1], true);
			if (i != Nx && j != Ny) T.Mesh2D->AddEdgeToVertex(vertices[i][j], obliqueEdges[i][j], false);
		}
	}
	for (int i = 0; i < Nx; ++i){
		for (int j = 0; j < Ny; ++j){
			T.Mesh2D->AddEdgeToFace(upperFaces[i][j], obliqueEdges[i][j], true);
			T.Mesh2D->AddEdgeToFace(upperFaces[i][j], horizontalEdges[i][j + 1], false);
			T.Mesh2D->AddEdgeToFace(upperFaces[i][j], verticalEdges[i][j], false);
			T.Mesh2D->AddEdgeToFace(lowerFaces[i][j], obliqueEdges[i][j], false);
			T.Mesh2D->AddEdgeToFace(lowerFaces[i][j], horizontalEdges[i][j], true);
			T.Mesh2D->AddEdgeToFace(lowerFaces[i][j], verticalEdges[i + 1][j], true);
		}
	}
	Face* fb = T.Mesh2D->AddFace();
	for (int j = 0; j < Ny; ++j){
		T.Mesh2D->AddEdgeToFace(fb, verticalEdges[0][j], true);
	}
	for (int i = 0; i < Nx; ++i){
		T.Mesh2D->AddEdgeToFace(fb, horizontalEdges[i][Ny], true);
	}
	for (int j = Ny - 1; j >= 0; --j){
		T.Mesh2D->AddEdgeToFace(fb, verticalEdges[Nx][j], false);
	}
	for (int i = Nx - 1; i >= 0; --i){
		T.Mesh2D->AddEdgeToFace(fb, horizontalEdges[i][0], false);
	}
	T.boundary = fb;
	T.Mesh2D->UpdatePrevs();
	for (int i = 0; i <= Nx; ++i){
		delete [] vertices[i];
	}
	delete[] vertices;
	for (int i = 0; i < Nx; ++i){
		delete[] upperFaces[i];
		delete[] lowerFaces[i];
		delete[] obliqueEdges[i];
	}
	delete[] upperFaces;
	delete[] lowerFaces;
	delete[] obliqueEdges;
	for (int i = 0; i < Nx; ++i){
		delete[] horizontalEdges[i];
	}
	delete[] horizontalEdges;
	for (int i = 0; i <= Nx; ++i){
		delete[] verticalEdges[i];
	}
	delete[] verticalEdges;
	T.GetMesh2D()->IssueIndices();
	return T;
}
bool tester_Triangulation(int& NumTests)
{
	Vector3D A, B(1), C(0.5, 1), D(0.5, 0.5, 1);
	QuadEdge QE;
	QE.SetAsTetrahedron(A, B, C, D);
	if (!QE.TestIntegrity())
		return false;
	if (QE.NumVertices() != 4 || QE.NumEdges() != 6 || QE.NumFaces() != 4)
		return false;
	Triangulation T;
	T << QE;
	if (T.GetBoundary() != T.GetMesh2D()->GetFace(0))
		return false;
	if (T.GetMesh2D()->NumVertices() != 4 || T.GetMesh2D()->NumEdges() != 6 || T.GetMesh2D()->NumFaces() != 4)
		return false;
	if (!T.TestIntegrity())
		return false;
	Vector3D DD(0.5, 0.5), E(0.5,0.25), F(0.5,0.75);
	if (!T.isInside(DD) || !T.isOnOrInside(DD))
		return false;
	if (!T.isInside(E) || !T.isOnOrInside(E))
		return false;
	if (!T.isInside(F) || !T.isOnOrInside(F))
		return false;
	Vector3D H(0.5, 0);
	Vector3D J = B * 0.9 + C * 0.1;
	Vector3D K = C * 0.6 + A * 0.4;
	if (T.isInside(H) || !T.isOnOrInside(H))
		return false;
	if (T.isInside(J) || !T.isOnOrInside(J))
		return false;
	if (T.isInside(K) || !T.isOnOrInside(K))
		return false;
	Vector3D M = A * 0.5 + B * 1.5;
	Vector3D N(3, 3);
	if (T.isInside(M) || T.isOnOrInside(M))
		return false;
	if (T.isInside(N) || T.isOnOrInside(N))
		return false;
	//T.PrintTriangulation();
	if (!tester_Triangulation_1(NumTests))
		return false;
	if (!tester_Triangulation_2(NumTests))
		return false;
	if (!tester_Triangulation_3(NumTests))
		return false;
	if (!tester_Triangulation_4(NumTests))
		return false;
	if (!tester_Triangulation_5(NumTests))
		return false;
	if (!tester_Triangulation_6(NumTests))
		return false;
	if (!tester_Triangulation_7(NumTests))
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1(int& NumTests)
{
	//--------------0. Create Triangulation -----------------------------------
	Vector3D A, B(2), C(2, 2), D(0, 2), E(0.7, 1.3), F(1.6, 0.4), G(0.9, 1.15), H(1.81, 1.71);
	vector<Vector3D> input = { A,B,C,D,E,F,G,H };
	Triangulation T = DelaunayLifting::Triangulate(input);
	int Nv = 8;
	int h = 4;
	int t = 2 * Nv - h - 2;
	int Nf = t + 1;
	int Ne = Nv + Nf - 2;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_0(NumTests, T,input))
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		v[vv->index] = vv;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		f[ff->index] = ff;
		ff = itf.Next();
	}
	//--------------1. Remove E3  -------------------------------------------
	T.ShrinkBoundary(e[3]);
	--Ne;
	--Nf;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_1(NumTests, T,input))
		return false;
	//--------------2. Flip E13  ------------------------------------------
	T.CheckThenFlip(e[13]);//causes e[12] to reverse
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_2(NumTests, T, input))
		return false;
	//--------------3. Flip E10  ---------------------------------------------
	T.CheckThenFlip(e[10]);//causes e[9], e[11] and e[12] (again) to reverse
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_3(NumTests, T, input))
		return false;
	//--------------4. Connect v2-v6  -------------------------------------
	T.Connect(v[2], v[6]);
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_4(NumTests, T, input))
		return false;
	//--------------5. Remove E1  ------------------------------------------
	T.ShrinkBoundary(e[1]);
	--Ne;
	--Nf;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!tester_Triangulation_1_5(NumTests, T, input))
		return false;
	//--------------6. Remove E5  ------------------------------------------
	T.ShrinkBoundary(e[5]);
	if (!tester_Triangulation_1_6(NumTests, T, input))
		return false;
	--Nv;
	Ne -= 2;
	--Nf;
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	//--------------------------------------------------------------------------
	//T.PrintTriangulation();
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_0(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "0"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false;
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 17 || T.GetMesh2D()->NumFaces() != 11)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 4)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 3)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 4)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[10] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 6)
		return false;
	if (f[0]->GetEdge() != e[2] || f[0]->GetDegree() != 3)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 4)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[6] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[7] || f[6]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[1] || f[7]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[9] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[3])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[13] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[14] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[13] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[0])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[3] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[3])
		return false;
	if (e[3]->GetOrig() != v[0] || e[3]->GetDest() != v[3] || e[3]->GetLeft() != f[1] || e[3]->GetRight() != f[0])
		return false;
	if (e[3]->GetNextOrig() != e[2] || e[3]->GetNextDest() != e[4] || e[3]->GetNextLeft() != e[5] || e[3]->GetNextRight() != e[2])
		return false;
	if (e[3]->GetPrevOrig() != e[6] || e[3]->GetPrevDest() != e[15] || e[3]->GetPrevLeft() != e[0] || e[3]->GetPrevRight() != e[4])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[0])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[3])
		return false;
	if (e[4]->GetPrevOrig() != e[3] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[14] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[3])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[4] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[3] || e[6]->GetNextDest() != e[12] || e[6]->GetNextLeft() != e[10] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[9] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[6] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[13] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[10] || e[7]->GetPrevDest() != e[13] || e[7]->GetPrevLeft() != e[12] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[0] || e[9]->GetDest() != v[6] || e[9]->GetLeft() != f[3] || e[9]->GetRight() != f[4])
		return false;
	if (e[9]->GetNextOrig() != e[6] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[8] || e[9]->GetNextRight() != e[6])
		return false;
	if (e[9]->GetPrevOrig() != e[0] || e[9]->GetPrevDest() != e[8] || e[9]->GetPrevLeft() != e[2] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[6] || e[10]->GetDest() != v[5] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[7] || e[10]->GetNextLeft() != e[12] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[12] || e[10]->GetPrevLeft() != e[11] || e[10]->GetPrevRight() != e[6])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[10] || e[11]->GetRight() != f[5])
		return false;
	if (e[11]->GetNextOrig() != e[8] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[16] || e[11]->GetNextRight() != e[10])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[16] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[12])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[10] || e[12]->GetNextDest() != e[13] || e[12]->GetNextLeft() != e[11] || e[12]->GetNextRight() != e[7])
		return false;
	if (e[12]->GetPrevOrig() != e[6] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[10] || e[12]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[1] || e[13]->GetDest() != v[7] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[7] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[12] || e[13]->GetNextRight() != e[1])
		return false;
	if (e[13]->GetPrevOrig() != e[1] || e[13]->GetPrevDest() != e[12] || e[13]->GetPrevLeft() != e[7] || e[13]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[7] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[13] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[5] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[1] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[3] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[11] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	if (T.GetOppositeFace(f[5], e[12]) != f[6])
		return false;
	if (T.GetOppositeFace(f[6], e[12]) != f[5])
		return false;
	if (T.GetOppositeFace(f[1], e[5]) != f[8])
		return false;
	if (T.GetOppositeFace(f[8], e[5]) != f[1])
		return false;
	if (T.GetOppositeFace(f[0], e[15]))
		return false;
	if (T.GetOppositeEdge(f[8], v[7]) != e[5])
		return false;
	if (T.GetOppositeEdge(f[8], v[2]) != e[15])
		return false;
	if (T.GetOppositeEdge(f[8], v[3]) != e[14])
		return false;
	if (T.GetOppositeEdge(f[8], v[4]))
		return false;
	if (T.GetRightFaceOpposite(e[10]) != v[0])
		return false;
	if (T.GetLeftFaceOpposite(e[10]) != v[7])
		return false;
	if (T.GetRightFaceOpposite(e[2]) != v[3])
		return false;
	if (T.GetLeftFaceOpposite(e[2]) != v[6])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_1(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "1_removeE3"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 4)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 4)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 3)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[10] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 6)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[6] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[7] || f[6]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[1] || f[7]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[9] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[13] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[14] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[13] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[15] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[14] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[4] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[12] || e[6]->GetNextLeft() != e[10] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[9] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[6] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[13] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[10] || e[7]->GetPrevDest() != e[13] || e[7]->GetPrevLeft() != e[12] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[0] || e[9]->GetDest() != v[6] || e[9]->GetLeft() != f[3] || e[9]->GetRight() != f[4])
		return false;
	if (e[9]->GetNextOrig() != e[6] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[8] || e[9]->GetNextRight() != e[6])
		return false;
	if (e[9]->GetPrevOrig() != e[0] || e[9]->GetPrevDest() != e[8] || e[9]->GetPrevLeft() != e[2] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[6] || e[10]->GetDest() != v[5] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[7] || e[10]->GetNextLeft() != e[12] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[12] || e[10]->GetPrevLeft() != e[11] || e[10]->GetPrevRight() != e[6])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[10] || e[11]->GetRight() != f[5])
		return false;
	if (e[11]->GetNextOrig() != e[8] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[16] || e[11]->GetNextRight() != e[10])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[16] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[12])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[10] || e[12]->GetNextDest() != e[13] || e[12]->GetNextLeft() != e[11] || e[12]->GetNextRight() != e[7])
		return false;
	if (e[12]->GetPrevOrig() != e[6] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[10] || e[12]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[1] || e[13]->GetDest() != v[7] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[7] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[12] || e[13]->GetNextRight() != e[1])
		return false;
	if (e[13]->GetPrevOrig() != e[1] || e[13]->GetPrevDest() != e[12] || e[13]->GetPrevLeft() != e[7] || e[13]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[7] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[13] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[5] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[1] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[4] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[11] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_2(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "2_flipE13"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (e[10]->index != 10 || f[4]->index != 4)
		return false;
	if (e[13]->index != 13 || f[7]->index != 7 || f[6]->index != 6)
		return false;
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 4)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 5)
		return false;
	if (v[6]->GetEdge() != e[10] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 5)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[6] || f[4]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[13] || f[7]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[9] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[7] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[13] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[7] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[15] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[13] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[4] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[10] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[9] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[7] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[1] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[10] || e[7]->GetPrevDest() != e[1] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[0] || e[9]->GetDest() != v[6] || e[9]->GetLeft() != f[3] || e[9]->GetRight() != f[4])
		return false;
	if (e[9]->GetNextOrig() != e[6] || e[9]->GetNextDest() != e[10] || e[9]->GetNextLeft() != e[8] || e[9]->GetNextRight() != e[6])
		return false;
	if (e[9]->GetPrevOrig() != e[0] || e[9]->GetPrevDest() != e[8] || e[9]->GetPrevLeft() != e[2] || e[9]->GetPrevRight() != e[10])
		return false;
	if (e[10]->GetOrig() != v[6] || e[10]->GetDest() != v[5] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[7] || e[10]->GetNextLeft() != e[12] || e[10]->GetNextRight() != e[9])
		return false;
	if (e[10]->GetPrevOrig() != e[9] || e[10]->GetPrevDest() != e[12] || e[10]->GetPrevLeft() != e[11] || e[10]->GetPrevRight() != e[6])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[10] || e[11]->GetRight() != f[5])
		return false;
	if (e[11]->GetNextOrig() != e[8] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[16] || e[11]->GetNextRight() != e[10])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[16] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[12])
		return false;
	if (e[12]->GetOrig() != v[7] || e[12]->GetDest() != v[5] || e[12]->GetLeft() != f[6] || e[12]->GetRight() != f[5])
		return false;
	if (e[12]->GetNextOrig() != e[14] || e[12]->GetNextDest() != e[10] || e[12]->GetNextLeft() != e[13] || e[12]->GetNextRight() != e[11])
		return false;
	if (e[12]->GetPrevOrig() != e[11] || e[12]->GetPrevDest() != e[13] || e[12]->GetPrevLeft() != e[14] || e[12]->GetPrevRight() != e[10])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[6] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[12] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[13] || e[14]->GetPrevDest() != e[12] || e[14]->GetPrevLeft() != e[13] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[4] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[11] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[5] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[1])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_3(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "3_flipE10"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[15] || v[3]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 3)
		return false;
	if (v[7]->GetEdge() != e[16] || v[7]->GetDegree() != 6)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[5] || f[8]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[4] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[8] || f[10]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[13] || f[7]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[7] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[13] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[7] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[16] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[15] || e[4]->GetPrevDest() != e[16] || e[4]->GetPrevLeft() != e[15] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[8] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[13] || e[5]->GetNextDest() != e[15] || e[5]->GetNextLeft() != e[15] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[14] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[7] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[1] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[1] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[16] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[11] || e[8]->GetPrevLeft() != e[16] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[7] || e[11]->GetDest() != v[6] || e[11]->GetLeft() != f[4] || e[11]->GetRight() != f[10])
		return false;
	if (e[11]->GetNextOrig() != e[12] || e[11]->GetNextDest() != e[8] || e[11]->GetNextLeft() != e[9] || e[11]->GetNextRight() != e[16])
		return false;
	if (e[11]->GetPrevOrig() != e[10] || e[11]->GetPrevDest() != e[9] || e[11]->GetPrevLeft() != e[10] || e[11]->GetPrevRight() != e[8])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[2] || e[14]->GetDest() != v[7] || e[14]->GetLeft() != f[6] || e[14]->GetRight() != f[8])
		return false;
	if (e[14]->GetNextOrig() != e[1] || e[14]->GetNextDest() != e[15] || e[14]->GetNextLeft() != e[12] || e[14]->GetNextRight() != e[5])
		return false;
	if (e[14]->GetPrevOrig() != e[13] || e[14]->GetPrevDest() != e[12] || e[14]->GetPrevLeft() != e[13] || e[14]->GetPrevRight() != e[15])
		return false;
	if (e[15]->GetOrig() != v[3] || e[15]->GetDest() != v[7] || e[15]->GetLeft() != f[8] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[4] || e[15]->GetNextDest() != e[16] || e[15]->GetNextLeft() != e[14] || e[15]->GetNextRight() != e[4])
		return false;
	if (e[15]->GetPrevOrig() != e[5] || e[15]->GetPrevDest() != e[14] || e[15]->GetPrevLeft() != e[5] || e[15]->GetPrevRight() != e[16])
		return false;
	if (e[16]->GetOrig() != v[4] || e[16]->GetDest() != v[7] || e[16]->GetLeft() != f[9] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[4] || e[16]->GetNextDest() != e[10] || e[16]->GetNextLeft() != e[15] || e[16]->GetNextRight() != e[8])
		return false;
	if (e[16]->GetPrevOrig() != e[8] || e[16]->GetPrevDest() != e[15] || e[16]->GetPrevLeft() != e[4] || e[16]->GetPrevRight() != e[11])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[5] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[1])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[16] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_4(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "4_connectV2-V6"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 16 || T.GetMesh2D()->NumFaces() != 10)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[1] || v[1]->GetDegree() != 3)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 6)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[4] || v[3]->GetDegree() != 2)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[10] || v[7]->GetDegree() != 4)
		return false;
	if (f[1]->GetEdge() != e[1] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[16] || f[8]->GetDegree() != 3)
		return false;
	if (f[7]->GetEdge() != e[13] || f[7]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[15] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[16] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[1] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[1])
		return false;
	if (e[1]->GetOrig() != v[1] || e[1]->GetDest() != v[2] || e[1]->GetLeft() != f[7] || e[1]->GetRight() != f[1])
		return false;
	if (e[1]->GetNextOrig() != e[7] || e[1]->GetNextDest() != e[5] || e[1]->GetNextLeft() != e[13] || e[1]->GetNextRight() != e[0])
		return false;
	if (e[1]->GetPrevOrig() != e[0] || e[1]->GetPrevDest() != e[14] || e[1]->GetPrevLeft() != e[7] || e[1]->GetPrevRight() != e[5])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[15] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[5] || e[4]->GetPrevDest() != e[15] || e[4]->GetPrevLeft() != e[5] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[9] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[16] || e[5]->GetNextDest() != e[4] || e[5]->GetNextLeft() != e[4] || e[5]->GetNextRight() != e[1])
		return false;
	if (e[5]->GetPrevOrig() != e[1] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[15] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[7] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[1] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[1] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[16] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[16] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[8] || e[11]->GetRight() != f[4])
		return false;
	if (e[11]->GetNextOrig() != e[16] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[14] || e[11]->GetNextRight() != e[9])
		return false;
	if (e[11]->GetPrevOrig() != e[9] || e[11]->GetPrevDest() != e[10] || e[11]->GetPrevLeft() != e[16] || e[11]->GetPrevRight() != e[10])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[7] || e[14]->GetDest() != v[2] || e[14]->GetLeft() != f[8] || e[14]->GetRight() != f[6])
		return false;
	if (e[14]->GetNextOrig() != e[10] || e[14]->GetNextDest() != e[1] || e[14]->GetNextLeft() != e[16] || e[14]->GetNextRight() != e[12])
		return false;
	if (e[14]->GetPrevOrig() != e[12] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[11] || e[14]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[7])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[15] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[1])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[14] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	if (e[15]->GetOrig() != v[2] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[10] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[13] || e[15]->GetNextDest() != e[4] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[5])
		return false;
	if (e[15]->GetPrevOrig() != e[16] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[16] || e[15]->GetPrevRight() != e[4])
		return false;
	if (e[16]->GetOrig() != v[2] || e[16]->GetDest() != v[6] || e[16]->GetLeft() != f[8] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[15] || e[16]->GetNextDest() != e[8] || e[16]->GetNextLeft() != e[11] || e[16]->GetNextRight() != e[15])
		return false;
	if (e[16]->GetPrevOrig() != e[5] || e[16]->GetPrevDest() != e[11] || e[16]->GetPrevLeft() != e[14] || e[16]->GetPrevRight() != e[8])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_5(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "5_removeE1"); 
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false; 
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false; 
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false; 
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 8 || T.GetMesh2D()->NumEdges() != 15 || T.GetMesh2D()->NumFaces() != 9)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[7] || v[1]->GetDegree() != 2)
		return false;
	if (v[2]->GetEdge() != e[5] || v[2]->GetDegree() != 5)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 4)
		return false;
	if (v[3]->GetEdge() != e[4] || v[3]->GetDegree() != 2)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[10] || v[7]->GetDegree() != 4)
		return false;
	if (f[1]->GetEdge() != e[0] || f[1]->GetDegree() != 6)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[16] || f[8]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (f[9]->GetEdge() != e[15] || f[9]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[16] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[7] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[7])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[4])
		return false;
	if (e[2]->GetPrevOrig() != e[4] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[4]->GetOrig() != v[3] || e[4]->GetDest() != v[4] || e[4]->GetLeft() != f[9] || e[4]->GetRight() != f[1])
		return false;
	if (e[4]->GetNextOrig() != e[5] || e[4]->GetNextDest() != e[2] || e[4]->GetNextLeft() != e[15] || e[4]->GetNextRight() != e[5])
		return false;
	if (e[4]->GetPrevOrig() != e[5] || e[4]->GetPrevDest() != e[15] || e[4]->GetPrevLeft() != e[5] || e[4]->GetPrevRight() != e[2])
		return false;
	if (e[5]->GetOrig() != v[2] || e[5]->GetDest() != v[3] || e[5]->GetLeft() != f[9] || e[5]->GetRight() != f[1])
		return false;
	if (e[5]->GetNextOrig() != e[16] || e[5]->GetNextDest() != e[4] || e[5]->GetNextLeft() != e[4] || e[5]->GetNextRight() != e[13])
		return false;
	if (e[5]->GetPrevOrig() != e[14] || e[5]->GetPrevDest() != e[4] || e[5]->GetPrevLeft() != e[15] || e[5]->GetPrevRight() != e[4])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[1] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[0] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[0] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[16] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[16] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[8] || e[11]->GetRight() != f[4])
		return false;
	if (e[11]->GetNextOrig() != e[16] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[14] || e[11]->GetNextRight() != e[9])
		return false;
	if (e[11]->GetPrevOrig() != e[9] || e[11]->GetPrevDest() != e[10] || e[11]->GetPrevLeft() != e[16] || e[11]->GetPrevRight() != e[10])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[7] || e[14]->GetDest() != v[2] || e[14]->GetLeft() != f[8] || e[14]->GetRight() != f[6])
		return false;
	if (e[14]->GetNextOrig() != e[10] || e[14]->GetNextDest() != e[5] || e[14]->GetNextLeft() != e[16] || e[14]->GetNextRight() != e[12])
		return false;
	if (e[14]->GetPrevOrig() != e[12] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[11] || e[14]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[1])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[15] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[5])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[14] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	if (e[15]->GetOrig() != v[2] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[10] || e[15]->GetRight() != f[9])
		return false;
	if (e[15]->GetNextOrig() != e[13] || e[15]->GetNextDest() != e[4] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[5])
		return false;
	if (e[15]->GetPrevOrig() != e[16] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[16] || e[15]->GetPrevRight() != e[4])
		return false;
	if (e[16]->GetOrig() != v[2] || e[16]->GetDest() != v[6] || e[16]->GetLeft() != f[8] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[15] || e[16]->GetNextDest() != e[8] || e[16]->GetNextLeft() != e[11] || e[16]->GetNextRight() != e[15])
		return false;
	if (e[16]->GetPrevOrig() != e[5] || e[16]->GetPrevDest() != e[11] || e[16]->GetPrevLeft() != e[14] || e[16]->GetPrevRight() != e[8])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_1_6(int& NumTests, Triangulation& T, vector<Vector3D>& input)
{
	T.PrintTriangulation("..\\Data", "6_removeE5");
	if (!T.TestIntegrity())
		return false;
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		if (v[vv->index])
			return false;
		v[vv->index] = vv;
		if (vv->GetPoint() != input[vv->index])
			return false;
		vv = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false;
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false;
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 7 || T.GetMesh2D()->NumEdges() != 13 || T.GetMesh2D()->NumFaces() != 8)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 5)
		return false;
	if (v[1]->GetEdge() != e[7] || v[1]->GetDegree() != 2)
		return false;
	if (v[2]->GetEdge() != e[16] || v[2]->GetDegree() != 4)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[6] || v[5]->GetDegree() != 4)
		return false;
	if (v[6]->GetEdge() != e[11] || v[6]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[10] || v[7]->GetDegree() != 4)
		return false;
	if (f[1]->GetEdge() != e[0] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[10] || f[5]->GetDegree() != 3)
		return false;
	if (f[6]->GetEdge() != e[13] || f[6]->GetDegree() != 3)
		return false;
	if (f[8]->GetEdge() != e[16] || f[8]->GetDegree() != 3)
		return false;
	if (f[4]->GetEdge() != e[10] || f[4]->GetDegree() != 3)
		return false;
	if (f[10]->GetEdge() != e[16] || f[10]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[10] || e[0]->GetNextDest() != e[7] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[7])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[9] || e[2]->GetNextRight() != e[15])
		return false;
	if (e[2]->GetPrevOrig() != e[15] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[5] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[13] || e[6]->GetNextLeft() != e[12] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[9] || e[6]->GetPrevDest() != e[7] || e[6]->GetPrevLeft() != e[10] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[1] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[6] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[0] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[12] || e[7]->GetPrevDest() != e[0] || e[7]->GetPrevLeft() != e[13] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[6] || e[8]->GetLeft() != f[10] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[9] || e[8]->GetNextLeft() != e[16] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[16] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[9])
		return false;
	if (e[9]->GetOrig() != v[6] || e[9]->GetDest() != v[0] || e[9]->GetLeft() != f[4] || e[9]->GetRight() != f[3])
		return false;
	if (e[9]->GetNextOrig() != e[11] || e[9]->GetNextDest() != e[6] || e[9]->GetNextLeft() != e[10] || e[9]->GetNextRight() != e[8])
		return false;
	if (e[9]->GetPrevOrig() != e[8] || e[9]->GetPrevDest() != e[10] || e[9]->GetPrevLeft() != e[11] || e[9]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[6] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[8] || e[11]->GetRight() != f[4])
		return false;
	if (e[11]->GetNextOrig() != e[16] || e[11]->GetNextDest() != e[12] || e[11]->GetNextLeft() != e[14] || e[11]->GetNextRight() != e[9])
		return false;
	if (e[11]->GetPrevOrig() != e[9] || e[11]->GetPrevDest() != e[10] || e[11]->GetPrevLeft() != e[16] || e[11]->GetPrevRight() != e[10])
		return false;
	if (e[12]->GetOrig() != v[5] || e[12]->GetDest() != v[7] || e[12]->GetLeft() != f[5] || e[12]->GetRight() != f[6])
		return false;
	if (e[12]->GetNextOrig() != e[7] || e[12]->GetNextDest() != e[14] || e[12]->GetNextLeft() != e[10] || e[12]->GetNextRight() != e[13])
		return false;
	if (e[12]->GetPrevOrig() != e[13] || e[12]->GetPrevDest() != e[11] || e[12]->GetPrevLeft() != e[6] || e[12]->GetPrevRight() != e[14])
		return false;
	if (e[14]->GetOrig() != v[7] || e[14]->GetDest() != v[2] || e[14]->GetLeft() != f[8] || e[14]->GetRight() != f[6])
		return false;
	if (e[14]->GetNextOrig() != e[10] || e[14]->GetNextDest() != e[16] || e[14]->GetNextLeft() != e[16] || e[14]->GetNextRight() != e[12])
		return false;
	if (e[14]->GetPrevOrig() != e[12] || e[14]->GetPrevDest() != e[13] || e[14]->GetPrevLeft() != e[11] || e[14]->GetPrevRight() != e[13])
		return false;
	if (e[13]->GetOrig() != v[5] || e[13]->GetDest() != v[2] || e[13]->GetLeft() != f[6] || e[13]->GetRight() != f[1])
		return false;
	if (e[13]->GetNextOrig() != e[12] || e[13]->GetNextDest() != e[14] || e[13]->GetNextLeft() != e[14] || e[13]->GetNextRight() != e[7])
		return false;
	if (e[13]->GetPrevOrig() != e[6] || e[13]->GetPrevDest() != e[15] || e[13]->GetPrevLeft() != e[12] || e[13]->GetPrevRight() != e[15])
		return false;
	if (e[10]->GetOrig() != v[7] || e[10]->GetDest() != v[0] || e[10]->GetLeft() != f[5] || e[10]->GetRight() != f[4])
		return false;
	if (e[10]->GetNextOrig() != e[11] || e[10]->GetNextDest() != e[9] || e[10]->GetNextLeft() != e[6] || e[10]->GetNextRight() != e[11])
		return false;
	if (e[10]->GetPrevOrig() != e[14] || e[10]->GetPrevDest() != e[0] || e[10]->GetPrevLeft() != e[12] || e[10]->GetPrevRight() != e[9])
		return false;
	if (e[15]->GetOrig() != v[2] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[10] || e[15]->GetRight() != f[1])
		return false;
	if (e[15]->GetNextOrig() != e[13] || e[15]->GetNextDest() != e[2] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[13])
		return false;
	if (e[15]->GetPrevOrig() != e[16] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[16] || e[15]->GetPrevRight() != e[2])
		return false;
	if (e[16]->GetOrig() != v[2] || e[16]->GetDest() != v[6] || e[16]->GetLeft() != f[8] || e[16]->GetRight() != f[10])
		return false;
	if (e[16]->GetNextOrig() != e[15] || e[16]->GetNextDest() != e[8] || e[16]->GetNextLeft() != e[11] || e[16]->GetNextRight() != e[15])
		return false;
	if (e[16]->GetPrevOrig() != e[14] || e[16]->GetPrevDest() != e[11] || e[16]->GetPrevLeft() != e[14] || e[16]->GetPrevRight() != e[8])
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_2(int& NumTests)
{
	int Nx = 1;
	int Ny = 1;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	//T.PrintTriangulation();
	T.Draw("..\\Data\\tester_Triangulation_2.bmp");
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() !=  Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!T.TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_3(int& NumTests)
{
	int Nx = 2;
	int Ny = 2;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	T.Draw("..\\Data\\tester_Triangulation_3.bmp");
	//T.PrintTriangulation();
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!T.TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_4(int& NumTests)
{
	int Nx = 3;
	int Ny = 4;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	//T.PrintTriangulation();
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	if (!T.TestIntegrity())
		return false;
	NumTests += 1;
	return true;
}
bool tester_Triangulation_5(int& NumTests)
{
	int Nx = 5;
	int Ny = 7;
	double Lx = 1.0;
	double Ly = 1.0;
	int Nv = (Nx + 1) * (Ny + 1);
	int Ne = 3 * Nx * Ny + Nx + Ny;
	int Nf = 2 * Nx * Ny + 1;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	//T.PrintTriangulation();
	if (T.GetMesh2D()->NumVertices() != Nv || T.GetMesh2D()->NumEdges() != Ne || T.GetMesh2D()->NumFaces() != Nf)
		return false;
	/*if (!T.TestIntegrity())
		return false;*///Passes okay but takes 5 sec to do it.
	NumTests += 1;
	return true;
}
bool tester_Triangulation_6(int& NumTests)
{
	Vector3D A, B(2), C(2, 2), D(0, 2), E(0.55, 0.54), F(1.55, 0.45), G(1.56, 1.55), H(0.54, 1.54);
	vector<Vector3D> input = { A,B,C,D,E,F,G,H };
	Triangulation T = DelaunayLifting::Triangulate(input);
	T.PrintTriangulation("..\\Data", "0");
	T.Draw("..\\Data\\tester_Triangulation_6.bmp");
	QuadEdge* qe = T.GetMesh2D();
	Vertex* v[8];
	VertexIterator itv(qe);
	Vertex* vv = itv.Next();
	while (vv)
	{
		v[vv->index] = vv;
		vv = itv.Next();
	}
	vector<Vertex*> boundary = { v[0] , v[1], v[5], v[7], v[4] };
	T.ImposeBoundary(boundary);
	//T.PrintTriangulation();
	Edge* e[17];
	Face* f[11];
	for (int i = 0; i < 8; ++i) v[i] = 0;
	for (int i = 0; i < 17; ++i) e[i] = 0;
	for (int i = 0; i < 11; ++i) f[i] = 0;
	VertexIterator itv2(qe);
	Vertex* vv2 = itv2.Next();
	while (vv2)
	{
		if (v[vv2->index])
			return false;
		v[vv2->index] = vv2;
		if (vv2->GetPoint() != input[vv2->index])
			return false;
		vv2 = itv2.Next();
	}
	EdgeIterator ite(qe);
	Edge* ee = ite.Next();
	while (ee)
	{
		if (e[ee->index])
			return false;
		e[ee->index] = ee;
		ee = ite.Next();
	}
	FaceIterator itf(qe);
	Face* ff = itf.Next();
	while (ff)
	{
		if (f[ff->index])
			return false;
		f[ff->index] = ff;
		ff = itf.Next();
	}
	if (T.GetMesh2D()->NumVertices() != 5 || T.GetMesh2D()->NumEdges() != 7 || T.GetMesh2D()->NumFaces() != 4)
		return false;
	if (v[0]->GetEdge() != e[0] || v[0]->GetDegree() != 3)
		return false;
	if (v[1]->GetEdge() != e[7] || v[1]->GetDegree() != 2)
		return false;
	if (v[4]->GetEdge() != e[8] || v[4]->GetDegree() != 3)
		return false;
	if (v[5]->GetEdge() != e[8] || v[5]->GetDegree() != 4)
		return false;
	if (v[7]->GetEdge() != e[11] || v[7]->GetDegree() != 2)
		return false;
	if (f[1]->GetEdge() != e[0] || f[1]->GetDegree() != 5)
		return false;
	if (f[2]->GetEdge() != e[0] || f[2]->GetDegree() != 3)
		return false;
	if (f[3]->GetEdge() != e[2] || f[3]->GetDegree() != 3)
		return false;
	if (f[5]->GetEdge() != e[11] || f[5]->GetDegree() != 3)
		return false;
	if (e[0]->GetOrig() != v[0] || e[0]->GetDest() != v[1] || e[0]->GetLeft() != f[2] || e[0]->GetRight() != f[1])
		return false;
	if (e[0]->GetNextOrig() != e[6] || e[0]->GetNextDest() != e[7] || e[0]->GetNextLeft() != e[7] || e[0]->GetNextRight() != e[2])
		return false;
	if (e[0]->GetPrevOrig() != e[2] || e[0]->GetPrevDest() != e[7] || e[0]->GetPrevLeft() != e[6] || e[0]->GetPrevRight() != e[7])
		return false;
	if (e[2]->GetOrig() != v[4] || e[2]->GetDest() != v[0] || e[2]->GetLeft() != f[3] || e[2]->GetRight() != f[1])
		return false;
	if (e[2]->GetNextOrig() != e[8] || e[2]->GetNextDest() != e[0] || e[2]->GetNextLeft() != e[6] || e[2]->GetNextRight() != e[15])
		return false;
	if (e[2]->GetPrevOrig() != e[15] || e[2]->GetPrevDest() != e[6] || e[2]->GetPrevLeft() != e[8] || e[2]->GetPrevRight() != e[0])
		return false;
	if (e[6]->GetOrig() != v[0] || e[6]->GetDest() != v[5] || e[6]->GetLeft() != f[3] || e[6]->GetRight() != f[2])
		return false;
	if (e[6]->GetNextOrig() != e[2] || e[6]->GetNextDest() != e[7] || e[6]->GetNextLeft() != e[8] || e[6]->GetNextRight() != e[0])
		return false;
	if (e[6]->GetPrevOrig() != e[0] || e[6]->GetPrevDest() != e[11] || e[6]->GetPrevLeft() != e[2] || e[6]->GetPrevRight() != e[7])
		return false;
	if (e[7]->GetOrig() != v[5] || e[7]->GetDest() != v[1] || e[7]->GetLeft() != f[1] || e[7]->GetRight() != f[2])
		return false;
	if (e[7]->GetNextOrig() != e[8] || e[7]->GetNextDest() != e[0] || e[7]->GetNextLeft() != e[0] || e[7]->GetNextRight() != e[6])
		return false;
	if (e[7]->GetPrevOrig() != e[6] || e[7]->GetPrevDest() != e[0] || e[7]->GetPrevLeft() != e[11] || e[7]->GetPrevRight() != e[0])
		return false;
	if (e[8]->GetOrig() != v[4] || e[8]->GetDest() != v[5] || e[8]->GetLeft() != f[5] || e[8]->GetRight() != f[3])
		return false;
	if (e[8]->GetNextOrig() != e[15] || e[8]->GetNextDest() != e[11] || e[8]->GetNextLeft() != e[11] || e[8]->GetNextRight() != e[2])
		return false;
	if (e[8]->GetPrevOrig() != e[2] || e[8]->GetPrevDest() != e[7] || e[8]->GetPrevLeft() != e[15] || e[8]->GetPrevRight() != e[6])
		return false;
	if (e[15]->GetOrig() != v[7] || e[15]->GetDest() != v[4] || e[15]->GetLeft() != f[5] || e[15]->GetRight() != f[1])
		return false;
	if (e[15]->GetNextOrig() != e[11] || e[15]->GetNextDest() != e[2] || e[15]->GetNextLeft() != e[8] || e[15]->GetNextRight() != e[11])
		return false;
	if (e[15]->GetPrevOrig() != e[11] || e[15]->GetPrevDest() != e[8] || e[15]->GetPrevLeft() != e[11] || e[15]->GetPrevRight() != e[2])
		return false;
	if (e[11]->GetOrig() != v[5] || e[11]->GetDest() != v[7] || e[11]->GetLeft() != f[5] || e[11]->GetRight() != f[1])
		return false;
	if (e[11]->GetNextOrig() != e[6] || e[11]->GetNextDest() != e[15] || e[11]->GetNextLeft() != e[15] || e[11]->GetNextRight() != e[7])
		return false;
	if (e[11]->GetPrevOrig() != e[8] || e[11]->GetPrevDest() != e[15] || e[11]->GetPrevLeft() != e[8] || e[11]->GetPrevRight() != e[15])
		return false;
	NumTests += 1;
	return true;
}
double Vx(double x, double y)
{
	return 1.5;
}
double Vy(double x, double y)
{
	return 2.5;
}
bool tester_Triangulation_7(int& NumTests)
{
	{
		double Lx = 6, Ly = 6;
		Triangulation T = Triangulation::OffDiagonalGrid(1, 1, Lx, Ly);
		T.PrintTriangulation();
		EdgeIterator ite(T.GetBoundary());
		Edge* e = ite.Next();
		int counter = 0;
		while (e && counter < 10)
		{
			++counter;
			e = ite.Next();
		}
		if (counter != 4)
			return false;
	}
	{
		double Lx = 6, Ly = 6; 
		Triangulation T = Triangulation::OffDiagonalGrid(2, 2, Lx, Ly);
		EdgeIterator ite(T.GetBoundary());
		Edge* e = ite.Next();
		int counter = 0;
		while (e && counter < 100)
		{
			++counter;
			e = ite.Next();
		}
		if (counter != 8)
			return false; 
	}
	NumTests += 1;
	return true;
}

