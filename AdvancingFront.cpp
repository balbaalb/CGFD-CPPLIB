#include <map>
#include "AdvancingFront.h"
#include "bPolygon.h"
#include "LineSegment2D.h"
#include "DelaunayLifting.h"
#include "QuadEdge.h"
#include "bData.h"
#include "MathUtils.h"
int AdvancingFront::counter = 0;
ofstream AdvancingFront::file;
bool AdvancingFront::printSteps = false;
bool AdvancingFront::AdvanceOneStep = false;
FrontNode::FrontNode(const Vector3D& P)
{
	this->point = new Vector3D(P);
	this->next = this->prev = 0;
}
FrontNode::FrontNode(FrontNode* rhs)
{
	this->point = new Vector3D(*(rhs->point));
	this->next = this->prev = 0;
}
FrontNode::~FrontNode()
{
	//cut the chain first:
	if (this->prev)
		this->prev->next = 0;
	if (this->next)
		delete this->next;
	this->prev = 0;
	this->next = 0;
	if (this->point)
		delete this->point;
	this->point = 0;
}
void AdvancingFront::AddToTail(const Vector3D& P)
{
	if (!this->head){
		this->head = new FrontNode(P);
		this->head->next = this->head->prev = this->head;
	}
	else{
		FrontNode* prv = this->head->prev;
		this->head->prev = new FrontNode(P);
		this->head->prev->next = this->head;
		this->head->prev->prev = prv;
		prv->next = this->head->prev;
	}
	++this->size;
}
void AdvancingFront::AddAfterHead(const Vector3D& P)
{
	if (!this->head)
	{
		this->head = new FrontNode(P);
		this->head->next = this->head->prev = this->head;
	}
	else
	{
		FrontNode* nxt = this->head->next;
		this->head->next = new FrontNode(P);
		this->head->next->next = nxt;
		this->head->next->prev = this->head;
		nxt->prev = this->head->next;
	}
	++this->size;
}
AdvancingFront::AdvancingFront(const vector<Vector3D>& input)
{
	this->BoundaryPoints.resize(input.size());
	this->head = 0;
	this->size = 0;
	for (int i = 0; i < input.size(); ++i)
	{
		this->BoundaryPoints[i] = input[i];
		this->AddToTail(input[i]);
	}
	this->LengthScale = -1;
	FrontNode* n = this->head;
	do {
		Vector3D A = *(n->point);
		Vector3D B = *(n->next->point);
		Vector3D AB = B - A;
		double L = AB.abs();
		if (this->LengthScale < 0)
			this->LengthScale = L;
		else if (L < this->LengthScale)
			this->LengthScale = L;
		n = n->next;
	} while (n != this->head);
	this->child1 = 0;
	this->child2 = 0;
	this->Completed = false;
}
AdvancingFront::AdvancingFront(FrontNode* h, double LengthScale_) : LengthScale(LengthScale_)
{
	this->head = 0;
	this->size = 0;
	FrontNode* n = h;
	do{
		if (!n)
			return;
		Vector3D P = *(n->point);
		this->AddToTail(P);
		n = n->next;
	} while (n != h);
	this->child1 = 0;
	this->child2 = 0;
	this->Completed = false;
}
AdvancingFront::~AdvancingFront()
{
	if (this->head)
		delete this->head;
	this->head = 0;
}
Vector3D AdvancingFront::SuggestPoint() const
{
	//IMPROVE: choose the minimum length edge from front
	Vector3D A = *(this->head->point);
	Vector3D B = *(this->head->next->point);
	Vector3D M = (A + B) / 2.0;
	Vector3D AB = B - A;
	Vector3D t = AB / AB.abs();
	Vector3D k(0, 0, 1);
	Vector3D n = k && t;
	double h = AB.abs() * (sqrt(3.0) / 2.0);
	Vector3D C = M + n * h;
	C(2) = 0;//to avoid float precision problems
	return C;
}
bool AdvancingFront::IsNewPointAcceptable(const Vector3D& C) const
{
	Vector3D A = *(this->head->point);
	Vector3D B = *(this->head->next->point);
	Vector3D AB = B - A;
	if (!this->IsInsideFront(C))//O(n)
		return false;
	if (this->IsIntersectingFront(C))//O(n)
		return false;
	if (IsTooCloseToFront(C, this->LengthScale * 0.9))//O(n)//0.9 IS ARBITARY, NEEDS IMPROVEMENT
		return false;
	return true;
}
void AdvancingFront::AdmitNewPointOverHead(const Vector3D& C)
{
	this->AddAfterHead(C);
	this->InsertedNodes.push_back(C);
	this->head = this->head->next->next;
}
bool AdvancingFront::IsInsideFront(const Vector3D& C) const//O(n^3)
{
	shared_ptr<bPolygon> FrontPolygon = this->CreatePolygonFromFront();
	bool isInside = FrontPolygon->isInside(C);//O(n^3)
	return isInside;
}
bool AdvancingFront::IsOnOrInsideFront(const Vector3D& C) const
{
	shared_ptr<bPolygon> FrontPolygon = this->CreatePolygonFromFront();
	bool isOnOrInside = FrontPolygon->isOnOrInside(C);
	return isOnOrInside;
}
bool AdvancingFront::IsIntersectingFront(const Vector3D& C, bool existingNode) const//O(n)
{
	//Are the lines AC and BC cut any other exisiting front segments? Here A is head and B is head->next.
	FrontNode* a = this->head;
	FrontNode* b = this->head->next;
	Vector3D A = *(a->point);
	Vector3D B = *(b->point);
	Vector3D AP = *(a->prev->point);
	Vector3D BN = *(b->next->point);
	LineSegment2D AC(A, C);
	LineSegment2D BC(B, C);
	FrontNode* n = this->head->next;
	do{
		FrontNode* m = n->next;
		Vector3D N = *(n->point);
		Vector3D M = *(m->point);
		LineSegment2D NM(N, M);
		if ((!existingNode || C != AP) && AC.isIntersecting(NM))
			return true;
		if ((!existingNode || C != BN) && BC.isIntersecting(NM))
			return true;
		n = m;
	} while (n != this->head);
	return false;
}
bool AdvancingFront::IsTooCloseToFront(const Vector3D& C, double LengthScale_) const//O(n)
{
	FrontNode* n = this->head->next;
	while (n != this->head){
		Vector3D P = *(n->point);
		Vector3D Q = *(n->next->point);
		Vector3D PQ = Q - P;
		if (P.distance(C) < LengthScale_){
			return true;
		}
		if (Q.distance(C) < LengthScale_){
			return true;
		}
		n = n->next;
	}
	return false;
}
shared_ptr<bPolygon> AdvancingFront::CreatePolygonFromFront() const//O(n)
{
	shared_ptr<bPolygon> FrontPolygon = make_shared<bPolygon>();
	FrontNode* n = this->head;
	do {
		FrontPolygon->AddVertex(*(n->point));
		n = n->next;
	} while (n != this->head);
	FrontPolygon->Close();
	return FrontPolygon;
}
shared_ptr<bPolygon> AdvancingFront::CreateLastChildPolygonFromFront() const
{
	shared_ptr<bPolygon> FrontPolygon;
	if (this->child1)
	{
		return this->child1->CreateLastChildPolygonFromFront();
	}
	if (this->child2)
	{
		return this->child2->CreateLastChildPolygonFromFront();
	}
	return this->CreatePolygonFromFront();
}
AdvancingFront* AdvancingFront::LastChild()
{
	shared_ptr<bPolygon> FrontPolygon;
	if (this->child1)
	{
		return this->child1->LastChild();
	}
	if (this->child2)
	{
		return this->child2->LastChild();
	}
	return this;
}
bool AdvancingFront::IsCompleted() const
{
	return this->Completed;
}
Vector3D AdvancingFront::GetHeadPoint() const
{
	return *(this->head->point);
}
Vector3D AdvancingFront::GetHeadNextPoint() const
{
	return *(this->head->next->point);
}
FrontNode* AdvancingFront::NodeClosestTo(const Vector3D& C)//O(n) Omega(n^2)
{
	FrontNode* n = this->head->next->next;
	Vector3D A = *(this->head->point);
	Vector3D B = *(this->head->next->point);
	Vector3D AB = B - A;
	multimap<double, FrontNode*> distances;
	do{
		Vector3D N = *(n->point);
		bool acceptable = true;
		if (acceptable && n == this->head->next->next)
		{
			Vector3D BN = N - B;
			acceptable = (AB.CrossProduct2D(BN) > 1.0e-10);
		}
		if (acceptable && n == this->head->prev)
		{
			Vector3D NA = A - N;
			acceptable = (NA.CrossProduct2D(AB) > 1.0e-10);
		}
		if (acceptable)
		{
			double L = n->point->distance(C);
			distances.insert(pair<double, FrontNode*>(L, n));
		}
		n = n->next;
	} while (n != this->head);
	for (auto it = distances.begin(); it != distances.end(); ++it)
	{
		Vector3D N = *(it->second->point);
		bool acceptable = !this->IsIntersectingFront(N, true);//O(n)
		if (acceptable && it->second != this->head->prev && it->second != this->head->next->next)
		{
			Vector3D ANM = (A + N) / 2.0;
			Vector3D BNM = (B + N) / 2.0;
			acceptable = this->IsInsideFront(ANM) && this->IsInsideFront(BNM);//O(n)
		}
		if (acceptable)
		{
			return it->second;
		}
	}
	return 0;
}
void AdvancingFront::RemoveNode(FrontNode* n)
{
	FrontNode* nxt = n->next;
	FrontNode* prv = n->prev;
	nxt->prev = prv;
	prv->next = nxt;
	if (n == this->head)
		this->head = nxt;
	n->next = n->prev = 0;
	delete n;
}
void AdvancingFront::CombineInsertedNodes(const AdvancingFront& Front)//O(n)
{
	for (int i = 0; i < Front.InsertedNodes.size(); ++i){
		this->InsertedNodes.push_back(Front.InsertedNodes[i]);
	}
}
void AdvancingFront::DivideFront(FrontNode* m)
{
	FrontNode* m2 = new FrontNode(m);
	m2->prev = m->prev;
	m2->prev->next = m2;
	m2->next = this->head->next;
	FrontNode* h2 = this->head->next;
	h2->prev = m2;
	this->head->next = m;
	m->prev = this->head;
	this->child1 = make_unique<AdvancingFront>(this->head,this->LengthScale);
	this->child2 = make_unique<AdvancingFront>(h2,this->LengthScale);
}
void AdvancingFront::AdvanceChildren1Step()
{
	if (this->child1)
	{
		this->child1->Advance();
		if (this->child1->Completed)
		{
			this->CombineInsertedNodes(*this->child1);//O(n)
			this->child1.reset();
		}
		return;
	}
	if (!this->child1 && this->child2)
	{
		this->child2->Advance();
		if (this->child2->Completed)
		{
			this->CombineInsertedNodes(*this->child2);//O(n)
			this->child2.reset();
			this->Completed = true;
		}
		return;
	}
	return;
}
void AdvancingFront::MakeComplete()
{
	if (!this->Completed)
	{
		delete this->head;
		this->head = 0;
		this->Completed = true;
	}
}
bool AdvancingFront::FrontIntegrity() const
{
	shared_ptr<bPolygon> FrontPolygon = this->CreatePolygonFromFront();
	if (!FrontPolygon->TestIntegrity())
		return false;
	return true;
}
void AdvancingFront::UpdateSize()
{
	int frontLength = 0;
	FrontNode* n = this->head;
	while (n)
	{
		++frontLength;
		n = n->next;
		if (n == this->head)
			n = 0;
	}
	this->size = frontLength;
}
double AdvancingFront::GetMaxFrontLength() const
{
	FrontNode* n = this->head;//For Debugging only
	double Lmax = 0;
	while (n)
	{
		FrontNode* m = n->next;
		Vector3D N = *(n->point);
		Vector3D M = *(m->point);
		Vector3D NM = M - N;
		double L = NM.abs();
		Lmax = fmax(L,Lmax);
		n = n->next;
		if (n = this->head)
			n = 0;
	}
	return Lmax;
}
void AdvancingFront::report()
{
	if (AdvancingFront::printSteps)
	{
		this->UpdateSize();
		double maxFrontSide = this->GetMaxFrontLength();
		AdvancingFront::file << endl << AdvancingFront::counter << "," << this << "," << this->size << "," << maxFrontSide;
		FrontNode* n = this->head;
		int i = 0;
		while (n)
		{
			Vector3D* P = n->point;
			AdvancingFront::file << ",P#" << i << "," << (*P)(0) << "," << (*P)(1);
			++i;
			n = n->next;
			if (n == this->head)
				n = 0;
		}
	}
	++AdvancingFront::counter;
}
void AdvancingFront::Advance()//O(n^3)
{
	this->report();//reports the status of the caller
	if (AdvancingFront::AdvanceOneStep && (this->child1 || this->child2))
	{
		this->AdvanceChildren1Step();
		return;
	}
	if (!this->head || this->head->next->next->next == this->head)
	{
		this->MakeComplete();
		return;
	}
	Vector3D C = this->SuggestPoint();
	if (this->IsNewPointAcceptable(C))
	{
		this->AdmitNewPointOverHead(C);
		if (AdvancingFront::AdvanceOneStep)
			return;
		this->Advance();
		return;
	}
	FrontNode* m = this->NodeClosestTo(C);//O(n^2) ChokePoint //could be unified with IsTooCloseToFront(...)
	if (!m || m == this->head || m == this->head->next)
		return;
	if (m == this->head->next->next){
		this->RemoveNode(this->head->next);
		this->head = m;
		if (AdvancingFront::AdvanceOneStep)
			return;
		this->Advance();
		return;
	}
	if (m == this->head->prev){
		this->RemoveNode(this->head);
		if (AdvancingFront::AdvanceOneStep)
			return;
		this->Advance();
		return;
	}
	this->DivideFront(m);
	if (AdvancingFront::AdvanceOneStep)
		return;
	this->child1->Advance();
	this->child2->Advance();
	this->CombineInsertedNodes(*this->child1);//O(n)
	this->CombineInsertedNodes(*this->child2);//O(n)
	return;
}
void AdvancingFront::GetPoints(vector<Vector3D>& points)
{
	points.resize(this->BoundaryPoints.size() + this->InsertedNodes.size());
	for (int i = 0; i < this->BoundaryPoints.size(); ++i)
	{
		points[i] = this->BoundaryPoints[i];
	}
	for (int i = this->BoundaryPoints.size(); i < this->BoundaryPoints.size() + this->InsertedNodes.size(); ++i)
	{
		points[i] = this->InsertedNodes[i - this->BoundaryPoints.size()];
	}
}
Triangulation AdvancingFront::Tessellate(const vector<Vector3D>& BoundaryPoints)//static
{
	AdvancingFront::AdvanceOneStep = false; 
	AdvancingFront::counter = 0;
	if (AdvancingFront::printSteps)
	{
		string path = bData::FolderPath + "\\AdvancingFrontReport.txt";
		AdvancingFront::file.open(path.c_str(),fstream::out);
		AdvancingFront::file <<"counter,this,front size, max front side length";
	}
	AdvancingFront AFT(BoundaryPoints);
	AFT.Advance();
	AFT.report();//last status
	vector<Vector3D> points;
	AFT.GetPoints(points);
	vector<int> boundary;
	boundary.resize(BoundaryPoints.size());
	for (int i = 0; i < BoundaryPoints.size(); ++i)
	{
		boundary[i] = i;
	}
	Triangulation T = DelaunayLifting::Triangulate(points,boundary);
	if (AdvancingFront::printSteps)
	{
		AdvancingFront::file.close();
	}
	return T;
}
bool tester_AdvancingFront(int& NumTests)
{
	if (!tester_AdvancingFront_1(NumTests))
		return false;
	if (!tester_AdvancingFront_2(NumTests))
		return false;
	if (!tester_AdvancingFront_3(NumTests))
		return false;
	if (!tester_AdvancingFront_4(NumTests))
		return false;
	if (!tester_AdvancingFront_5(NumTests))
		return false;
	if (!tester_AdvancingFront_6(NumTests))
		return false;
	if (!tester_AdvancingFront_7(NumTests))
		return false;
	if (!tester_AdvancingFront_9(NumTests))
		return false;
	if (!tester_AdvancingFront_10(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront_1(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 1;
	int Ny = 1;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront_2(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 2;
	int Ny = 2;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;/**/
	++NumTests;
	return true;
}
bool tester_AdvancingFront_3(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 3;
	int Ny = 3;
	double Lx = 3;
	double Ly = 3;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;/**/
	++NumTests;
	return true;
}
bool tester_AdvancingFront_4(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 4;
	int Ny = 4;
	double Lx = 4;
	double Ly = 4;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront_5(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 5;
	int Ny = 5;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	Triangulation T = AdvancingFront::Tessellate(input); 
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;/**/
	++NumTests;
	return true;
}
bool tester_AdvancingFront_6(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 6;
	int Ny = 6;
	double Lx = 6;
	double Ly = 6;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	Triangulation T = AdvancingFront::Tessellate(input); 
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;/**/
	++NumTests;
	return true;
}
bool tester_AdvancingFront_7(int& NumTests)
{
	//Like test 5 but  input[2:3] are shifted upward
	vector<Vector3D> input;
	int Nx = 5;
	int Ny = 5;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i) {
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i) {
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i) {
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i) {
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	input.erase(input.begin() + 19);
	input.erase(input.begin() + 18);
	input.erase(input.begin() + 17);
	input.erase(input.begin() + 16);
	input[13](1) *= 2.0;
	input.erase(input.begin() + 11);
	input.erase(input.begin() + 9);
	input.erase(input.begin() + 8);
	input.erase(input.begin() + 7);
	input.erase(input.begin() + 6);
	input[2](1) = 0.5;
	input[3](1) = 0.5;
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	ofstream file;
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront_9(int& NumTests)
{
	//Like test 5 but  input[2:3] are shifted upward
	vector<Vector3D> input;
	double R1 = 1, R2 = 2;
	int Nr = 5;
	int Nt = 5;
	double dr = (R2 - R1) / double(Nr);
	double dt = (pi / 2.0) / double(Nt);
	int n = -1;
	Vector3D buffer;
	for (int i = 0; i < Nr; ++i) {
		buffer(0) = R1 + double(i) * dr; buffer(1) = 0;
		input.push_back(buffer);
	}
	for (int i = 0; i < 2 * Nt; ++i) {
		double theta = double(i) * dt / 2.0;
		buffer(0) = R2 * cos(theta); buffer(1) = R2 * sin(theta);
		input.push_back(buffer);
	}
	for (int i = 0; i < Nr; ++i) {
		buffer(0) = 0; buffer(1) = R2 - double(i) * dr;
		input.push_back(buffer);
	}
	for (int i = 0; i < Nt; ++i) {
		double theta = pi / 2.0 - double(i) * dt;
		buffer(0) = R1 * cos(theta); buffer(1) = R1 * sin(theta);
		input.push_back(buffer);
	}
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	ofstream file;
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront_10(int& NumTests)
{
	int N = 4;
	vector<Vector3D> input;
	Vector3D A, B(1, 0), C(2, 1);
	Vector3D P = A;
	Vector3D dP = B - A;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	P = B;
	dP = C - B;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	P = C;
	dP = A - C;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	Triangulation Triang = AdvancingFront::Tessellate(input);
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	if (fabs(alpha_max - 135) > 0.01)
		return false;
	int Nv = Triang.NumVertices();
	if (Nv != 13)
		return false;
	++NumTests;
	return true;
}