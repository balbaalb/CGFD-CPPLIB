#include "NodeComposite.h"
#include "GeoGraphObject.h"
#include "Face.h"
#include "Vertex.h"
#include "Edge.h"
#include "FaceIterator.h"
#include "VertexIterator.h"
#include "EdgeIterator.h"
#include "QuadEdge.h"
#include "bGrid.h"
#include "SquareMatrix.h"
#include "Vector.h"
#include "LineSegment2D.h"
#include "Triangle.h"
#include "ShapeFunction.h"
void NodeComposite::CopyBody(const NodeComposite& rhs)
{
	//this assumes that this->Reset has been called.
	this->qe = rhs.qe;
	this->fb = rhs.fb;
	this->type = rhs.type;
	this->index = this->index2 = 0;
	//implement code for this->constValue and this->constNormalGrad ...
	this->nodes.resize(rhs.nodes.size(), 0);
	for (int i = 0; i < rhs.nodes.size(); ++i)
	{
		if (rhs.nodes[i])
		{
			GeoGraphObject* geoPtr = rhs.nodes[i]->GetHost();
			this->AddNode(geoPtr);
			*(this->nodes[i]) = *(rhs.nodes[i]);
		}
	}
	this->DependentNodes.resize(rhs.DependentNodes.size(), 0);
	for (int i = 0; i < rhs.DependentNodes.size(); ++i)
	{
		if (rhs.DependentNodes[i])
		{
			GeoGraphObject* geoPtr = rhs.DependentNodes[i]->GetHost();
			this->AddDependentNode(geoPtr);
			*(this->DependentNodes[i]) = *(rhs.DependentNodes[i]);
		}
	}
	this->K = 0;
	if (rhs.K)
		this->K = new SquareMatrix(*rhs.K);
	this->C = 0;
	if (rhs.C)
		this->C = new Vector(*rhs.C);
	this->X = 0;
	if (rhs.X)
		this->X = new Vector(*rhs.X);
	this->EdgeCondition = 0;
	if(rhs.EdgeCondition)
		this->SetEdgeCondition(rhs.qe, rhs.EdgeCondition);
}
void NodeComposite::Reset()
{
	this->qe = 0;
	this->fb = 0;
	this->type = CELLS;
	this->index = 0;
	this->index2 = 0;
	this->NodeHash.clear();
	this->isDependent.clear();
	this->constValues.clear();
	this->constNormalGrad.clear();
	for (int i = 0; i < this->nodes.size(); ++i) {
		if (this->nodes[i])
			delete this->nodes[i];
		this->nodes[i] = 0;
	}
	this->nodes.resize(0);
	for (int i = 0; i < this->DependentNodes.size(); ++i) {
		if (this->DependentNodes[i])
			delete this->DependentNodes[i];
		this->DependentNodes[i] = 0;
	}
	this->DependentNodes.resize(0);
	if (this->K)
		delete this->K;
	this->K = 0;
	if (this->C)
		delete this->C;
	this->C = 0;
}
void NodeComposite::AddNode(GeoGraphObject* gPtr)
{
	this->nodes[index] = new Node(gPtr);
	this->nodes[index]->index = index;
	this->NodeHash[gPtr] = this->nodes[index];
	this->isDependent[gPtr] = false;
	++index;
}
void NodeComposite::AddDependentNode(GeoGraphObject* gPtr)
{
	this->DependentNodes[index2] = new Node(gPtr);
	this->DependentNodes[index2]->index = index2;
	this->NodeHash[gPtr] = this->DependentNodes[index2];
	this->isDependent[gPtr] = true;
	++index2;
}
void NodeComposite::populate_basedOn_CELLS()
{
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		Node* n = this->NodeHash[e]; 
		auto constValueIter = this->constValues.find(e);
		if (constValueIter != this->constValues.end())
		{
			n->value = constValueIter->second;
		}
		else
		{
			Face* fL = e->GetLeft();
			Face* fR = e->GetRight();
			int counter = 0;
			double val = 0;
			if (fL != this->fb)
			{
				Node* nL = this->NodeHash[fL];
				val += nL->value;
				++counter;
			}
			if (fR != this->fb)
			{
				Node* nR = this->NodeHash[fR];
				val += nR->value;
				++counter;
			}
			val /= counter ? double(counter) : 1;
			n->value = val;
		}
		e = ite.Next();
	}
	VertexIterator itv(qe);
	Vertex* v = itv.Next();
	while (v)
	{
		Node* n = this->NodeHash[v]; 
		auto constValueIter = this->constValues.find(v);
		if (constValueIter != this->constValues.end())
		{
			n->value = constValueIter->second;
		}
		else
		{
			int counter = 0;
			double val = 0;
			FaceIterator itf(v);
			Face* f = itf.Next();
			while (f)
			{
				if (f != this->fb)
				{
					Node* nf = this->NodeHash[f];
					val += nf->value;
					++counter;
				}
				f = itf.Next();
			}
			val /= counter ? double(counter) : 1;
			n->value = val;
		}
		v = itv.Next();
	}
}
void NodeComposite::populate_basedOn_EDGES()
{
	FaceIterator itf(this->qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != this->fb)
		{
			Node* n = this->NodeHash[f];
			auto constValueIter = this->constValues.find(f);
			if (constValueIter != this->constValues.end())
			{
				n->value = constValueIter->second;
			}
			else
			{
				int counter = 0;
				double val = 0;
				EdgeIterator ite(f);
				Edge* e = ite.Next();
				while (e)
				{
					bool condition = this->EdgeCondition ? (*this->EdgeCondition)[e] : true;
					if (condition)
					{
						Node* ne = this->NodeHash[e];
						val += ne->value;
						++counter;
					}
					e = ite.Next();
				}
				val /= counter ? double(counter) : 1;
				n->value = val;
			}
		}
		f = itf.Next();
	}
	VertexIterator itv(qe);
	Vertex* v = itv.Next();
	while (v)
	{
		if (this->type == EDGES || (this->type == EDGES_AND_BOUNDARY && !this->IsIndependent(v)))
		{
			Node* n = this->NodeHash[v];
			auto constValueIter = this->constValues.find(v);
			if (constValueIter != this->constValues.end())
			{
				n->value = constValueIter->second;
			}
			else
			{
				int counter = 0;
				double val = 0;
				EdgeIterator ite(v);
				Edge* e = ite.Next();
				while (e)
				{
					bool condition = this->EdgeCondition ? (*this->EdgeCondition)[e] : true;
					if (condition)
					{
						Node* ne = this->NodeHash[e];
						val += ne->value;
						++counter;
					}
					e = ite.Next();
				}
				val /= counter ? double(counter) : 1;
				n->value = val;
			}
		}
		v = itv.Next();
	}
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		bool condition = this->EdgeCondition ? (*this->EdgeCondition)[e] : true;
		if (!condition)
		{
			Node* n = this->NodeHash[e];
			auto constValueIter = this->constValues.find(e);
			if (constValueIter != this->constValues.end())
			{
				n->value = constValueIter->second;
			}
			else
			{
				Vertex* vO = e->GetOrig();
				Vertex* vD = e->GetDest();
				Node* nO = this->NodeHash[vO];
				double val = 0.5 * nO->value;
				Node* nD = this->NodeHash[vD];
				val += 0.5 * nD->value;
				n->value = val;
			}
		}
		e = ite.Next();
	}
}
void NodeComposite::populate_basedOn_VERTICES()
{
	FaceIterator itf(this->qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != this->fb)
		{
			Node* n = this->NodeHash[f]; 
			auto constValueIter = this->constValues.find(f);
			if (constValueIter != this->constValues.end())
			{
				n->value = constValueIter->second;
			}
			else
			{
				int counter = 0;
				double val = 0;
				VertexIterator itv(f);
				Vertex* v = itv.Next();
				while (v)
				{
					Node* nv = this->NodeHash[v];
					val += nv->value;
					++counter;
					v = itv.Next();
				}
				val /= counter ? double(counter) : 1;
				n->value = val;
			}
		}
		f = itf.Next();
	}
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		Node* n = this->NodeHash[e];
		auto constValueIter = this->constValues.find(e);
		if (constValueIter != this->constValues.end())
		{
			n->value = constValueIter->second;
		}
		else
		{
			Vertex* vO = e->GetOrig();
			Vertex* vD = e->GetDest();
			Node* nO = this->NodeHash[vO];
			double val = 0.5 * nO->value;
			Node* nD = this->NodeHash[vD];
			val += 0.5 * nD->value;
			n->value = val;
		}
		e = ite.Next();
	}
}
void NodeComposite::populate_basedOn_CELLS_AND_BOUNDARY()
{
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		Node* n = this->NodeHash[e];
		auto constValueIter = this->constValues.find(e);
		if (constValueIter != this->constValues.end())
		{
			n->value = constValueIter->second;
		}
		else
		{
			Face* fL = e->GetLeft();
			Face* fR = e->GetRight();
			if (fR != fb && fL != fb)
			{
				Node* nL = this->NodeHash[fL];
				double val = 0.5 * nL->value;
				Node* nR = this->NodeHash[fR];
				val += 0.5 * nR->value;
				n->value = val;
			}
		}
		e = ite.Next();
	}
	VertexIterator itv(qe);
	Vertex* v = itv.Next();
	while (v)
	{
		Node* n = this->NodeHash[v];
		auto constValueIter = this->constValues.find(v);
		if (constValueIter != this->constValues.end())
		{
			n->value = constValueIter->second;
		}
		else
		{
			int counter = 0;
			double val = 0;
			bool onBoundary = false;
			FaceIterator itf(v);
			Face* f = itf.Next();
			while (f)
			{
				if (f != this->fb)
				{
					Node* nf = this->NodeHash[f];
					val += nf->value;
					++counter;
				}
				else
				{
					onBoundary = true;
					break;
				}
				f = itf.Next();
			}
			if (onBoundary)
			{
				counter = 0;
				val = 0;
				EdgeIterator itev(v);
				Edge* ev = itev.Next();
				while (ev) {
					if (ev->GetRight() == fb || ev->GetLeft() == fb)
					{
						Node* ne = this->NodeHash[ev];
						val += ne->value;
						++counter;
					}
					ev = itev.Next();
				}
				if (counter == 2)
				{
					EdgeIterator itev2(v);
					Edge* e1 = itev2.Next();
					Edge* e2 = itev2.Next();
					FaceIterator itfv(v);
					Face* f0 = itfv.Next();
					if(f0 == this->fb)
						f0 = itfv.Next();
					Vector3D P0 = f0->GetPoint();
					Vector3D P1 = e1->GetPoint();
					Vector3D P2 = e2->GetPoint();
					Vector3D P = v->GetPoint();
					Triangle t(P0, P1, P2);
					ShapeFunction N(t);
					double T0 = this->GetValue(f0);
					double T1 = this->GetValue(e1);
					double T2 = this->GetValue(e2);
					val = N.GetValue(0, P) * T0 + N.GetValue(1, P) * T1 + N.GetValue(2, P) * T2;
					counter = 1;
				}
			}
			val /= counter ? double(counter) : 1;
			n->value = val;
		}
		v = itv.Next();
	}
}
void NodeComposite::populate_basedOn_EDGES_AND_BOUNDARY()
{
	this->populate_basedOn_EDGES();
}
void NodeComposite::ApplyBC_internal(string caller)
{
	for (auto it = this->constValues.begin(); it != this->constValues.end(); ++it)
	{
		if (this->IsIndependent(it->first));
		{
			if (caller == "Solve")
				this->SetValueInEquations(it->first, it->second);
			else if (caller == "Apply") 
			{
				Node* n = this->NodeHash[it->first];
				n->value = it->second;
			}

		}
	}
	if (this->type == NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY)
	{
		for (auto it = this->constNormalGrad.begin(); it != this->constNormalGrad.end(); ++it)
		{
			//Works only for the grid. The cross-diffusion term is not considered.
			Edge* e = dynamic_cast<Edge*>(it->first);
			if (e && this->IsIndependent(it->first))
			{
				Face* f = e->GetRight() == fb ? e->GetLeft() : e->GetRight();
				Node* ne = this->NodeHash[e];
				int ie = ne->index;
				this->K->SetRowZero(ie);
				Vector3D Pe = e->GetPoint();
				Vector3D Pf = f->GetPoint();
				Vector3D Pef = Pe - Pf;
				double dist = Pef.abs();
				if (caller == "Solve")
				{
					this->AddToK(e, e, 1.0);
					this->AddToK(e, f, -1.0);
					this->AddToC(e, it->second * dist);
				}
				else if (caller == "Apply")
				{
					Node* nf = this->NodeHash[f];
					ne->value = nf->value + dist * it->second;
				}
			}
		}
	}
	if (this->type == NODE_COMPOSITE_TYPE::EDGES)
	{
		for (auto it = this->constNormalGrad.begin(); it != this->constNormalGrad.end(); ++it)
		{
			//Works only for the grid. The cross-diffusion term is not considered.
			Edge* e = dynamic_cast<Edge*>(it->first);
			if (e && this->IsIndependent(it->first))
			{
				Face* f = e->GetRight() == fb ? e->GetLeft() : e->GetRight();
				Edge* e2 = e->GetNext(f)->GetNext(f);
				Node* ne = this->NodeHash[e];
				int ie = ne->index;
				this->K->SetRowZero(ie);
				Vector3D Pe = e->GetPoint();
				Vector3D Pe2 = e2->GetPoint();
				Vector3D Pe_e2 = Pe - Pe2;
				double dist = Pe_e2.abs();
				if (caller == "Solve")
				{
					this->AddToK(e, e, 1.0);
					this->AddToK(e, e2, -1.0);
					this->AddToC(e, it->second * dist);
				}
				else if (caller == "Apply")
				{
					Node* ne2 = this->NodeHash[e2];
					ne->value = ne2->value + dist * it->second;
				}
			}
		}
	}
	if (this->type == NODE_COMPOSITE_TYPE::EDGES_AND_BOUNDARY)
	{
		for (auto it = this->constNormalGrad.begin(); it != this->constNormalGrad.end(); ++it)
		{
			//Works only for the grid. The cross-diffusion term is not considered.
			Vertex* v = dynamic_cast<Vertex*>(it->first);
			if (v && this->IsIndependent(it->first))
			{
				EdgeIterator ite(v);
				Edge* e = ite.Next();
				Edge* e2 = 0;
				while (e)
				{
					if (e->GetRight() != this->fb && e->GetLeft() != this->fb) {
						e2 = e;
						break;
					}
					e = ite.Next();
				}
				if (e2)
				{
					Node* nv = this->NodeHash[v];
					int iv = nv->index;
					this->K->SetRowZero(iv);
					Vector3D Pv = v->GetPoint();
					Vector3D Pe2 = e2->GetPoint();
					Vector3D Pv_e2 = Pv - Pe2;
					double dist = Pv_e2.abs();
					if (caller == "Solve")
					{
						this->AddToK(v, v, 1.0);
						this->AddToK(v, e2, -1.0);
						this->AddToC(v, it->second * dist);
					}
					else if (caller == "Apply")
					{
						Node* ne2 = this->NodeHash[e2];
						nv->value = ne2->value + dist * it->second;
					}
				}
			}
		}
	}
}
NodeComposite::NodeComposite()
{
	this->qe = 0;
	this->fb = 0;
	this->type = CELLS;
	this->index = 0;
	this->index2 = 0;
	this->K = 0;
	this->C = 0;
	this->X = 0;
	this->EdgeCondition = 0;
}
NodeComposite::NodeComposite(const NodeComposite& rhs)
{
	this->CopyBody(rhs);
}
NodeComposite::~NodeComposite()
{
	this->Reset();
}
void NodeComposite::operator=(const NodeComposite& rhs)
{
	this->Reset();
	this->CopyBody(rhs);
}
void NodeComposite::SetEdgeCondition(QuadEdge* QE, EdgeConditionMap* f)
{
	//this->EdgeCondition = f;
	delete this->EdgeCondition;
	this->EdgeCondition = new EdgeConditionMap;
	this->EdgeCondition->clear();
	this->EdgeCondition->rehash(f->size());
	EdgeIterator ite(QE);
	Edge* e = ite.Next();
	while (e)
	{
		(*this->EdgeCondition)[e] = (*f)[e];
		e = ite.Next();
	}
}
void NodeComposite::Initialize(QuadEdge* QE, Face* Fb, NODE_COMPOSITE_TYPE Type)
{
	this->Reset();
	this->qe = QE;
	this->fb = Fb;
	this->type = Type;
	int Nf = qe->NumFaces();
	int Nv = qe->NumVertices();
	int Ne = qe->NumEdges();
	this->NodeHash.rehash(Nf + Nv + Ne);
	this->isDependent.rehash(Nf + Nv + Ne);
	this->constValues.rehash(Nf + Nv + Ne);
	this->constNormalGrad.rehash(Nf + Nv + Ne);
	nodes.resize(Nf + Nv + Ne, 0);
	DependentNodes.resize(Nf + Nv + Ne, 0);
	this->index = 0;
	this->index2 = 0;
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			if (type == CELLS || type == CELLS_AND_BOUNDARY)
				this->AddNode(f);
			else
				this->AddDependentNode(f);
		}
		f = itf.Next();
	}
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		bool isBoundary = (e->GetLeft() == fb) || (e->GetRight() == fb);
		bool condition = this->EdgeCondition ? (*this->EdgeCondition)[e] : true;
		if ((type == EDGES && condition) || (isBoundary && type == CELLS_AND_BOUNDARY)
			|| (type == EDGES_AND_BOUNDARY && condition))
			this->AddNode(e);
		else
			this->AddDependentNode(e);
		e = ite.Next();
	}
	VertexIterator itv(qe);
	Vertex* v = itv.Next();
	while (v)
	{
		if (type == VERTICES || (type == EDGES_AND_BOUNDARY && this->IsAtBoundary(v)))
			this->AddNode(v);
		else
			this->AddDependentNode(v);
		v = itv.Next();
	}
	this->AddDependentNode(fb);
}
void NodeComposite::InitializeEquations()
{
	if (this->K)
		delete this->K;
	this->K = 0;
	if (this->C)
		delete this->C;
	this->C = 0; 
	if (this->X)
		delete this->X;
	this->X = 0;
	this->K = new SquareMatrix(index);
	this->C = new Vector(index);
	this->X = new Vector(index);
}
void NodeComposite::AddToK(GeoGraphObject* row, GeoGraphObject* col, double value)
{
	if (!this->isDependent[row] && !this->isDependent[col])
	{
		int i = this->NodeHash[row]->index;
		int j = this->NodeHash[col]->index;
		(*K)(i, j) += value;
	}
}
void NodeComposite::AddToC(GeoGraphObject* row, double value)
{
	if (!this->isDependent[row])
	{ 
		int i = this->NodeHash[row]->index;
		(*C)(i) += value; 
	}
}
void NodeComposite::SetValueInEquations(GeoGraphObject* row, double value)
{
	if (!this->isDependent[row])
	{
		int i = this->NodeHash[row]->index;
		for (int j = 0; j < index; ++j)
			(*K)(i, j) = 0;
		(*C)(i) = value;
		(*K)(i, i) = 1;
	}
}
void NodeComposite::SetValue(GeoGraphObject* objPtr, double value)
{
	Node* n = objPtr ? this->GetNode(objPtr) : 0;
	if(n)
		n->value = value;
}
void NodeComposite::SetValue(GeoGraphObject* objPtr, function<double(const Vector3D& P)> f)
{
	Node* n = objPtr ? this->GetNode(objPtr) : 0;
	if (n)
	{
		Vector3D P = n->GetPoint();
		n->value = f(P);
	}
}
void NodeComposite::SetConstantValue(GeoGraphObject* objPtr, double value)
{
	this->constValues[objPtr] = value;
	this->SetValue(objPtr, value);
}
void NodeComposite::SetConstantValue(GeoGraphObject* objPtr, function<double(const Vector3D& P)> f)
{
	Vector3D P = objPtr->GetPoint();
	double value = f(P);
	this->SetConstantValue(objPtr, value);
}
void NodeComposite::SetConstantNormalGradient(GeoGraphObject* objPtr, double value)
{
	this->constNormalGrad[objPtr] = value;
}
void NodeComposite::SetConstantNormalGradient(GeoGraphObject* objPtr, function<double(const Vector3D& P)> f)
{
	Vector3D P = objPtr->GetPoint();
	double value = f(P);
	this->SetConstantNormalGradient(objPtr, value);
}
void NodeComposite::SetValueAllNodes(function<double(const Vector3D& P)> f)
{
	VertexIterator itv(this->qe);
	Vertex* v = itv.Next();
	while (v) {
		Node* n = this->GetNode(v);
		if (n) {
			Vector3D P = n->GetPoint();
			double val = f(P);
			n->value = val;
		}
		v = itv.Next();
	}
	this->populate_basedOn_VERTICES();
}
void NodeComposite::Solve()
{
	this->ApplyBC_internal("Solve");
	*(this->X) = this->K->Solve(*C);
}
double NodeComposite::SolveAndUpdate(double underRelaxation)
{
	this->Solve();
	double convergence = 0;
	for (int i = 0; i < index; ++i)
	{
		if (fabs(this->nodes[i]->value - (*this->X)(i)) > convergence)
			convergence = fabs(this->nodes[i]->value - (*this->X)(i));
		this->nodes[i]->value = this->nodes[i]->value * (1 - underRelaxation) + (*this->X)(i) * underRelaxation;
	}
	return convergence;
}
double NodeComposite::SolveAndAdd(double underRelaxation)
{
	this->Solve();
	double convergence = 0;
	for (int i = 0; i < index; ++i)
	{
		if (fabs((*this->X)(i)) > convergence)
			convergence = fabs((*this->X)(i));
		this->nodes[i]->value += (*this->X)(i) * underRelaxation;
	}
	return convergence;
}
void NodeComposite::Populate()
{
	if (this->type == CELLS)
		this->populate_basedOn_CELLS();
	else if (this->type == EDGES)
		this->populate_basedOn_EDGES();
	else if (this->type == VERTICES)
		this->populate_basedOn_VERTICES();
	else if (this->type == CELLS_AND_BOUNDARY)
		this->populate_basedOn_CELLS_AND_BOUNDARY();
	else if (this->type == EDGES_AND_BOUNDARY)
		this->populate_basedOn_EDGES_AND_BOUNDARY();
}
Node* NodeComposite::GetNode(GeoGraphObject* objPtr)
{
	return this->NodeHash[objPtr];
}
double NodeComposite::GetValue(GeoGraphObject* objPtr)
{
	return this->NodeHash[objPtr]->value;
}
double NodeComposite::GetValue(const Vector3D& P)//O(n)
{
	VertexIterator itv(this->qe);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D Pv = v->GetPoint();
		if (Pv == P)
		{
			return this->GetValue(v);
		}
		v = itv.Next();
	}
	EdgeIterator ite(this->qe);
	Edge* e = ite.Next();
	while (e)
	{
		Vertex* Vo = e->GetOrig();
		Vertex* Vd = e->GetDest();
		Vector3D Po = Vo->GetPoint();
		Vector3D Pd = Vd->GetPoint();
		LineSegment2D LS(Po, Pd);
		if (LS.isOn(P))
		{
			Vector3D OD = Pd - Po;
			Vector3D OP = P - Po;
			if (OP.abs() < 1.0e-10) {
				return this->GetValue(Vo);
			}
			double alpha = OP.abs() / OD.abs();
			return ((1.0 - alpha) * this->GetValue(Vo) + alpha * this->GetValue(Vd));
		}
		e = ite.Next();
	}
	FaceIterator itf(this->qe);
	Face* f = itf.Next();
	while (f)
	{
		//if f is a quadrilateral, this works only if f is convex.
		if (f != this->fb)
		{
			VertexIterator itvf(f);
			Vertex* v0 = itvf.Next();
			Vertex* v1 = itvf.Next();
			Vertex* v2 = itvf.Next();
			Vertex* v3 = itvf.Next();
			Vector3D p0 = v0->GetPoint();
			Vector3D p1 = v1->GetPoint();
			Vector3D p2 = v2->GetPoint();
			Triangle t(p0, p1, p2);
			if (t.isOnOrInside(P))
			{
				ShapeFunction N(t);
				double T0 = this->GetValue(v0);
				double T1 = this->GetValue(v1);
				double T2 = this->GetValue(v2);
				double N0 = N.GetValue(0, P);
				double N1 = N.GetValue(1, P);
				double N2 = N.GetValue(2, P);
				return (T0 * N0 + T1 * N1 + T2 * N2);
			}
			if (v3)
			{
				Vector3D p3 = v3->GetPoint();
				Triangle t2(p0, p2, p3);
				if (t2.isOnOrInside(P))
				{
					ShapeFunction N(t2);
					double T0 = this->GetValue(v0);
					double T2 = this->GetValue(v2);
					double T3 = this->GetValue(v3);
					double N0 = N.GetValue(0, P);
					double N2 = N.GetValue(1, P);
					double N3 = N.GetValue(2, P);
					return (T0 * N0 + T2 * N2 + T3 * N3);
				}
			}
		}
		f = itf.Next();
	}
	return -1.0;
}
int NodeComposite::GetNodeNumber(GeoGraphObject* objPtr)
{
	return this->NodeHash[objPtr]->index;
}
void NodeComposite::GetResults_Vertices(vector<Node*>& results)
{
	results.resize(0);
	VertexIterator itv(this->qe);
	Vertex* v = itv.Next();
	while (v)
	{
		Node* n = this->GetNode(v);
		results.push_back(0);
		results.back() = new Node(*n);
		v = itv.Next();
	}
}
void NodeComposite::GetResults_Edges(vector<Node*>& results)
{
	results.resize(0);
	EdgeIterator ite(this->qe);
	Edge* e = ite.Next();
	while (e)
	{
		Node* n = this->GetNode(e);
		results.push_back(0);
		results.back() = new Node(*n);
		e = ite.Next();
	}
}
void NodeComposite::GetResults_Cells(vector<Node*>& results)
{
	results.resize(0);
	FaceIterator itf(this->qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != this->fb)
		{
			Node* n = this->GetNode(f);
			results.push_back(0);
			results.back() = new Node(*n);
		}
		f = itf.Next();
	}
}
bool NodeComposite::IsIndependent(GeoGraphObject* objPtr)
{
	return !this->isDependent[objPtr];
}
bool NodeComposite::IsIndependent(Node* n)
{
	return this->IsIndependent(n->GetHost());
}
double NodeComposite::GetSolution(GeoGraphObject* objPtr)
{
	Node* n = this->GetNode(objPtr);
	if (this->IsIndependent(n))
	{
		int i = n->index;
		if (this->X)
			return (*this->X)(i);
	}
	return 0.0;
}
void NodeComposite::GetK(SquareMatrix* Kr)
{
	if (Kr)
		delete Kr;
	Kr = new SquareMatrix(*this->K);
}
void NodeComposite::GetC(Vector* Cr)
{
	if (Cr)
		delete Cr;
	Cr = new Vector(*this->C);
}
void NodeComposite::GetX(Vector* Xr)
{
	if (Xr)
		delete Xr;
	Xr = new Vector(*this->X);
}
void NodeComposite::PrintEquations(ofstream& f, string label)
{
	f << endl << "K" << label << " = ";
	if (this->K)
		this->K->print(f);
	else
		f << "NULL";
	f << endl << "C" << label << " = ";
	if (this->C)
		this->C->print(f);
	else
		f << "NULL";
	f << endl << "X" << label << " = ";
	if (this->X)
		this->X->print(f);
	else
		f << "NULL";
}
bool NodeComposite::IsAtBoundary(Vertex* v)
{
	FaceIterator itf(v);
	Face* f = itf.Next();
	while (f) {
		if (f == this->fb)
			return true;
		f = itf.Next();
	}
	return false;
}
void NodeComposite::ApplyBC()
{
	this->ApplyBC_internal("Apply");
}
bool tester_NodeComposite(int& NumTests) 
{
	double A = 1.65, B = 8.75, T0 = 10.0;//T = Ax + By + T0
	bGrid grid(3, 2, 3.0, 2.0);
	QuadEdge* qe = grid.GetMesh2D();
	Face* fb = grid.GetBoundary();
	NodeComposite TempNodes;
	TempNodes.Initialize(qe, fb, CELLS_AND_BOUNDARY);
	//---- Setting Boundary conditions ---------------------
	EdgeIterator iteb(fb);
	Edge* eb = iteb.Next();
	while (eb)
	{
		double x = eb->GetPoint()(0);
		double y = eb->GetPoint()(1);
		double T = T0 + A * x + B * y;
		TempNodes.SetConstantValue(eb, T);
		eb = iteb.Next();
	}
	//---- Write Equations, boundary conditions will be applied automatically from constant values. 
	TempNodes.InitializeEquations();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb) {
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			double a = 0;
			while (e) 
			{
				Face* f1 = e->GetOtherFace(f);
				if (f1 != fb)
				{
					TempNodes.AddToK(f, f1, -5.0);
					a += 5.0;
				}
				else 
				{
					TempNodes.AddToK(f, e, -10.0);
					a += 10;
				}
				e = ite.Next();
			}
			TempNodes.AddToK(f, f, a);
		}
		f = itf.Next();
	}
	//------- Solve ---------------------------------------------
	TempNodes.SolveAndUpdate();
	TempNodes.Populate();
	//------ Post Processing ------------------------------------
	Vertex* v7 = qe->GetVertex(7);
	Node* nv7 = TempNodes.GetNode(v7);
	if (!nv7)
		return false;
	double Tv7 = nv7->value;
	Vector3D Pv7 = v7->GetPoint();
	Vector3D Pnv7 = nv7->GetPoint();
	if (Pv7 != Pnv7)
		return false;
	double Tv7_correct = A * Pv7(0) + B * Pv7(1) + T0;
	double error_v7 = fabs(Tv7 - Tv7_correct) / Tv7_correct * 100;
	if (error_v7 > 0.01)
		return false;
	double Tv7a = TempNodes.GetValue(Pnv7);
	if (fabs(A * Pnv7(0) + B * Pnv7(1) + T0 - Tv7a) > 1.0e-10)
		return false;
	Edge* e3 = qe->GetEdge(3);
	Vector3D PA = e3->GetOrig()->GetPoint();
	Vector3D PB = e3->GetDest()->GetPoint();
	Vector3D C = (PA * 2.0 / 3.0) + (PB * 1.0 / 3.0);
	double TA = TempNodes.GetValue(e3->GetOrig());
	double TB = TempNodes.GetValue(e3->GetDest());
	double TC = TempNodes.GetValue(C);
	if (fabs(A * C(0) + B * C(1) + T0 - TC) > 1.0e-10)
		return false;
	Vector3D D(1.2, 0.8);
	double TD = TempNodes.GetValue(D);
	if (fabs(A * D(0) + B * D(1) + T0 - TD) > 1.0e-10)
		return false;
	Vector3D E(2.5, 1.25);
	double TE = TempNodes.GetValue(E);
	if (fabs(A * E(0) + B * E(1) + T0 - TE) > 1.0e-10)
		return false;
	Vector3D F(2.5, 1.75);
	double TF = TempNodes.GetValue(F);
	if (fabs(A * F(0) + B * F(1) + T0 - TF) > 1.0e-10)//TF becomes 28.785 because of bug in calculating corner vertices
		return false;
	++NumTests;
	return true;
}

