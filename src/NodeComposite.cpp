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
	this->initilized = rhs.initilized;
	this->NeighborsIndex.resize(rhs.NeighborsIndex.size());
	for (int i = 0; i < rhs.NeighborsIndex.size(); ++i)
	{
		
		for (int j = 0; j < 3; ++j)
		{
			this->NeighborsIndex[i].ind[j] = rhs.NeighborsIndex[i].ind[j];
		}
	}
	this->GaussSeidelTolerance = rhs.GaussSeidelTolerance;
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
	this->initilized = false;
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
		if (this->IsIndependent(it->first))
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
void NodeComposite::Solve_CellBasedTriangulation()
{
	double alpha = 0.5;
	Vector convergence(this->K->GetDim());
	int counter = 0;
	do
	{
		Vector X_old = (*this->X);
		for (int i = 0; i < this->K->GetDim(); ++i)
		{
			double sum = (*(this->C))(i);
			int n[3]{ 0,0,0 };

			for (int jj = 0; jj < 3; ++jj)
			{
				int j = this->NeighborsIndex[i].ind[jj];
				sum -= (j != -1) ? (*this->K)(i, j) * (*this->X)(j) : 0;
			}
			(*this->X)(i) = sum / (*this->K)(i, i);
		}
		(*this->X) = (*this->X) * alpha + X_old * (1.0 - alpha);
		convergence = (*this->X) - X_old;
		++counter;
	} while (counter < 1000 && (counter < 1 || convergence.abs() > this->K->GetGaussSeidelTolerance()));
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
NodeComposite& NodeComposite::operator=(const NodeComposite& rhs)
{
	this->Reset();
	this->CopyBody(rhs);
	return *this;
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
	this->initilized = true;
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
void NodeComposite::AddToValue(GeoGraphObject* objPtr, double value)
{
	Node* n = objPtr ? this->GetNode(objPtr) : 0;
	if (n)
		n->value += value;
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
void NodeComposite::SetValueAllNodes(function<double(const Vector3D& P)> func, bool basedOnVertices)
{
	VertexIterator itv(this->qe);
	Vertex* v = itv.Next();
	while (v) {
		Node* n = this->GetNode(v);
		if (n) {
			Vector3D P = n->GetPoint();
			double val = func(P);
			n->value = val;
		}
		v = itv.Next();
	}
	if(basedOnVertices)
		this->populate_basedOn_VERTICES();
	else
	{
		EdgeIterator ite(this->qe);
		Edge* e = ite.Next();
		while (e)
		{
			Node* n = this->GetNode(e);
			if (n) {
				Vector3D P = n->GetPoint();
				double val = func(P);
				n->value = val;
			}
			e = ite.Next();
		}
		FaceIterator itf(this->qe);
		Face* face = itf.Next();
		while (face)
		{
			Node* n = this->GetNode(face);
			if (n) {
				Vector3D P = n->GetPoint();
				double val = func(P);
				n->value = val;
			}
			face = itf.Next();
		}
	}
}
void NodeComposite::SetSolveMethod(METHOD value, double tolerance)
{
	this->solveMethod = value;
	if (this->solveMethod == METHOD::TRIANGULATION_CELL_BASED
		&& this->type != NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY
		&& this->type != NODE_COMPOSITE_TYPE::CELLS)
	{
		this->solveMethod = METHOD::GAUSS_SEIDEL;
	}
	if (this->solveMethod == METHOD::GAUSS_SEIDEL || this->solveMethod == METHOD::TRIANGULATION_CELL_BASED)
	{
		this->GaussSeidelTolerance = tolerance;
	}
	if (this->solveMethod == METHOD::TRIANGULATION_CELL_BASED)
	{
		this->NeighborsIndex.resize(this->index);
		FaceIterator itf(this->qe);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->fb)
			{
				int n = this->NodeHash[f]->index;
				EdgeIterator ite(f);
				for (int i = 0; i < 3; ++i)
				{
					Edge* e = ite.Next();
					Face* f1 = e->GetOtherFace(f);
					if (f1 != fb)
					{
						this->NeighborsIndex[n].ind[i] = this->NodeHash[f1]->index;
					}
					else if(this->type == NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY)
					{
						this->NeighborsIndex[n].ind[i] = this->NodeHash[e]->index;
					}
				}
			}
			f = itf.Next();
		}
		if (this->type == NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY)
		{
			EdgeIterator ite(this->fb);
			Edge* e = ite.Next();
			while (e)
			{
				int n = this->NodeHash[e]->index;
				Face* f1 = e->GetOtherFace(this->fb);
				this->NeighborsIndex[n].ind[0] = this->NodeHash[f1]->index;
				e = ite.Next();
			}
			int test = 0;
		}
	}
}
void NodeComposite::StabilizeK()
{
	/*To fix round up erros where aP is supposed to be equal or higher(in abs value) 
		than a_nb as per the numerical reciepe of the discretization method.*/
	if (this->K)
	{
		for (int i = 0; i < this->K->GetDim(); ++i)
		{
			double aP = fabs((*(this->K))(i, i));
			double sum = 0;
			for (int j = 0; j < this->K->GetDim(); ++j)
			{
				if (i != j)
				{
					sum += fabs((*(this->K))(i, j));
				}
			}
			double delta = fabs(aP - sum);
			if (aP < sum)
			{
				double alpha = 1.01;
				(*(this->K))(i, i) += (*(this->K))(i, i) > 0 ? alpha * delta : -alpha * delta;
			}
		}
	}
}
void NodeComposite::Solve()
{
	this->ApplyBC_internal("Solve");
	bool fixRoundOffError = true;
	if (this->solveMethod == METHOD::GAUSS_SEIDEL && this->K->CanUse_GaussSeidel(fixRoundOffError))
	{
		this->K->SetGaussSeidelTolerance(this->GaussSeidelTolerance);
		*(this->X) = this->K->Solve(*C, SquareMatrix::METHOD::GAUSS_SEIDEL);
	}
	if (this->solveMethod == METHOD::TRIANGULATION_CELL_BASED && this->K->CanUse_GaussSeidel(fixRoundOffError))
	{
		this->Solve_CellBasedTriangulation();
	}
	else
	{
		*(this->X) = this->K->Solve(*C, SquareMatrix::METHOD::GAUSS_ELIMINATION);
	}
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
double NodeComposite::GetValue(GeoGraphObject* objPtr, bool& isValid)
{
	isValid = false;
	if (this->IsInitialized())
	{
		Node* n = this->NodeHash[objPtr];
		if (n)
		{
			isValid = true;
			return n->value;
		}
	}
	return 0.0;
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
double NodeComposite::GetK(GeoGraphObject* row, GeoGraphObject* col)
{
	if (!this->isDependent[row] && !this->isDependent[col])
	{
		int i = this->NodeHash[row]->index;
		int j = this->NodeHash[col]->index;
		return (*K)(i, j);
	}
	return 0;
}
void NodeComposite::SetK(const SquareMatrix& Kr)
{
	*(this->K) = Kr;
}
void NodeComposite::SetK(const NodeComposite& NC)
{
	*(this->K) = *(NC.K);
}
double NodeComposite::GetC(GeoGraphObject* row)
{
	if (!this->isDependent[row])
	{
		int i = this->NodeHash[row]->index;
		return (*C)(i);
	}
	return 0.0;
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
bool NodeComposite::IsInitialized() const
{
	return this->initilized;
}
bool NodeComposite::IsConstValue(GeoGraphObject* objPtr)
{
	auto constValueIter = this->constValues.find(objPtr);
	return constValueIter != this->constValues.end();
}