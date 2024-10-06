#include "EdgeList.h"
#include "Edge.h"
EdgeContainer::EdgeContainer(Edge* e)
{
	this->edge = e;
	this->next = 0;
	this->prev = 0;
	this->index = -1;
}
EdgeContainer::EdgeContainer(const EdgeContainer& E)
{
	this->edge = E.edge;
	this->next = 0;
	this->prev = 0;
	this->index = -1;
}
EdgeContainer& EdgeContainer::operator=(const EdgeContainer& E)
{
	this->edge = E.edge;
	this->index = -1;
	return *this;
}
void EdgeList::CopyBody(const EdgeList& rhs)
{
	this->edgeContainerHash.rehash(rhs.edgeContainerHash.size() * 2);
	if (rhs.container0)
	{
		EdgeContainer* E = rhs.container0;
		do
		{
			this->AddContainer(E->edge);
			E = E->next;
		} while (E != rhs.container0);
	}
}
EdgeList::EdgeList(int Order)
{
	this->container0 = 0;
	this->NumEdgeContainers = 0;
	this->edgeContainerHash.rehash(Order);
}
EdgeList::EdgeList(const EdgeList& rhs)
{
	this->container0 = 0;
	this->NumEdgeContainers = 0; 
	this->CopyBody(rhs);
}
EdgeList::~EdgeList()
{
	this->Reset();
}
EdgeList& EdgeList::operator=(const EdgeList& rhs)
{
	this->Reset();
	this->CopyBody(rhs);
	return *this;
}
EdgeContainer* EdgeList::GetContainer(int ei) const
{
	if (ei < 0 || (ei > 0 && ei >= this->NumEdgeContainers))
		throw "EdgeList::GetContainer() Error";
	EdgeContainer* E = this->container0;
	for (int i = 0; i < ei; ++i)
	{
		if (!E)
			throw "EdgeList::GetContainer() Error";
		E = E->next;
	}
	return E;
}
int EdgeList::GetNumEdgeContainers() const
{
	return this->NumEdgeContainers;
}
void EdgeList::AddContainer(Edge* e)
{
	++this->NumEdgeContainers;
	if (!this->container0)
	{
		this->container0 = new EdgeContainer(e);
		this->edgeContainerHash[e] = this->container0;
		this->container0->next = this->container0;
		this->container0->prev = this->container0;
	}
	else
	{
		EdgeContainer* E1 = this->container0->prev;
		this->container0->prev = new EdgeContainer(e);
		this->edgeContainerHash[e] = this->container0->prev;
		this->container0->prev->next = this->container0;
		this->container0->prev->prev = E1;
		E1->next = this->container0->prev;
	}
}
void EdgeList::RemoveContainer(EdgeContainer* E)
{
	EdgeContainer* En = E->next;
	EdgeContainer* Ep = E->prev;
	if (this->container0 == E)
		this->container0 = En;//what if En == E?
	delete E;
	--this->NumEdgeContainers;
	if (!this->NumEdgeContainers)
		this->container0 = 0;
	else
	{
		Ep->next = En;
		En->prev = Ep;
	}
}
void EdgeList::RemoveContainerAndDeleteEdge(EdgeContainer* E)
{
	delete E->edge;
	E->edge = 0;
	this->RemoveContainer(E);
}
void EdgeList::Reset()
{
	if (this->container0 && this->container0->prev)
		this->container0->prev->next = 0;
	while (this->container0)
	{
		EdgeContainer* E1 = this->container0->next;
		delete this->container0;
		this->container0 = E1;
	}
	this->NumEdgeContainers = 0;
	this->container0 = 0;
	this->edgeContainerHash.clear();
}
void EdgeList::ResetAndDeleteEdges()
{
	if (this->container0)
	{
		EdgeContainer* container = this->container0;
		do
		{
			Edge* e = container->edge;
			delete e;
			container = container->next;
		} while (container != this->container0);
	}
	this->Reset();
}
EdgeContainer* EdgeList::FindContainer(Edge* e)
{
	if (this->container0)
	{
		return this->edgeContainerHash[e];
	}
	return 0;
}
void EdgeList::MoveToEnd(EdgeContainer* E)
{
	if (E != this->container0->prev)
	{
		if (E == this->container0)
		{
			this->container0 = this->container0->next;
		}
		else
		{
			EdgeContainer* En = E->next;
			EdgeContainer* Ep = E->prev;
			Ep->next = En;
			En->prev = Ep;
			EdgeContainer* E0p = this->container0->prev;
			this->container0->prev = E;
			E->next = this->container0;
			E->prev = E0p;
			E0p->next = E;
		}
	}
}
bool tester_EdgeList(int& NumTests)
{
	EdgeList eL;
	EdgeBasic e0, e1, *e2;
	e2 = new EdgeBasic;
	eL.AddContainer(&e0);
	eL.AddContainer(&e1);
	eL.AddContainer(e2);
	if (eL.GetNumEdgeContainers() != 3)
		return false;
	EdgeContainer* E0 = eL.GetContainer(0);
	EdgeContainer* E1 = eL.GetContainer(1);
	EdgeContainer* E2 = eL.GetContainer(2);
	if (E0->edge != &e0 || E0->next != E1 || E0->prev != E2)
		return false;
	if (E1->edge != &e1 || E1->next != E2 || E1->prev != E0)
		return false;
	if (E2->edge != e2 || E2->next != E0 || E2->prev != E1)
		return false;
	if (eL.FindContainer(&e0) != E0 || eL.FindContainer(&e1) != E1 || eL.FindContainer(e2) != E2)
		return false;
	//---------------------------------------------------------
	EdgeList eL2(eL), eL3;
	eL3 = eL2;
	if (eL3.GetNumEdgeContainers() != 3)
		return false;
	EdgeContainer* E0p = eL3.GetContainer(0);
	EdgeContainer* E1p = eL3.GetContainer(1);
	EdgeContainer* E2p = eL3.GetContainer(2);
	if (E0p->edge != &e0 || E0p->next != E1p || E0p->prev != E2p)
		return false;
	if (E1p->edge != &e1 || E1p->next != E2p || E1p->prev != E0p)
		return false;
	if (E2p->edge != e2 || E2p->next != E0p || E2p->prev != E1p)
		return false;
	//---------------------------------------------------------
	eL.RemoveContainer(E1);
	if (eL.GetContainer(0) != E0)
		return false;
	if (eL.GetNumEdgeContainers() != 2)
		return false;
	if (E0->edge != &e0 || E0->next != E2 || E0->prev != E2)
		return false;
	if (E2->edge != e2 || E2->next != E0 || E2->prev != E0)
		return false;
	//---------------------------------------------------------
	eL.RemoveContainerAndDeleteEdge(E2);
	if (eL.GetContainer(0) != E0)
		return false;
	if (eL.GetNumEdgeContainers() != 1)
		return false;
	if (E0->edge != &e0 || E0->next != E0 || E0->prev != E0)
		return false;
	//---------------------------------------------------------
	eL.RemoveContainer(E0);
	if (eL.GetContainer(0))
		return false;
	if (eL.GetNumEdgeContainers())
		return false;
	//---------------------------------------------------------
	e2 = new EdgeBasic;
	eL.AddContainer(&e0);
	eL.AddContainer(&e1);
	eL.AddContainer(e2);
	if (eL.GetNumEdgeContainers() != 3)
		return false;
	E0 = eL.GetContainer(0);
	E1 = eL.GetContainer(1);
	E2 = eL.GetContainer(2);
	//---------------------------------------------------------
	eL.RemoveContainerAndDeleteEdge(E2);
	if (eL.GetContainer(0) != E0)
		return false;
	if (eL.GetNumEdgeContainers() != 2)
		return false;
	if (E0->edge != &e0 || E0->next != E1 || E0->prev != E1)
		return false;
	if (E1->edge != &e1 || E1->next != E0 || E1->prev != E0)
		return false;
	//---------------------------------------------------------
	eL.Reset();
	if (eL.GetContainer(0))
		return false;
	if (eL.GetNumEdgeContainers())
		return false;
	//---------------------------------------------------------
	EdgeList eL4;
	e2 = new EdgeBasic;
	EdgeBasic e3;
	eL4.AddContainer(&e0);
	eL4.AddContainer(&e1);
	eL4.AddContainer(e2);
	eL4.AddContainer(&e3);
	if (eL4.GetNumEdgeContainers() != 4)
		return false;
	E0 = eL4.GetContainer(0);
	E1 = eL4.GetContainer(1);
	E2 = eL4.GetContainer(2);
	EdgeContainer* E3 = eL4.GetContainer(3);
	eL4.MoveToEnd(E1);
	if (eL4.GetContainer(0) != E0 || eL4.GetContainer(1) != E2 || eL4.GetContainer(2) != E3 || eL4.GetContainer(3) != E1)
		return false;
	eL4.MoveToEnd(E1);
	if (eL4.GetContainer(0) != E0 || eL4.GetContainer(1) != E2 || eL4.GetContainer(2) != E3 || eL4.GetContainer(3) != E1)
		return false;
	eL4.MoveToEnd(E0);
	if (eL4.GetContainer(0) != E2 || eL4.GetContainer(1) != E3 || eL4.GetContainer(2) != E1 || eL4.GetContainer(3) != E0)
		return false;
	eL4.MoveToEnd(E1);
	if (eL4.GetContainer(0) != E2 || eL4.GetContainer(1) != E3 || eL4.GetContainer(2) != E0 || eL4.GetContainer(3) != E1)
		return false;
	eL4.MoveToEnd(E3);
	if (eL4.GetContainer(0) != E2 || eL4.GetContainer(1) != E0 || eL4.GetContainer(2) != E1 || eL4.GetContainer(3) != E3)
		return false;
	NumTests += 1;
	return true;
}