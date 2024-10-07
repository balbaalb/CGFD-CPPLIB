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