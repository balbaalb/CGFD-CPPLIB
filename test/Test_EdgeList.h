#ifndef Test_EdgeListH
#define Test_EdgeListH
#include "../src/EdgeList.h"
#include "../src/Edge.h"
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
#endif