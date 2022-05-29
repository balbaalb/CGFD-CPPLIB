#include "Node.h"
#include "GeoGraphObject.h"
#include "Vertex.h" 
Node::Node(GeoGraphObject* Host)
{
	this->host = Host;
	this->index = -1;
	this->value = 0;
}
Node::Node(const Node& rhs)
{
	(*this) = rhs;
}
void Node::operator=(const Node& rhs)
{
	this->host = rhs.host; 
	this->index = rhs.index;
	this->value = rhs.value;
	this->type = rhs.type;
}
GeoGraphObject* Node::GetHost()
{
	return this->host;
}
void Node::SetHost(GeoGraphObject* Host)
{
	this->host = Host;
}
Vector3D Node::GetPoint()
{
	return this->host->GetPoint();
}
NodeT::NodeT(GeoGraphObject* Host) : Node(Host)
{
	this->T = 0;
}
NodeT::NodeT(const NodeT& rhs) : Node(rhs)
{
	this->T = rhs.T;
}
NodeT::NodeT(const Node& rhs) : Node(rhs)
{
	this->T = 0;
}
void NodeT::operator=(const NodeT& rhs)
{
	this->Node::operator=(rhs);
	this->T = rhs.T;
}
NodeV::NodeV(GeoGraphObject* Host) : Node(Host)
{
	this->V = 0;
	this->d = 0;
	this->Horizontal = false;
}
NodeV::NodeV(const NodeV& rhs) : Node(rhs)
{
	this->V = rhs.V;
	this->d = rhs.d;
	this->Horizontal = rhs.Horizontal;
}
void NodeV::operator=(const NodeV& rhs)
{
	this->Node::operator=(rhs);
	this->V = rhs.V;
	this->d = rhs.d;
	this->Horizontal = rhs.Horizontal;
}
NodeTV::NodeTV(GeoGraphObject* Host) : Node(Host)
{
	this->T = 0;
	this->vx = 0;
	this->vy = 0;
	this->Tprime = 0;
}
NodeTV::NodeTV(const NodeTV& rhs) : Node(rhs)
{
	this->T = rhs.T;
	this->vx = rhs.vx;
	this->vy = rhs.vy;
	this->Tprime = rhs.Tprime;
}
void NodeTV::operator=(const NodeTV& rhs)
{
	this->Node::operator=(rhs);
	this->T = rhs.T;
	this->vx = rhs.vx;
	this->vy = rhs.vy;
	this->Tprime = rhs.Tprime;
}
bool tester_Node(int& NumTests)
{
	Vector3D P(1.5, 3.5);
	VertexBasic* v = new VertexBasic(P);
	Node n0(v);
	Node n1(n0), n2(0);
	n2 = n1;
	if (n2.GetPoint() != P)
		return false;
	NodeT n3(v);
	n3.T = 10.5;
	NodeT n4(n3), n5(0);
	n5 = n4;
	if (n5.T != 10.5)
		return false;
	//---------------
	NodeTV n6(v);
	n6.T = 10.5;
	n6.vx = -1.5;
	n6.vy = 8.13;
	NodeTV n7(n6), n8(0);
	n8 = n7;
	if (n8.T != 10.5)
		return false;
	if (n8.vx != -1.5)
		return false;
	if (n8.vy != 8.13)
		return false;
	++NumTests;
	return true;
}