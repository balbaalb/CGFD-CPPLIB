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
Node& Node::operator=(const Node& rhs)
{
	this->host = rhs.host; 
	this->index = rhs.index;
	this->value = rhs.value;
	this->type = rhs.type;
	return *this;
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
NodeT& NodeT::operator=(const NodeT& rhs)
{
	this->Node::operator=(rhs);
	this->T = rhs.T;
	return *this;
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
NodeV& NodeV::operator=(const NodeV& rhs)
{
	this->Node::operator=(rhs);
	this->V = rhs.V;
	this->d = rhs.d;
	this->Horizontal = rhs.Horizontal;
	return *this;
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
NodeTV& NodeTV::operator=(const NodeTV& rhs)
{
	this->Node::operator=(rhs);
	this->T = rhs.T;
	this->vx = rhs.vx;
	this->vy = rhs.vy;
	this->Tprime = rhs.Tprime;
	return *this;
}