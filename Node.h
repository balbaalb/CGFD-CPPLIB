#ifndef NodeH
#define NodeH
#include <string>
#include "Vector3D.h"
using namespace std;
class GeoGraphObject;
class Node
{
	void virtual dummy() {};
	GeoGraphObject* host;
public:
	int index;
	double value;
	string type;
	Node(GeoGraphObject* Host = 0);
	Node(const Node& rhs);
	Node& operator=(const Node& rhs);
	GeoGraphObject* GetHost();
	void SetHost(GeoGraphObject* Host);
	Vector3D GetPoint();
};
class NodeT : public virtual Node
{
public:
	double T;
	NodeT(GeoGraphObject* Host = 0);
	NodeT(const NodeT& rhs);
	NodeT(const Node& rhs);
	NodeT& operator=(const NodeT& rhs);
};
class NodeV : public virtual Node
{
public:
	double V;
	bool Horizontal;
	double d;//for SIMPLE equations 
	NodeV(GeoGraphObject* Host = 0);
	NodeV(const NodeV& rhs);
	NodeV& operator=(const NodeV& rhs);
};
class NodeTV : public virtual Node
{
public:
	double T;
	double vx;
	double vy;
	double Tprime;
	NodeTV(GeoGraphObject* Host = 0);
	NodeTV(const NodeTV& rhs);
	NodeTV& operator=(const NodeTV& rhs);
};
bool tester_Node(int& NumTests);
#endif