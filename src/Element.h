#ifndef ElementH
#define ElementH
#include "Node.h"
#include "Face.h"
class Element2D
{
public:
	Face* face;//does not own
	Element2D(Face* f = 0);
	Element2D(const Element2D& rhs);
	void operator=(const Element2D& rhs);
};
class Cell2D1N : public Element2D
{
public:
	NodeT center;
	Cell2D1N(Face* f = 0);
	Cell2D1N(const Element2D& rhs);
	Cell2D1N(const Cell2D1N& rhs);
	void operator=(const Element2D& rhs);
	void operator=(const Cell2D1N& rhs);
};
bool tester_Element(int& NumTests);
#endif