#include "Element.h"
Element2D::Element2D(Face* f)
{
	this->face = f;
}
Element2D::Element2D(const Element2D& rhs)
{
	this->face = rhs.face;
}
void Element2D::operator=(const Element2D& rhs)
{
	this->face = rhs.face;
}
Cell2D1N::Cell2D1N(Face* f) : Element2D(f)
{

}
Cell2D1N::Cell2D1N(const Element2D& rhs) : Element2D(rhs)
{

}
Cell2D1N::Cell2D1N(const Cell2D1N& rhs) : Element2D(rhs)
{
	this->center = rhs.center;
}
void Cell2D1N::operator=(const Element2D& rhs)
{
	this->Element2D::operator=(rhs);
}
void Cell2D1N::operator=(const Cell2D1N& rhs)
{
	this->Element2D::operator=(rhs);
	this->center = rhs.center;
}
bool tester_Element(int& NumTests)
{
	Face* f = new FaceBasic;
	Element2D e0(f);
	Element2D e1(e0), e2;
	e2 = e1;
	if (e2.face != f)
		return false;
	Cell2D1N c0(e2);
	c0.center.T = 10.5;
	Cell2D1N c1(c0), c2;
	c2 = c1;
	if (c2.center.T != 10.5)
		return false;
	++NumTests;
	return true;
}