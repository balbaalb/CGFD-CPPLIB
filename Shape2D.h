#ifndef Shape2DH
#define Shape2DH
#include "Vector3D.h"
class Shape2D
{
public:
	virtual double GetArea() const = 0;
	virtual double GetPerimeter() const = 0;
};
#endif