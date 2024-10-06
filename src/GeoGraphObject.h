#ifndef GeoGraphObjectH
#define GeoGraphObjectH
#include "Vector3D.h"
class GeoGraphObject
{
public:
	virtual Vector3D GetPoint() = 0;
};
#endif
