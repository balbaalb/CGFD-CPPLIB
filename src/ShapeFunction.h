#ifndef ShapeFunctionH
#define ShapeFunctionH
#include "Triangle.h"
#include "Vector3D.h"
#include "SquareMatrix.h"
class ShapeFunction
{
	enum SHAPE_FUNCTION_TYPE
	{
		LINEAR,
		EXPONENTIAL,
		PURE_CONVECTION
	} type;
	SquareMatrix* N;
	double Xmax;
	Vector3D uh;
	Vector3D nh;
	Vector3D Centeriod;
	double Pe;
	void CopyBody(const ShapeFunction& rhs);
	void DeleteBody();
public:
	ShapeFunction(const Triangle& T);
	ShapeFunction(const Triangle& T, double Pe_,const Vector3D& uh);
	ShapeFunction(const Triangle& T, const Vector3D& uh);//for Pe --> inifinity
	ShapeFunction(const ShapeFunction& rhs);
	~ShapeFunction();
	ShapeFunction& operator=(const ShapeFunction& rhs);
	double GetValue(int i, const Vector3D& P) const;
	Vector3D Grad(int i, const Vector3D& P) const;
};
#endif