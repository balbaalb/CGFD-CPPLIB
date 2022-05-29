#ifndef CircleH
#define CircleH
#include "Vector3D.h"
#include "Line.h"
class Circle //2D
{
	Vector3D center;
	double radius;
	void CopyBody(const Circle& c);
public:
	Circle();
	Circle(Vector3D center, double radius);
	Circle(const Circle& c);
	void operator=(const Circle& c);
	double Area() const;
	double Perimeter() const;
	void SetCenter(const Vector3D& Center);
	Vector3D GetCenter() const;
	void SetRadius(double R);
	double GetRadius() const;
	void Set3Points(const Vector3D& P1, const Vector3D& P2, const Vector3D& P3);
	void SetCenterAndTangent(const Vector3D& Center, const Line& tangent);
	bool operator==(const Circle& c) const;
	bool operator!=(const Circle& c) const;
	int NumIntersectionPoints(const Circle& c, Vector3D* P1 = 0, Vector3D* P2 = 0) const;
	int tangent(const Vector3D& point, Line* t1 = 0, Line* t2 = 0) const;
	bool isInside(const Vector3D& p, double tolerance = 1.0e-10) const;//if on the circle will return false
};
bool tester_Circle(int& NumTests);
#endif