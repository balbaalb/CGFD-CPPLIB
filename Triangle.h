#ifndef TriangleH
#define TriangleH
#include "Shape2D.h"
#include "LineSegment2D.h"
#include "SquareMatrix.h"
class Triangle : public Shape2D
{
	Vector3D points[3];
	void CopyBody(const Triangle& rhs);
public:
	Triangle();
	Triangle(const Vector3D& A, const Vector3D& B, const Vector3D& C);
	Triangle(const Triangle& rhs);
	void operator=(const Triangle& rhs);
	double GetArea() const;
	double GetPerimeter() const;
	Vector3D GetPoint(int i) const;
	void SetPoint(int i, const Vector3D& p);
	void GetLineSegments(LineSegment2D& e0, LineSegment2D& e1, LineSegment2D& e2) const;
	bool isInside(const Vector3D& p) const;
	bool isOnOrInside(const Vector3D& p) const;
	bool isIntersecting(const Triangle& T) const;
	bool operator==(const Triangle& rhs) const;
	bool operator!=(const Triangle& rhs) const;
	bool isCongruent(const Triangle& rhs) const;
	double GetShapeFactor() const;
	double GetMaxAngle() const;
	double GetMinAngle() const;
	double GetPiAngles() const;
	double GetIncircleRadius() const;
};
bool tester_Triangle(int& NumTests);
#endif
