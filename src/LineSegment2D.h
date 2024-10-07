#ifndef LineSegment2DH
#define LineSegment2DH
#include "Shape2D.h"
enum POINT_LINESEGMENT2D_RELATION
{
	ON_LEFT,
	ON_RIGHT,
	ON_EXTENSION,
	ON_LINE_SEGMENT
};
class LineSegment2D : public Shape2D
{
	Vector3D points[2];
	void CopyBody(const LineSegment2D& rhs);
public:
	LineSegment2D();
	LineSegment2D(const Vector3D& A, const Vector3D& B);
	LineSegment2D(const LineSegment2D& rhs);
	LineSegment2D& operator=(const LineSegment2D& rhs);
	double GetArea() const;
	double GetPerimeter() const;
	Vector3D GetPoint(int i) const;
	void SetPoint(int i, const Vector3D& p);
	bool isOn(const Vector3D& P) const;
	bool isIntersecting(const LineSegment2D& L) const;
	bool isIntersectingHorizontalRayFrom(const Vector3D& P);//Does intersect the line from P(x,y) to (inf,y)
	bool isIntersectingRayFrom(const Vector3D& P, double theta);//Does intersect the line from P(x,y) to (inf,y)
	bool operator==(const LineSegment2D& rhs) const;
	bool operator!=(const LineSegment2D& rhs) const;
	POINT_LINESEGMENT2D_RELATION GetRelationTo(const Vector3D& p) const;
};
#endif