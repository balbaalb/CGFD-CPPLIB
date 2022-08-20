#ifndef LineH
#define LineH
#include "Vector3D.h"
class Line : public Vector3D //start Point + dx, dy, dz, infinitely extended in both directions
{
	Vector3D Start;
protected:
	void CopyBody(const Line& L);
public:
	Line();
	Line(const Vector3D& x0, const Vector3D& x1);
	Line(const Line& L);
	Line& operator=(const Line& rhs);
	void SetTwoPoints(const Vector3D& start, const Vector3D& end);
	void SetPointAndVector(const Vector3D& start, const Vector3D& dx);
	Vector3D GetStartPoint() const;
	bool IsOnSamePlane(const Line& L2) const;
	Vector3D project(const Vector3D& P) const;
	double distance(const Vector3D& P) const;
	Vector3D NearestPointTo(const Line& L2) const;
	double distance(const Line& L2) const;
	bool DoesIntersect(const Line& L2) const;
	Vector3D* intersection(const Line& L2) const;
	bool operator==(const Line& L2) const;
	bool operator!=(const Line& L2) const;
	//Note: Do NOT create Line(const Vector3D& P) 
};
bool tester_Line(int& NumTests);
#endif