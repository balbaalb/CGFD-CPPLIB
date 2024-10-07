#ifndef PlaneH
#define PlaneH
#include "Vector3D.h"
#include "Line.h"
class Plane : public Vector3D //normal to plane; 
{
	Vector3D Ref;//Refernce Point
protected:
	void CopyBody(const Plane& P);
public:
	Plane();
	Plane(const Plane& P);
	Plane(const Vector3D& n, double D);
	Plane(const Vector3D& n, const Vector3D& ref);
	Plane(const Vector3D& A, const Vector3D& B, const Vector3D& D);
	Plane& operator=(const Plane& P);
	void SetNormalAndConst(const Vector3D& n, double C);
	void SetNormalAndRef(const Vector3D& n, const Vector3D& ref);
	void SetThreePoints(const Vector3D& A, const Vector3D& B, const Vector3D& C);
	Vector3D GetRefPoint() const;
	Vector3D normal() const;
	Vector3D tangent() const;//OverWriting
	double GetConst() const;
	bool IsParallel(const Line& L) const;
	bool IsPerpendicular(const Line& L) const;
	bool IsParallel(const Plane& P2) const;
	bool IsPerpendicular(const Plane& P2) const;
	bool operator==(const Plane& P) const;
	bool operator!=(const Plane& P) const;
	bool IsOnPlane(const Vector3D& Point) const;
	Line* intersection(const Plane& P2) const;
	double distance(const Vector3D& Point) const;
	Vector3D project(const Vector3D& Point) const;
	Vector3D* intersection(const Line& L) const;
	Line* project(const Line& L) const;
	double distance(const Line& L) const;
	double distance(const Plane& P) const;
};
#endif

