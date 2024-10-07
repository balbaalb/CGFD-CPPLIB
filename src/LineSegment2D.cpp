#include "MathUtils.h"
#include "LineSegment2D.h"
void LineSegment2D::CopyBody(const LineSegment2D& rhs)
{
	for (int i = 0; i < 2; ++i)
	{
		this->points[i] = rhs.points[i];
	}
}
LineSegment2D::LineSegment2D()
{

}
LineSegment2D::LineSegment2D(const Vector3D& A, const Vector3D& B)
{
	this->points[0] = A;
	this->points[1] = B;
}
LineSegment2D::LineSegment2D(const LineSegment2D& rhs)
{
	this->CopyBody(rhs);
}
LineSegment2D& LineSegment2D::operator=(const LineSegment2D& rhs)
{
	this->CopyBody(rhs);
	return *this;
}
double LineSegment2D::GetArea() const
{
	return 0;
}
double LineSegment2D::GetPerimeter() const
{
	return this->points[0].distance(this->points[1]);
}
Vector3D LineSegment2D::GetPoint(int i) const
{
	return this->points[i];
}
void LineSegment2D::SetPoint(int i, const Vector3D& p)
{
	this->points[i] = p;
}
bool LineSegment2D::isOn(const Vector3D& P) const
{
	if (this->GetRelationTo(P) == ON_LINE_SEGMENT)
		return true;
	return false;
}
bool LineSegment2D::isIntersecting(const LineSegment2D& L) const
{
	if (this->operator==(L))
		return true;
	/*Vector3D V = this->points[1] - this->points[0];
	Vector3D L0 = L.points[0] - this->points[0];
	Vector3D L1 = L.points[1] - this->points[0];
	double c0 = V.CrossProduct2D(L0);
	double c1 = V.CrossProduct2D(L1);
	Vector3D W = L.points[1] - L.points[0];
	Vector3D S0 = this->points[0] - L.points[0];
	Vector3D S1 = this->points[1] - L.points[0];
	double c2 = W.CrossProduct2D(S0);
	double c3 = W.CrossProduct2D(S1);
	if (fabs(c0) > 1.0e-10 && fabs(c1) > 1.0e-10 && c0 * c1 < 0 && c2 * c3 < 0)
		return true;
	//Needs expansion*/
	Vector3D A = this->GetPoint(0);
	Vector3D B = this->GetPoint(1);
	Vector3D C = L.GetPoint(0);
	Vector3D D = L.GetPoint(1);
	if (this->GetRelationTo(C) == ON_LINE_SEGMENT && this->GetRelationTo(D) == ON_LINE_SEGMENT)
		return true;
	if (L.GetRelationTo(A) == ON_LINE_SEGMENT && L.GetRelationTo(B) == ON_LINE_SEGMENT)
		return true;
	if (this->GetRelationTo(C) == ON_LINE_SEGMENT && C != A && C != B)
		return true;
	if (this->GetRelationTo(D) == ON_LINE_SEGMENT && D != A && D != B)
		return true;
	if (L.GetRelationTo(A) == ON_LINE_SEGMENT && A != C && A != D)
		return true;
	if (L.GetRelationTo(B) == ON_LINE_SEGMENT && B != C && B != D)
		return true;
	if (this->GetRelationTo(C) != ON_LINE_SEGMENT && this->GetRelationTo(D) != ON_LINE_SEGMENT
		&& L.GetRelationTo(A) != ON_LINE_SEGMENT && L.GetRelationTo(B) != ON_LINE_SEGMENT
		&& this->GetRelationTo(C) != this->GetRelationTo(D) 
		&& L.GetRelationTo(A) != L.GetRelationTo(B))
		return true;
	//Error: 1 on line segment, another on extension but covers the line Example: A(0,0), B(1,0), C = A & D = (2,0)
	return false;
}
bool LineSegment2D::isIntersectingHorizontalRayFrom(const Vector3D& P)
{
	return this->isIntersectingRayFrom(P,0);
}
bool LineSegment2D::isIntersectingRayFrom(const Vector3D& P, double theta)
{//Does intersect the line from P(x,y)
	if (theta > 2.0 * pi)
		theta -= floor(theta / 2.0 / pi) * 2.0 * pi;
	if (theta == 2.0 * pi)
		theta = 0;
	Vector3D PA = this->points[0] - P;
	Vector3D PB = this->points[1] - P;
	double maxDist = std::fmax(PA.abs(), PB.abs());
	double farDistance = maxDist * 2.0;
	Vector3D Q;
	if (theta == pi / 2.0)
	{
		Q(0) = P(0);
		Q(1) = P(1) + farDistance;
	}
	else if (theta == 30 * pi / 2.0)
	{
		Q(0) = P(0);
		Q(1) = P(1) - farDistance;
	}
	else
	{
		double m = tan(theta);
		double b = P(1) - m * P(0);
		Q(0) = P(0) + maxDist * cos(theta);
		Q(1) = P(1) + maxDist * sin(theta);
	}
	LineSegment2D L2(P, Q);
	return this->isIntersecting(L2);
}
bool LineSegment2D::operator==(const LineSegment2D& rhs) const
{
	if (this->points[0] == rhs.points[0] && this->points[1] == rhs.points[1])
		return true;
	if (this->points[0] == rhs.points[1] && this->points[1] == rhs.points[0])
		return true;
	return false;
}
bool LineSegment2D::operator!=(const LineSegment2D& rhs) const
{
	return !(this->operator==(rhs));
}
POINT_LINESEGMENT2D_RELATION LineSegment2D::GetRelationTo(const Vector3D& P) const
{
	Vector3D L = this->points[1] - this->points[0];
	Vector3D LP = P - this->points[0];
	double LxP = L.CrossProduct2D(LP);
	if (LxP > 1.0e-10)
		return ON_LEFT;
	if (LxP < -1.0e-10)
		return ON_RIGHT;
	if (P == this->points[0] || P == this->points[1])
		return ON_LINE_SEGMENT;
	if ((L || LP) > -1.0e-10 && LP.abs() <= L.abs())
		return ON_LINE_SEGMENT;
	return ON_EXTENSION;
}