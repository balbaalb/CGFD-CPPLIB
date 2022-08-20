#ifndef FaceH
#define FaceH
#include "GeoGraphObject.h"
class Edge;
#include "Vector3D.h"
class Face : public GeoGraphObject
{
protected:
	Edge* edge;
	void CopyBody(const Face& f);
	int degree;
public:
	int type;//For debugging only
	int index;//For debugging only
	Face* next;
	Face* prev;
	Face();
	Face(const Face& f);
	Face& operator=(const Face& f);
	Edge* GetEdge() const;
	void SetEdge(Edge* e);
	int GetDegree() const;
	int IncrementDegree();
	int DecrementDegree();
	int SetDegree(int d);
	Vector3D GetCenteriod();
	Vector3D GetPoint();
	virtual Face* Clone() const = 0;
};
class FaceBasic : public Face
{
public:
	FaceBasic();
	Face* Clone() const;
};
bool tester_Face(int& NumTests);
#endif