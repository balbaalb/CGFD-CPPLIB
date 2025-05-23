#ifndef Test_VertexH
#define Test_VertexH
#include "../src/Vertex.h"
#include "../src/MathUtils.h"
#include "../src/Edge.h"
bool tester_Vertex(int& NumTests)
{
	EdgeBasic e1;
	VertexBasic v1(11.5, 8.5, -135.6);
	VertexBasic v2(v1), v3;
	v3 = v2;
	v3.SetEdge(&e1);
	if (!v3.GetEdge() || v3.GetEdge() != &e1)
		return false;
	Vector3D u3 = v3.GetPoint();
	if (IsNotEqual(u3(0), 11.5) || IsNotEqual(u3(1), 8.5) || IsNotEqual(u3(2), -135.6))
		return false;
	Vector3D p(12.3, 4.5, 6);
	VertexBasic v4(p);
	VertexBasic v5(v4), v6;
	v6 = v5;
	if (v4.GetEdge() || v5.GetEdge() || v6.GetEdge())
		return false;
	VertexBasic v7;
	v7 << p;
	if (v7.GetEdge())
		return false;
	if (v7.GetPoint() != v6.GetPoint() || !(v6.GetPoint() == v7.GetPoint()))
		return false;
	Vector3D q1 = v6.GetPoint();
	Vector3D q2 = v7.GetPoint();
	if (p != q1 || p != q2 || q1 != q2)
		return false;
	VertexBasic v8;
	if (v8.GetDegree())
		return false;
	if (v8.IncreaseDegree() != 1 || v8.IncreaseDegree() != 2)
		return false;
	if (v8.SetDegree(10) != 10 || v8.DecreaseDegree() != 9)
		return false;
	Vector3D q3(11.5, 13.5, -4.65);
	VertexBasic* v9 = new VertexBasic(q3);
	v9->SetDegree(5);
	VertexBasic* v10 = new VertexBasic(*v9);
	VertexBasic* v11 = new VertexBasic;
	*v11 = *v10;
	Vector3D q4 = v11->GetPoint();
	if (q4(0) != 11.5 || q4(1) != 13.5 || q4(2) != -4.65)
		return false;
	if (v11->GetDegree())
		return false;
	Vector3D q5(-1.25, 12.18, 13.175);
	(*v11) << q5;
	Vector3D q6 = v11->GetPoint();
	if (IsNotEqual(q6(0) , -1.25) || IsNotEqual(q6(1) , 12.18) || IsNotEqual(q6(2) , 13.175))
		return false;
	delete v9;
	delete v10;
	delete v11;
	NumTests += 1;
	return true;
}
#endif