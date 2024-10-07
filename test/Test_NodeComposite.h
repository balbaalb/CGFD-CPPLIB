#ifndef Test_NodeCompositeH
#define Test_NodeCompositeH
#include "../src/NodeComposite.h"
#include "../src/bGrid.h"
#include "../src/QuadEdge.h"
#include "../src/EdgeIterator.h"
#include "../src/Edge.h"
#include "../src/Face.h"
#include "../src/FaceIterator.h"
#include "../src/Vertex.h"
bool tester_NodeComposite(int& NumTests) 
{
	double A = 1.65, B = 8.75, T0 = 10.0;//T = Ax + By + T0
	bGrid grid(3, 2, 3.0, 2.0);
	QuadEdge* qe = grid.GetMesh2D();
	Face* fb = grid.GetBoundary();
	NodeComposite TempNodes;
	TempNodes.Initialize(qe, fb, CELLS_AND_BOUNDARY);
	//---- Setting Boundary conditions ---------------------
	EdgeIterator iteb(fb);
	Edge* eb = iteb.Next();
	while (eb)
	{
		double x = eb->GetPoint()(0);
		double y = eb->GetPoint()(1);
		double T = T0 + A * x + B * y;
		TempNodes.SetConstantValue(eb, T);
		eb = iteb.Next();
	}
	//---- Write Equations, boundary conditions will be applied automatically from constant values. 
	TempNodes.InitializeEquations();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb) {
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			double a = 0;
			while (e) 
			{
				Face* f1 = e->GetOtherFace(f);
				if (f1 != fb)
				{
					TempNodes.AddToK(f, f1, -5.0);
					a += 5.0;
				}
				else 
				{
					TempNodes.AddToK(f, e, -10.0);
					a += 10;
				}
				e = ite.Next();
			}
			TempNodes.AddToK(f, f, a);
		}
		f = itf.Next();
	}
	//------- Solve ---------------------------------------------
	TempNodes.SolveAndUpdate();
	TempNodes.Populate();
	//------ Post Processing ------------------------------------
	Vertex* v7 = qe->GetVertex(7);
	Node* nv7 = TempNodes.GetNode(v7);
	if (!nv7)
		return false;
	double Tv7 = nv7->value;
	Vector3D Pv7 = v7->GetPoint();
	Vector3D Pnv7 = nv7->GetPoint();
	if (Pv7 != Pnv7)
		return false;
	double Tv7_correct = A * Pv7(0) + B * Pv7(1) + T0;
	double error_v7 = fabs(Tv7 - Tv7_correct) / Tv7_correct * 100;
	if (error_v7 > 0.01)
		return false;
	double Tv7a = TempNodes.GetValue(Pnv7);
	if (fabs(A * Pnv7(0) + B * Pnv7(1) + T0 - Tv7a) > 1.0e-10)
		return false;
	Edge* e3 = qe->GetEdge(3);
	Vector3D PA = e3->GetOrig()->GetPoint();
	Vector3D PB = e3->GetDest()->GetPoint();
	Vector3D C = (PA * 2.0 / 3.0) + (PB * 1.0 / 3.0);
	double TA = TempNodes.GetValue(e3->GetOrig());
	double TB = TempNodes.GetValue(e3->GetDest());
	double TC = TempNodes.GetValue(C);
	if (fabs(A * C(0) + B * C(1) + T0 - TC) > 1.0e-10)
		return false;
	Vector3D D(1.2, 0.8);
	double TD = TempNodes.GetValue(D);
	if (fabs(A * D(0) + B * D(1) + T0 - TD) > 1.0e-10)
		return false;
	Vector3D E(2.5, 1.25);
	double TE = TempNodes.GetValue(E);
	if (fabs(A * E(0) + B * E(1) + T0 - TE) > 1.0e-10)
		return false;
	Vector3D F(2.5, 1.75);
	double TF = TempNodes.GetValue(F);
	if (fabs(A * F(0) + B * F(1) + T0 - TF) > 1.0e-10)//TF becomes 28.785 because of bug in calculating corner vertices
		return false;
	++NumTests;
	return true;
}
#endif