#include <fstream>
#include "FVM.h"
#include "MathUtils.h"
#include "Vector.h"
#include "Vector3D.h"
#include "SquareMatrix.h"
#include "QuadEdge.h"
#include "FaceIterator.h"
#include "VertexIterator.h"
#include "Vertex.h"
#include "Triangle.h"
#include "ShapeFunction.h"
#include "Edge.h"
#include "EdgeIterator.h"
#include "Face.h"
#include "ConductionConvectionProblem.h"
#include "TVD.h"
#include "DelaunayLifting.h"
void FVM::InitializeTNodes()
{
	if (!this->TNodes.IsInitialized())
	{
		QuadEdge* qe = this->triangulation->GetMesh2D();
		Face* fb = this->triangulation->GetBoundary();
		if (this->cell_geo == CELL_GEOMETRY::VERTEX_BASED)
			this->TNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::VERTICES);
		else
			this->TNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
		this->TNodes.InitializeEquations();
	}
}
void FVM::AddDiffusion_vertexBased(double conductivity)
{
	/*In Triangle ABC:
	T(x,y) = N_A(x,y) * T_A + N_B(x,y) * T_B + N_C(x,y) * T_C;
	N_i(x,y) = N_i0 * x + N_i1 * y + N_i2; where i = 0,1,2 corresponds to A,B,C
	dT/dx = N_00 * T_A + N_10 * T_B + N_20 * T_C;
	dT/dy = N_01 * T_A + N_11 * T_B + N_21 * T_C;
	grad(T) = dT/dx * ex + dT/dy * ey;
	grad(T) dotProduct n = dT/dx * nx + dT/dy * ny
		 = (nx * N_00 + ny * N_01) T_A + (nx * N_10 + ny * N_11) T_B + (nx * N_20 + ny * N_21) T_C
		 = (n dot N_A) * T_A + (n dot N_B) * T_B + (n dot N_C) * T_C
	*/
	double k = conductivity;
	Vector3D ez(0, 0, 1);
	int N = this->triangulation->NumVertices();
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			VertexIterator itv(f);
			Vertex* ver[3];
			ver[0] = itv.Next();
			ver[1] = itv.Next();
			ver[2] = itv.Next();
			Vector3D point[3];
			point[0] = ver[0]->GetPoint();
			point[1] = ver[1]->GetPoint();
			point[2] = ver[2]->GetPoint();
			Vector3D O = (point[0] + point[1] + point[2]) / 3.0;
			Triangle T(point[0], point[1], point[2]);
			ShapeFunction N(T);
			for (int i = 0; i < 3; ++i)
			{
				int ib = (i != 2) ? i + 1 : 0;
				int ic = (i != 0) ? i - 1 : 2;
				Vector3D AO = O - point[i];
				Vector3D P = (point[i] + point[ib]) / 2.0;
				Vector3D Q = (point[i] + point[ic]) / 2.0;
				Vector3D R = (O + P) / 2.0;
				Vector3D S = (O + Q) / 2.0;
				Vector3D AP = P - point[i];
				Vector3D OP = P - O;
				Vector3D OQ = Q - O;
				Vector3D nOP = ez && OP;
				Vector3D nOQ = OQ && ez;
				nOP(2) = 0;
				nOQ(2) = 0;
				nOP = nOP / nOP.abs();
				nOQ = nOQ / nOQ.abs();
				if (AO.CrossProduct2D(AP) > 0){
					nOP = nOP * (-1);
					nOQ = nOQ * (-1);
				}
				for (int ii = 0; ii < 3; ++ii)
				{
					double Kij = 0;
					Kij += k * (nOP || N.Grad(ii, R)) * OP.abs();
					Kij += k * (nOQ || N.Grad(ii, S)) * OQ.abs();
					this->TNodes.AddToK(ver[i], ver[ii], Kij);
				}
			}
		}
		f = itf.Next();
	}
}
void FVM::AddDiffusion_cellBased(NodeComposite& nodes, function<double(const Vector3D& P)> conductivity)
{
	Vector3D ez(0, 0, 1);
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			double aP = 0; 
			double b = 0;
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			while (e)
			{
				Vector3D M = e->GetMidPoint();
				Vector3D CM = M - C;
				Vector3D eta = e->GetVector();
				Vertex* Va = e->GetOrig();
				Vertex* Vb = e->GetDest();
				if (CM.CrossProduct2D(eta) < 0)
				{
					eta = eta * (-1);
					std::swap(Va, Vb);
				}
				double Phi_a = nodes.GetValue(Va);
				double Phi_b = nodes.GetValue(Vb);
				double etaL = eta.abs();
				Vector3D e_eta = eta / etaL;
				Face* f1 = e->GetOtherFace(f);
				Vector3D C1;
				GeoGraphObject* colPtr;
				if (f1 != fb)
				{
					C1 = f1->GetCenteriod();
					colPtr = f1;
				}
				else
				{
					C1 = (e->GetOrig()->GetPoint() + e->GetDest()->GetPoint()) / 2.0;
					colPtr = e;
				}
				Vector3D ksi = C1 - C;
				double ksiL = ksi.abs();
				Vector3D e_ksi = ksi / ksiL;
				Vector3D n = e_eta && ez;
				n(2) = 0;//avoiding numerical issues.
				double k = conductivity(M);
				double De = k / ksiL;
				double Diff = De * etaL / (n || e_ksi);
				nodes.AddToK(f, colPtr, Diff);
				aP += Diff;
				double cross_diff_coef = -k * (e_ksi || e_eta) / (n || e_ksi);
				double source_cross_diffusion = cross_diff_coef * (Phi_b - Phi_a);
				b += source_cross_diffusion;
				e = ite.Next();
			}
			nodes.AddToK(f, f, -aP);
			nodes.AddToC(f, -b);
		}
		f = itf.Next();
	}
}
void FVM::AddDiffusionAndConvection_vertexBased(function<double(const Vector3D& P)> conductivity, double density, function<double(const Vector3D& P)> V[2])
{
	/*In Triangle ABC:
	T(x,y) = N_A(x,y) * T_A + N_B(x,y) * T_B + N_C(x,y) * T_C;
	N_i(x,y) = N_i0 * exp(Pe*(X-Xmax)) + N_i1 * Y + N_i2; where i = 0,1,2 corresponds to A,B,C
	Pe = density * U_average / conductivity
	X = (r - ro) dot u_hat
	r = x ex + y ey
	ro = location of triangle centeroid
	u_hat: velocity direction
	Y = (r - ro) dot n
	n = ez cross u_hat
	Equation to solve:
	density * v dot grad(T) - div(conductivity * grad(T)) = 0
	*/
	Vector3D ez(0, 0, 1);
	double rho = density;
	int N = this->triangulation->NumVertices();
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			VertexIterator itv(f);
			Vertex* ver[3];
			ver[0] = itv.Next();
			ver[1] = itv.Next();
			ver[2] = itv.Next();
			Vector3D point[3];
			point[0] = ver[0]->GetPoint();
			point[1] = ver[1]->GetPoint();
			point[2] = ver[2]->GetPoint();
			Vector3D O = (point[0] + point[1] + point[2]) / 3.0;
			Triangle T(point[0], point[1], point[2]);
			double vx = (V[0](point[0]) + V[0](point[1]) + V[0](point[2])) / 3.0;
			double vy = (V[1](point[0]) + V[1](point[1]) + V[1](point[2])) / 3.0;
			Vector3D v(vx, vy);
			Vector3D vh = v / v.abs();
			double k = (conductivity(point[0]) + conductivity(point[1]) + conductivity(point[2])) / 3.0;
			ShapeFunction* N;
			if (k > 1.0e-10)
			{
				double Pe = rho * v.abs() / k;
				N = new ShapeFunction(T, Pe, vh);
			}
			else
			{
				N = new ShapeFunction(T, vh);
			}
			for (int i = 0; i < 3; ++i)
			{
				int ib = (i != 2) ? i + 1 : 0;
				int ic = (i != 0) ? i - 1 : 2;
				Vector3D AO = O - point[i];
				Vector3D P = (point[i] + point[ib]) / 2.0;
				Vector3D Q = (point[i] + point[ic]) / 2.0;
				Vector3D R = (O + P) / 2.0;
				Vector3D S = (O + Q) / 2.0;
				Vector3D AP = P - point[i];
				Vector3D OP = P - O;
				Vector3D OQ = Q - O;
				Vector3D nOP = ez && OP;
				Vector3D nOQ = OQ && ez;
				nOP(2) = 0;
				nOQ(2) = 0;
				nOP = nOP / nOP.abs();
				nOQ = nOQ / nOQ.abs();
				if (AO.CrossProduct2D(AP) > 0) {
					nOP = nOP * (-1);
					nOQ = nOQ * (-1);
				}
				Vector3D VO(V[0](O), V[1](O));
				Vector3D VP(V[0](P), V[1](P));
				Vector3D VQ(V[0](Q), V[1](Q));
				Vector3D VR(V[0](R), V[1](R));
				Vector3D VS(V[0](S), V[1](S));;
				double kO = conductivity(O);
				double kP = conductivity(P);
				double kQ = conductivity(Q);
				double kR = conductivity(R);
				double kS = conductivity(S);
				for (int ii = 0; ii < 3; ++ii)
				{
					double Kij = 0;
					Kij += -kO * (nOP || N->Grad(ii, O)) * OP.abs() / 6.0;
					Kij += -kR * (nOP || N->Grad(ii, R)) * OP.abs() / 6.0 * 4.0;
					Kij += -kP * (nOP || N->Grad(ii, P)) * OP.abs() / 6.0;
					Kij += -kO * (nOQ || N->Grad(ii, O)) * OQ.abs() / 6.0;
					Kij += -kS * (nOQ || N->Grad(ii, S)) * OQ.abs() / 6.0 * 4.0;
					Kij += -kQ * (nOQ || N->Grad(ii, Q)) * OQ.abs() / 6.0;
					Kij += rho * (nOP || VO) * N->GetValue(ii, O) * OP.abs() / 6.0;
					Kij += rho * (nOP || VR) * N->GetValue(ii, R) * OP.abs() / 6.0 * 4.0;
					Kij += rho * (nOP || VP) * N->GetValue(ii, P) * OP.abs() / 6.0;
					Kij += rho * (nOQ || VO) * N->GetValue(ii, O) * OQ.abs() / 6.0;
					Kij += rho * (nOQ || VS) * N->GetValue(ii, S) * OQ.abs() / 6.0 * 4.0;
					Kij += rho * (nOQ || VQ) * N->GetValue(ii, Q) * OQ.abs() / 6.0;
					this->TNodes.AddToK(ver[i], ver[ii], Kij);
				}
			}
			delete N;
		}
		f = itf.Next();
	}
}
Vector3D FVM::GetGradient(NodeComposite& nodes, Face* f)
{
	Vector3D C = f->GetCenteriod(); 
	bMatrix A(3, 2);
	Vector B(3);
	EdgeIterator ite(f);
	Edge* e = ite.Next();
	Face* fb = this->triangulation->GetBoundary();
	Node* n0 = nodes.GetNode(f);
	double Phi0 = n0->value;
	int row = 0;
	while (e)
	{
		Face* f1 = e->GetOtherFace(f); 
		Vector3D C1;
		GeoGraphObject* colPtr;
		if (f1 != fb)
		{
			C1 = f1->GetCenteriod();
			colPtr = f1;
		}
		else
		{
			C1 = (e->GetOrig()->GetPoint() + e->GetDest()->GetPoint()) / 2.0;
			colPtr = e;
		}
		A(row, 0) = C1(0) - C(0);
		A(row, 1) = C1(1) - C(1);
		double Phi = nodes.GetValue(colPtr);
		B(row) = Phi - Phi0;
		++row;
		e = ite.Next();
	}
	bMatrix ATA = A.Transpose() * A;
	double determinant = ATA(0, 0) * ATA(1, 1) - ATA(0, 1) * ATA(1, 0);
	bMatrix ATA_inv = ATA / determinant;
	std::swap(ATA_inv(1, 0), ATA_inv(0, 1));
	Vector X = ATA_inv * A.Transpose() * B;
	return Vector3D(X(0), X(1));
}
void FVM::AddConvection_cellBased(NodeComposite& nodes)
{
	double k = this->constantConductivity; 
	double rho = this->constantDensity;
	Vector3D ez(0, 0, 1);
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	TVD psi(TVD::MODE::VAN_ALBADA);
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			double b = 0;
			Vector3D grad_Phi = this->GetGradient(nodes, f);
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			while (e)
			{
				Vector3D M = e->GetMidPoint();
				Node* uxNode = this->VxNodes.GetNode(e);
				Node* uyNode = this->VyNodes.GetNode(e);
				double ux = uxNode ? uxNode->value : 0;
				double uy = uyNode ? uyNode->value : 0;
				Vector3D uM(ux, uy);
				Vector3D CM = M - C;
				Vector3D eta = e->GetVector();
				if (CM.CrossProduct2D(eta) < 0)
				{
					eta = eta * (-1);
				}
				double etaL = eta.abs();
				Vector3D e_eta = eta / etaL;
				Face* f1 = e->GetOtherFace(f);
				Vector3D C1;
				GeoGraphObject* colPtr;
				if (f1 != fb)
				{
					C1 = f1->GetCenteriod();
					colPtr = f1;
				}
				else
				{
					C1 = (e->GetOrig()->GetPoint() + e->GetDest()->GetPoint()) / 2.0;
					colPtr = e;
				}
				Vector3D ksi = C1 - C;
				Vector3D n = e_eta && ez;
				n(2) = 0;//avoiding numerical issues.
				double un = uM || n;
				double flux = rho * un * etaL;
				double Phi0 = nodes.GetNode(f)->value;
				double Phi1 = nodes.GetNode(colPtr)->value;
				double Phi_downstream{0}, Phi_upstream{0};
				if (flux > 0)
				{
					nodes.AddToK(f, f, -flux);
					Phi_downstream = Phi1;
					Phi_upstream = Phi0;
				}
				else
				{
					nodes.AddToK(f, colPtr, -flux);
					Phi_downstream = Phi0;
					Phi_upstream = Phi1;
				}
				double dPhi = Phi_downstream - Phi_upstream;
				double r = fabs(dPhi) > 1.0e-10 ? 2 * (grad_Phi || ksi) / dPhi - 1 : 0;
				b += -flux * psi(r) / 2.0 * dPhi;
				e = ite.Next();
			}
			nodes.AddToC(f, -b);
		}
		f = itf.Next();
	}
}
void FVM::SolveSIMPLE_v(int i)
{
	NodeComposite* nodes = (i == 0) ? &(this->VxNodes) : &(this->VyNodes);
	NodeComposite* aPnodes = (i == 0) ? &(this->aPu) : &(this->aPv);
	nodes->InitializeEquations();
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	double mu = this->liquid.viscosity;
	auto viscosityFunction = [mu](const Vector3D& P) {return mu; };
	this->AddDiffusion_cellBased(*nodes, viscosityFunction);
	this->AddConvection_cellBased(*nodes);
	//Adding the pressure source term:
	Vector3D ez(0, 0, 1);
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			double b = 0;
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			while (e)
			{
				Vector3D M = e->GetMidPoint();
				Vector3D CM = M - C;
				Vector3D eta = e->GetVector();
				if (CM.CrossProduct2D(eta) < 0)
				{
					eta = eta * (-1);
				}
				double p = this->PNodes.GetNode(e)->value;
				double etaL = eta.abs();
				Vector3D e_eta = eta / etaL;
				Vector3D n = e_eta && ez;
				n(2) = 0;//avoiding numerical issues.
				b += -n(i) * p * etaL;
				e = ite.Next();
			}
			Node* aPnode = aPnodes->GetNode(f);
			aPnode->value = nodes->GetK(f, f);
			nodes->AddToC(f, -b);
		}
		f = itf.Next();
	}
	this->ConvergenceV[i] = nodes->SolveAndUpdate(0.5);
	//nodes->Populate();
}
void FVM::SolveSIMPLE_p()
{
	//Based on Versteeg & Malalasekara equation (11.90)
	NodeComposite* nodes = &(this->PNodes);
	nodes->InitializeEquations();
	Vector3D ez(0, 0, 1);
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	Face* f0 = 0;
	while (f)
	{
		if (f != fb)
		{
			if (!f0)
				f0 = f;
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			double b = 0;
			double aPx = this->aPu.GetValue(f);
			double aPy = this->aPv.GetValue(f);
			double uxP = this->VxNodes.GetValue(f);
			double uyP = this->VyNodes.GetValue(f);
			Vector3D uP(uxP, uyP);
			double areaP = this->triangulation->GetArea(f);
			Vector3D gradPress_P = this->GetGradient(this->PNodes, f);
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			while (e)
			{
				Vector3D M = e->GetMidPoint();
				Vector3D CM = M - C;
				Vector3D eta = e->GetVector();
				if (CM.CrossProduct2D(eta) < 0)
				{
					eta = eta * (-1);
				}
				double etaL = eta.abs();
				Vector3D e_eta = eta / etaL;
				Face* f1 = e->GetOtherFace(f);
				Vector3D C1;
				GeoGraphObject* colPtr;
				double areaA = 0;
				Vector3D gradPress_A;
				if (f1 != fb)
				{
					C1 = f1->GetCenteriod();
					colPtr = f1;
					areaA = this->triangulation->GetArea(f1);
					gradPress_A = this->GetGradient(this->PNodes, f1);
				}
				else
				{
					C1 = (e->GetOrig()->GetPoint() + e->GetDest()->GetPoint()) / 2.0;
					colPtr = e;
				}
				Vector3D ksi = C1 - C;
				double ksiL = ksi.abs();
				Vector3D e_ksi = ksi / ksiL;
				Vector3D n = e_eta && ez;
				n(2) = 0;//avoiding numerical issues.
				double uxA = this->VxNodes.GetValue(colPtr);
				double uyA = this->VyNodes.GetValue(colPtr);
				Vector3D uA(uxA, uyA);
				double aAx = this->aPu.GetValue(colPtr);
				double aAy = this->aPv.GetValue(colPtr);
				double un_P = uP || n;
				double un_A = uA || n;
				b += -0.5 * (un_P + un_A) * etaL;
				double x1 = 0.5 * areaP;
				double denom = (aPx * n(0) + aPy * n(1));
				double x2 = x1 / denom;
				double x3 = x2 * (gradPress_P || e_ksi);
				double x4 = x3 * etaL;
				b += 0.5 * areaP * (n(0) / aPx + n(1) / aPy) * (gradPress_P || e_ksi) * etaL;
				b += (colPtr == f1) ? 0.5 * areaA * (n(0) / aAx + n(1) / aAy) * (gradPress_A || e_ksi) * etaL : 0;
				double k = 0.5 * areaP * (n(0) / aPx + n(1) / aPy) / ksiL * etaL;
				k += (colPtr == f1) ? 0.5 * areaA * (n(0) / aAx + n(1) / aAy) / ksiL * etaL : 0;
				nodes->AddToK(f, colPtr, k);
				aP += k;
				e = ite.Next();
			}
			nodes->AddToK(f, f, -aP);
			nodes->AddToC(f, -b);
		}
		f = itf.Next();
	}
	this->PNodes.SetValueInEquations(f0, 0);
	this->ConvergenceP = this->PNodes.SolveAndAdd(0.5);
	nodes->Populate();
}
void FVM::SolveSIMPLE_CorrectV()
{
	//Based on Versteeg & Malalasekara equation (11.90)
	Vector3D ez(0, 0, 1);
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		if (!this->VxNodes.IsConstValue(e));
			this->VxNodes.AddToValue(e, 0);
		if (!this->VyNodes.IsConstValue(e));
			this->VyNodes.AddToValue(e, 0);
		e = ite.Next();
	}
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			double b = 0;
			double aPx = this->aPu.GetValue(f);
			double aPy = this->aPv.GetValue(f);
			double uxP = this->VxNodes.GetValue(f);
			double uyP = this->VyNodes.GetValue(f);
			Vector3D uP(uxP, uyP);
			double areaP = this->triangulation->GetArea(f);
			Vector3D gradPress_P = this->GetGradient(this->PNodes, f);
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			while (e)
			{
				Vector3D M = e->GetMidPoint();
				Vector3D CM = M - C;
				Vector3D eta = e->GetVector();
				if (CM.CrossProduct2D(eta) < 0)
				{
					eta = eta * (-1);
				}
				double etaL = eta.abs();
				Vector3D e_eta = eta / etaL;
				Vector3D C1;
				Vector3D ksi = C1 - C;
				double ksiL = ksi.abs();
				Vector3D e_ksi = ksi / ksiL;
				Vector3D n = e_eta && ez;
				n(2) = 0;//avoiding numerical issues.
				double un_P = uP || n;
				double uf = 0.5 * un_P;
				uf += 0.5 * areaP * (n(0) / aPx + n(1) / aPy) / ksiL * this->PNodes.GetValue(f);
				uf -= 0.5 * areaP * (n(0) / aPx + n(1) / aPy) * (gradPress_P || e_ksi);
				if (!this->VxNodes.IsConstValue(e));
					this->VxNodes.AddToValue(e, uf * n(0));
				if (!this->VyNodes.IsConstValue(e));
					this->VyNodes.AddToValue(e, uf * n(1));
				e = ite.Next();
			}
		}
		f = itf.Next();
	}
}
void FVM::populateVelocities(function<double(const Vector3D& P)> V[2])
{
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	if (!this->VxNodes.IsInitialized())
		this->VxNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::VERTICES);
	if (!this->VyNodes.IsInitialized())
		this->VyNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::VERTICES);
	this->VxNodes.SetValueAllNodes(V[0], false);
	this->VyNodes.SetValueAllNodes(V[1], false);
}
void FVM::populateVertexVelocities()
{
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	VertexIterator itv(qe);
	Vertex* v = itv.Next();
	while (v)
	{
		bool isAtBoundary = false;
		double weightSum = 0;
		double vx = 0;
		double vy = 0;
		Vector3D Pv = v->GetPoint();
		EdgeIterator ite(v);
		Edge* e = ite.Next();
		while (e)
		{
			if (e->GetRight() == fb || e->GetLeft() == fb)
			{
				isAtBoundary = true;
				break;
			}
			Vector3D Pe = e->GetPoint();
			double dist = Pv.distance(Pe);
			double weight = 1.0 / dist;
			vx += this->VxNodes.GetValue(e) * weight;
			vy += this->VyNodes.GetValue(e) * weight;
			weightSum += weight;
			e = ite.Next();
		}
		if (!isAtBoundary)
		{
			vx /= weightSum;
			vy /= weightSum;
			this->VxNodes.SetValue(v, vx);
			this->VyNodes.SetValue(v, vy);
		}
		v = itv.Next();
	}
}
void FVM::Debugger()
{
	//---- DEBUGGING ------------
	Face* fb = this->triangulation->GetBoundary();
	QuadEdge* qe = this->triangulation->GetMesh2D();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	int counter = 1;
	while (f)
	{
		if (f != fb)
		{
			Vector3D P = f->GetPoint();
			double Px = P(0);
			double Py = P(1);
			Node* n = this->VxNodes.GetNode(f);
			int i = n->index;
			double vx = this->VxNodes.GetValue(f);
			double vy = this->VyNodes.GetValue(f);
			double p = this->PNodes.GetValue(f);
			++counter;
		}
		f = itf.Next();
		
	}
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		Vector3D P = e->GetPoint();
		double Px = P(0);
		double Py = P(1);
		Node* n = this->VxNodes.GetNode(e);
		int i = n->index;
		double vx = this->VxNodes.GetValue(e);
		double vy = this->VyNodes.GetValue(e);
		double p = this->PNodes.GetValue(e);
		e = ite.Next();
		++counter;
	}
	VertexIterator itv(qe);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D P = v->GetPoint();
		double Px = P(0);
		double Py = P(1);
		Node* n = this->VxNodes.GetNode(v);
		int i = n->index;
		double vx = this->VxNodes.GetValue(v);
		double vy = this->VyNodes.GetValue(v);
		double p = this->PNodes.GetValue(v);
		v = itv.Next();
		++counter;
	}
	//---------------------------
}
FVM::FVM(const Triangulation& T, CELL_GEOMETRY cell_geo_) : cell_geo(cell_geo_)
{
	this->triangulation = new Triangulation(T);
}
FVM::~FVM()
{
	if (this->triangulation)
		delete this->triangulation;
	this->triangulation = 0;
}
void FVM::AddDiffusion(double conductivity)
{
	this->InitializeTNodes();
	if (this->cell_geo == CELL_GEOMETRY::VERTEX_BASED)
		this->AddDiffusion_vertexBased(conductivity);
	else
		this->constantConductivity = conductivity;
}
void FVM::AddDiffusionAndConvection(function<double(const Vector3D& P)> conductivity, double density, function<double(const Vector3D& P)> V[2])
{
	this->InitializeTNodes();
	if (this->cell_geo == CELL_GEOMETRY::VERTEX_BASED)
		this->AddDiffusionAndConvection_vertexBased(conductivity, density, V);
	else
	{
		this->conductivityFunction = conductivity;
		this->constantDensity = density;
		this->vFunction[0] = V[0];
		this->vFunction[1] = V[1];
	}
}
void FVM::AddDiffusionAndConvection(double conductivity, double density, function<double(const Vector3D& P)> V[2])
{
	auto constant_conductivity = [conductivity](const Vector3D& P) {return conductivity; };
	this->AddDiffusionAndConvection(constant_conductivity, density, V);
}
void FVM::SetThermalBoundaryConditions(function<double(const Vector3D& P)> BC)
{
	this->mode = MODE::CONSERVATION;
	this->InitializeTNodes();
	Face* fb = this->triangulation->GetBoundary();
	VertexIterator itv(fb);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D P = v->GetPoint();
		double T = BC(P);
		this->TNodes.SetConstantValue(v, T);
		v = itv.Next();
	}
	if (this->cell_geo == CELL_GEOMETRY::CELL_BASED)
	{
		EdgeIterator ite(fb);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetMidPoint();
			double T = BC(P);
			this->TNodes.SetConstantValue(e, T);
			e = ite.Next();
		}
	}
}
void FVM::SetVxBoundaryConditions(function<double(const Vector3D& P)> BC)
{
	Face* fb = this->triangulation->GetBoundary();
	VertexIterator itv(fb);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D P = v->GetPoint();
		double vx = BC(P);
		this->VxNodes.SetConstantValue(v, vx);
		v = itv.Next();
	}
	if (this->cell_geo == CELL_GEOMETRY::CELL_BASED)
	{
		EdgeIterator ite(fb);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetMidPoint();
			double vx = BC(P);
			this->VxNodes.SetConstantValue(e, vx);
			e = ite.Next();
		}
	}
}
void FVM::SetVyBoundaryConditions(function<double(const Vector3D& P)> BC)
{
	Face* fb = this->triangulation->GetBoundary();
	VertexIterator itv(fb);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D P = v->GetPoint();
		double vy = BC(P);
		this->VyNodes.SetConstantValue(v, vy);
		v = itv.Next();
	}
	if (this->cell_geo == CELL_GEOMETRY::CELL_BASED)
	{
		EdgeIterator ite(fb);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetMidPoint();
			double vy = BC(P);
			this->VyNodes.SetConstantValue(e, vy);
			e = ite.Next();
		}
	}
}
void FVM::SetFlowProblem(const LiquidProperties& liq)
{
	this->mode = MODE::SIMPLE;
	this->liquid = liq;
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	this->VxNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->VyNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->PNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS);
	this->aPu.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->aPv.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
}
void FVM::SetFlowProblem(double Re)
{
	LiquidProperties liq;
	liq.density = 1.0;
	liq.viscosity = 1.0 / Re;
	this->SetFlowProblem(liq);
}
void FVM::Solve(vector<Node*>& results)
{
	if (this->mode == MODE::CONSERVATION)
	{
		if (this->cell_geo == CELL_GEOMETRY::VERTEX_BASED)
		{
			this->TNodes.SolveAndUpdate();
			this->TNodes.Populate();
		}
		else
		{
			int iter = -1;
			int max_iter = 50;
			double Min_convergence = 0.001;
			while (iter < max_iter)
			{
				++iter;
				this->TNodes.InitializeEquations();
				double k = this->constantConductivity;
				auto conductivityFunction = [k](const Vector3D& P) {return k; };
				this->AddDiffusion_cellBased(this->TNodes, conductivityFunction);
				if (this->vFunction[0] || this->vFunction[1])
				{
					this->populateVelocities(this->vFunction);
					this->AddConvection_cellBased(this->TNodes);
				}
				double convergence = this->TNodes.SolveAndUpdate();
				if (convergence < Min_convergence)
					break;
				this->TNodes.Populate();
			}
		}
		this->TNodes.GetResults_Vertices(results);
	}
}
double FVM::Get_TValue(const Vector& P)
{
	return this->TNodes.GetValue(P);
}
double FVM::GetTValue(GeoGraphObject* objPtr, bool& isValid)
{
	return this->TNodes.GetValue(objPtr, isValid);
}
double FVM::GetVxValue(GeoGraphObject* objPtr, bool& isValid)
{
	return this->VxNodes.GetValue(objPtr, isValid);
}
double FVM::GetVyValue(GeoGraphObject* objPtr, bool& isValid)
{
	return this->VyNodes.GetValue(objPtr, isValid);
}
double FVM::GetPValue(GeoGraphObject* objPtr, bool& isValid)
{
	return this->PNodes.GetValue(objPtr, isValid);
}
bool tester_FVM(int& NumTests)
{
	if (!tester_FVM_1(NumTests))
		return false;
	if (!tester_FVM_2(NumTests))
		return false;
	if (!tester_FVM_3(NumTests))
		return false;
	if (!tester_FVM_4(NumTests))
		return false;
	if (!tester_FVM_5(NumTests))
		return false;
	if (!tester_FVM_6(NumTests))
		return false;
	if (!tester_FVM_7(NumTests))
		return false;
	if (!tester_FVM_8(NumTests))
		return false;
	if (!tester_FVM_9(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_1(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(2, 2, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	fvm.AddDiffusion(k);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.Solve(results);
	double T4_actual = cond.T(1, 1, 0);
	double T5 = cond.T(1, 2, 0);
	double T4_expected = (T5 + 3.0) / 4.0;
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		if (i != 4 && i != 5 && fabs(T - 1) > 1.0e-10)
			return false;
		else if (i == 5 && fabs(T - T5) > 1.0e-10)
			return false;
		else if (i == 4 && fabs(T - T4_expected) > 1.0e-10)
			return false;
		Vector3D P = results[i]->GetPoint();
		double T_actual = cond.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error) 
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 25)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_2(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	fvm.AddDiffusion(k);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.Solve(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double T_actual = cond.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 1.6)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_3(int& NumTests)
{//conduction-convection Rectangular 2 x 2
	double Lx = 10, Ly = 10;
	double k = 10.0;
	double rho = 1.0;
	double vx = 1.5, vy = 2.5, T0 = 1.0, T1 = 1;
	double kapa = k / rho;
	convectionConstVelocity conv(T0, T1, vx, vy, Lx, Ly, kapa);
	Triangulation T = Triangulation::OffDiagonalGrid(20, 20, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k,rho,V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.Solve(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double T_actual = conv.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 0.5)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_4(int& NumTests)
{//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1985 Example 1 Case 1
	double Lx = sqrt(2.0), Ly = sqrt(2.0);
	double k = 1;
	double rho = 1.0;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	BaligaPatankarExample1 conv;
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.Solve(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	double error_percent_CenterNode = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double T_actual = conv.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
		if (i == 60)
			error_percent_CenterNode = error_percent;
	}
	if (max_error_percent > 0.2)
		return false;
	if (error_percent_CenterNode > 0.05)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_5(int& NumTests)
{//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1985 Example 1 Case 2
	double Lx = sqrt(2.0), Ly = sqrt(2.0);
	double rho = 1.0;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	BaligaPatankarExample1_case2 conv;
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	auto variableConductivity = [conv](const Vector3D& P) {return conv.variableConductivity(P); };
	fvm.AddDiffusionAndConvection(variableConductivity, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.Solve(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	double error_percent_CenterNode = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double T_actual = conv.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
		if (i == 60)
			error_percent_CenterNode = error_percent;
	}
	if (max_error_percent > 0.5)
		return false;
	if (error_percent_CenterNode > 0.02)
		return false; 
	++NumTests;
	return true;
}
bool tester_FVM_6(int& NumTests)
{
	if (!tester_FVM_6a(NumTests))
		return false;
	if (!tester_FVM_6b(NumTests))
		return false;
	if (!tester_FVM_6c(NumTests))
		return false;
	if (!tester_FVM_6d(NumTests))
		return false;
	if (!tester_FVM_6e(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6a(int& NumTests, vector<double>* verticalCenterlineResults)
{
	//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1985 Example 2
	BaligaPatankarExample2 conv;
	conv.alpha = 0;
	conv.averageOnBoundary = false;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double rho = 1;
	double k = 0;
	double T00 = conv.T(0, 0, 0);
	vector<Node*> results;
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	fvm.Solve(results);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	ofstream file;
	file.open("..\\Data\\tester_FVM_6.1.txt", fstream::out);
	file << "Node#, x, y, Vx, Vy, T_actual, T, error, erro_percent";
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = conv.T(x, y);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) * 100;
		file << endl << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
		if (fabs(x - conv.Lx / 2.0) < 1.0e-10)
		{
			if (abs_error > max_abs_error_VerticalCenterline)
				max_abs_error_VerticalCenterline = abs_error;
			if(verticalCenterlineResults)
				verticalCenterlineResults->push_back(T);
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 1.0)
		return false;
	if (max_error_percent > 25)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6b(int& NumTests, vector<double>* verticalCenterlineResults)
{
	BaligaPatankarExample2 conv; 
	conv.alpha = 0.3;
	conv.averageOnBoundary = false;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node*> results;
	FVM fvm2(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm2.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm2.SetThermalBoundaryConditions(temperatureField);
	fvm2.Solve(results);
	ofstream file; 
	file.open("..\\Data\\tester_FVM_6.2.txt", fstream::out);
	file << "Node#, x, y, Vx, Vy, T_actual, T, error, erro_percent";
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = conv.T(x, y);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) * 100;
		file << endl << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
		if (fabs(x - conv.Lx / 2.0) < 1.0e-10)
		{
			if (abs_error > max_abs_error_VerticalCenterline)
				max_abs_error_VerticalCenterline = abs_error;
			if (verticalCenterlineResults)
				verticalCenterlineResults->push_back(T);
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.24)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6c(int& NumTests, vector<double>* verticalCenterlineResults)
{
	BaligaPatankarExample2 conv; 
	conv.alpha = 0.5;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node*> results;
	FVM fvm3(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm3.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm3.SetThermalBoundaryConditions(temperatureField);
	fvm3.Solve(results);
	ofstream file;
	file.open("..\\Data\\tester_FVM_6.3.txt", fstream::out);
	file << "Node#, x, y, Vx, Vy, T_actual, T, error, erro_percent";
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = conv.T(x, y);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) * 100;
		file << endl << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
		if (fabs(x - conv.Lx / 2.0) < 1.0e-10)
		{
			if (abs_error > max_abs_error_VerticalCenterline)
				max_abs_error_VerticalCenterline = abs_error;
			if (verticalCenterlineResults)
				verticalCenterlineResults->push_back(T);
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.0001)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6d(int& NumTests, vector<double>* verticalCenterlineResults)
{
	BaligaPatankarExample2 conv; 
	conv.alpha = 0.8;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node*> results;
	FVM fvm4(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm4.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm4.SetThermalBoundaryConditions(temperatureField);
	fvm4.Solve(results);
	ofstream file;
	file.open("..\\Data\\tester_FVM_6.4.txt", fstream::out);
	file << "Node#, x, y, Vx, Vy, T_actual, T, error, erro_percent";
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = conv.T(x, y);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) * 100;
		file << endl << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
		if (fabs(x - conv.Lx / 2.0) < 1.0e-10)
		{
			if (abs_error > max_abs_error_VerticalCenterline)
				max_abs_error_VerticalCenterline = abs_error;
			if (verticalCenterlineResults)
				verticalCenterlineResults->push_back(T);
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.23)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6e(int& NumTests, vector<double>* verticalCenterlineResults)
{
	BaligaPatankarExample2 conv; 
	conv.alpha = 1.0;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node*> results;
	FVM fvm5(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm5.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm5.SetThermalBoundaryConditions(temperatureField);
	fvm5.Solve(results);
	ofstream file;
	file.open("..\\Data\\tester_FVM_6.5.txt", fstream::out);
	file << "Node#, x, y, Vx, Vy, T_actual, T, error, erro_percent";
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = conv.T(x, y);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) * 100;
		file << endl << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
		if (fabs(x - conv.Lx / 2.0) < 1.0e-10)
		{
			if (abs_error > max_abs_error_VerticalCenterline)
				max_abs_error_VerticalCenterline = abs_error;
			if (verticalCenterlineResults)
				verticalCenterlineResults->push_back(T);
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.3)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_7(int& NumTests)
{
	double Lx = 6, Ly = 6, T0 = 1.0, T1 = 10 / (exp(pi) - exp(-pi));
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(2, 2, Lx, Ly);
	FVM fvm(T);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	fvm.AddDiffusion(k);
	vector<Node*> results;
	fvm.Solve(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double T_actual = cond.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 2)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_8(int& NumTests)
{
	double Lx = 6, Ly = 6, T0 = 1.0, T1 = 10 / (exp(pi) - exp(-pi));
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	fvm.AddDiffusion(k);
	vector<Node*> results;
	fvm.Solve(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double T_actual = cond.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 2)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_9(int& NumTests)
{
	//conduction-convection Rectangular 2 x 2
	double Lx = 10, Ly = 10;
	double k = 10.0;
	double rho = 1.0;
	double vx = 1.5, vy = 2.5, T0 = 1.0, T1 = 1;
	double kapa = k / rho;
	convectionConstVelocity conv(T0, T1, vx, vy, Lx, Ly, kapa);
	Triangulation T = Triangulation::OffDiagonalGrid(6, 6, Lx, Ly);
	FVM fvm(T);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.Solve(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double T_actual = conv.T(P);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 0.1)
		return false;
	++NumTests;
	return true;
}
