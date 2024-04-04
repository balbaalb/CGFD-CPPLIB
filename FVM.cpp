#include <fstream>
#include <iostream>
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
		QuadEdge *qe = this->triangulation->GetMesh2D();
		Face *fb = this->triangulation->GetBoundary();
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
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			VertexIterator itv(f);
			Vertex *ver[3];
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
				if (AO.CrossProduct2D(AP) > 0)
				{
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
void FVM::AddDiffusion_cellBased(NodeComposite &nodes, function<double(const Vector3D &P)> conductivity)
{
	Vector3D ez(0, 0, 1);
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			double b = 0;
			EdgeIterator ite(f);
			Edge *e = ite.Next();
			while (e)
			{
				Vector3D M = e->GetMidPoint();
				Vector3D CM = M - C;
				Vector3D eta = e->GetVector();
				Vertex *Va = e->GetOrig();
				Vertex *Vb = e->GetDest();
				if (CM.CrossProduct2D(eta) < 0)
				{
					eta = eta * (-1);
					std::swap(Va, Vb);
				}
				double Phi_a = nodes.GetValue(Va);
				double Phi_b = nodes.GetValue(Vb);
				double etaL = eta.abs();
				Vector3D e_eta = eta / etaL;
				Face *f1 = e->GetOtherFace(f);
				Vector3D C1;
				GeoGraphObject *colPtr;
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
				n(2) = 0; // avoiding numerical issues.
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
void FVM::AddDiffusionAndConvection_vertexBased(function<double(const Vector3D &P)> conductivity, double density, function<double(const Vector3D &P)> V[2])
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
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			VertexIterator itv(f);
			Vertex *ver[3];
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
			ShapeFunction *N;
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
				if (AO.CrossProduct2D(AP) > 0)
				{
					nOP = nOP * (-1);
					nOQ = nOQ * (-1);
				}
				Vector3D VO(V[0](O), V[1](O));
				Vector3D VP(V[0](P), V[1](P));
				Vector3D VQ(V[0](Q), V[1](Q));
				Vector3D VR(V[0](R), V[1](R));
				Vector3D VS(V[0](S), V[1](S));
				;
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
Vector3D FVM::GetGradient(NodeComposite &nodes, Face *f)
{
	Vector3D C = f->GetCenteriod();
	bMatrix A(3, 2);
	Vector B(3);
	EdgeIterator ite(f);
	Edge *e = ite.Next();
	Face *fb = this->triangulation->GetBoundary();
	Node *n0 = nodes.GetNode(f);
	double Phi0 = n0->value;
	int row = 0;
	while (e)
	{
		Face *f1 = e->GetOtherFace(f);
		Vector3D C1;
		GeoGraphObject *colPtr;
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
Vector3D FVM::GetGradient_p(Face *f)
{
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	vector<Face *> f1_vec;
	{
		EdgeIterator ite(f);
		Edge *e = ite.Next();
		while (e)
		{
			Face *f1 = e->GetOtherFace(f);
			if (f1 != fb)
			{
				f1_vec.push_back(f1);
			}
			e = ite.Next();
		}
	}
	if (f1_vec.size() == 3)
	{
		return this->GetGradient(this->PNodes, f);
	}
	if (f1_vec.size() == 2)
	{
		Vector3D C = f->GetCenteriod();
		bMatrix A(2, 2);
		Vector B(2);
		double Phi0 = this->PNodes.GetValue(f);
		int row = 0;
		EdgeIterator ite(f);
		Edge *e = ite.Next();
		while (e)
		{
			Face *f1 = e->GetOtherFace(f);
			if (f1 != fb)
			{
				Vector3D C1 = f1->GetCenteriod();
				A(row, 0) = C1(0) - C(0);
				A(row, 1) = C1(1) - C(1);
				double Phi = this->PNodes.GetValue(f1);
				B(row) = Phi - Phi0;
				++row;
			}
			e = ite.Next();
		}
		bMatrix ATA = A.Transpose() * A;
		double determinant = ATA(0, 0) * ATA(1, 1) - ATA(0, 1) * ATA(1, 0);
		bMatrix ATA_inv = ATA / determinant;
		std::swap(ATA_inv(1, 0), ATA_inv(0, 1));
		Vector X = ATA_inv * A.Transpose() * B;
		return Vector3D(X(0), X(1));
	}
	if (f1_vec.size() == 1)
	{
		Face *f1 = f1_vec[0];
		Vector3D C0 = f->GetCenteriod();
		Vector3D C1 = f1->GetCenteriod();
		Vector3D Ksi = C1 - C0;
		double ds = Ksi.abs();
		Vector e_Ksi = Ksi / ds;
		double Phi0 = this->PNodes.GetValue(f);
		double Phi1 = this->PNodes.GetValue(f1);
		double dPhi = Phi1 - Phi0;
		double gradP = dPhi / ds;
		return e_Ksi * gradP;
	}
	return Vector3D();
}
void FVM::AddConvection_cellBased(NodeComposite &nodes)
{
	double k = this->constantConductivity;
	double rho = this->constantDensity;
	Vector3D ez(0, 0, 1);
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	TVD psi(this->TVD_mode);
	FaceIterator itf(qe);
	Face *f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			double b = 0;
			Vector3D grad_Phi = this->GetGradient(nodes, f);
			EdgeIterator ite(f);
			Edge *e = ite.Next();
			while (e)
			{
				Vector3D M = e->GetMidPoint();
				Node *uxNode = this->VxNodes.GetNode(e);
				Node *uyNode = this->VyNodes.GetNode(e);
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
				Face *f1 = e->GetOtherFace(f);
				Vector3D C1;
				GeoGraphObject *colPtr;
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
				n(2) = 0; // avoiding numerical issues.
				double un = uM || n;
				double flux = rho * un * etaL;
				double Phi0 = nodes.GetNode(f)->value;
				double Phi1 = nodes.GetNode(colPtr)->value;
				double Phi_downstream{0}, Phi_upstream{0};
				if (flux > 0)
				{
					Phi_downstream = Phi1;
					Phi_upstream = Phi0;
				}
				else
				{
					nodes.AddToK(f, colPtr, -flux);
					aP += -flux;
					Phi_downstream = Phi0;
					Phi_upstream = Phi1;
				}
				double dPhi = Phi_downstream - Phi_upstream;
				double r = fabs(dPhi) > 1.0e-10 ? 2 * (grad_Phi || ksi) / dPhi - 1 : 0;
				b += -flux * psi(r) / 2.0 * dPhi;
				e = ite.Next();
			}
			nodes.AddToK(f, f, -aP);
			nodes.AddToC(f, -b);
		}
		f = itf.Next();
	}
}
void FVM::SolveSIMPLE_v()
{
	// NodeComposite* nodes = (i == 0) ? &(this->VxNodes) : &(this->VyNodes);
	this->VxNodes.InitializeEquations();
	this->VyNodes.InitializeEquations();
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	double mu = this->liquid.viscosity;
	auto viscosityFunction = [mu](const Vector3D &P)
	{ return mu; };
	this->AddDiffusion_cellBased(this->VxNodes, viscosityFunction);
	this->AddConvection_cellBased(this->VxNodes);
	this->AddDiffusion_cellBased(this->VyNodes, viscosityFunction);
	this->AddConvection_cellBased(this->VyNodes);
	// Adding the pressure source term:
	Vector3D ez(0, 0, 1);
	FaceIterator itf(qe);
	Face *f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			Vector3D b(-this->gradPx.GetValue(f), -this->gradPy.GetValue(f));
			b = b * this->triangulation->GetArea(f);
			this->VxNodes.AddToC(f, -b(0));
			this->VyNodes.AddToC(f, -b(1));
		}
		f = itf.Next();
	}
	this->VxNodes.StabilizeK();
	this->VyNodes.StabilizeK();
	this->ConvergenceV[0] = this->VxNodes.SolveAndUpdate(this->alpha_v);
	this->ConvergenceV[1] = this->VyNodes.SolveAndUpdate(this->alpha_v);
}
void FVM::SolveSIMPLE_updateHelper()
{
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			double aPu = fabs(this->VxNodes.GetK(f, f));
			this->Helper.SetValue(f, aPu);
		}
		f = itf.Next();
	}
	EdgeIterator ite(qe);
	Edge *e = ite.Next();
	while (e)
	{
		Vector3D n = e->GetNormal_2D();
		Face *fR = e->GetRight();
		Face *fL = e->GetLeft();
		if (fR != fb && fL != fb) // This should be changed if a boundary with a known pressure but unknown velocity is considered.
		{
			Vector3D vR, vL;
			vR(0) = this->VxNodes.GetValue(fR);
			vL(0) = this->VxNodes.GetValue(fL);
			vR(1) = this->VyNodes.GetValue(fR);
			vL(1) = this->VyNodes.GetValue(fL);
			double uf = 0.5 * (vR || n) + 0.5 * (vL || n);
			Vector3D gradP_P(this->gradPx.GetValue(fL), this->gradPy.GetValue(fL));
			Vector3D gradP_A(this->gradPx.GetValue(fR), this->gradPy.GetValue(fR));
			Vector3D P = fL->GetPoint();
			Vector3D A = fR->GetPoint();
			Vector3D Xi = A - P;
			Vector3D e_xi = Xi / Xi.abs();
			double aP = this->Helper.GetValue(fL);
			double aA = this->Helper.GetValue(fR);
			double areaP = this->triangulation->GetArea(fL);
			double areaA = this->triangulation->GetArea(fR);
			uf -= 0.5 * areaP / aP * (gradP_P || e_xi);
			uf -= 0.5 * areaA / aA * (gradP_A || e_xi);
			this->Helper.SetValue(e, uf);
		}
		e = ite.Next();
	}
}
void FVM::SolveSIMPLE_p()
{
	// Based on Versteeg & Malalasekara equation (11.90)
	NodeComposite *nodes = &(this->PNodes);
	nodes->InitializeEquations();
	Vector3D ez(0, 0, 1);
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	Face *f0 = 0;
	while (f)
	{
		if (f != fb)
		{
			if (!f0)
				f0 = f;
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			double b = 0;
			double aPu = this->Helper.GetValue(f);
			double areaP = this->triangulation->GetArea(f);
			EdgeIterator ite(f);
			Edge *e = ite.Next();
			while (e)
			{
				Face *f1 = e->GetOtherFace(f);
				if (f1 != fb)
				{
					double uf = this->Helper.GetValue(e);
					Vector3D M = e->GetMidPoint();
					Vector3D CM = M - C;
					Vector3D eta = e->GetVector();
					double etaL = eta.abs();
					Vector3D n = e->GetNormal_2D();
					if ((n || CM) < 0)
					{
						uf *= -1;
					}
					b += uf * etaL;
					Vector3D C1 = f1->GetCenteriod();
					double areaA = this->triangulation->GetArea(f1);
					Vector3D ksi = C1 - C;
					double ksiL = ksi.abs();
					double aAu = this->Helper.GetValue(f1);
					double aA = 0.5 * areaP / aPu / ksiL * etaL;
					aA += 0.5 * areaA / aAu / ksiL * etaL;
					nodes->AddToK(f, f1, -aA);
					aP += aA;
				}
				e = ite.Next();
			}
			nodes->AddToK(f, f, aP);
			nodes->AddToC(f, -b);
		}
		f = itf.Next();
	}
	this->PNodes.SetValueInEquations(f0, 0);
	this->PNodes.StabilizeK();
	this->ConvergenceP = this->PNodes.SolveAndAdd(this->alpha_p);
	this->PNodes.Populate();
}
void FVM::SolveSIMPLE_gradP()
{
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D gradP = this->GetGradient_p(f);
			this->gradPx.SetValue(f, gradP(0));
			this->gradPy.SetValue(f, gradP(1));
		}
		f = itf.Next();
	}
}
void FVM::SolveSIMPLE_CorrectV()
{
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	EdgeIterator ite(qe);
	Edge *e = ite.Next();
	while (e)
	{
		Vector3D n = e->GetNormal_2D();
		Face *fR = e->GetRight();
		Face *fL = e->GetLeft();
		if (fR != fb && fL != fb) // This should be changed if a boundary with a known pressure but unknown velocity is considered.
		{
			Vector3D vR, vL;
			vR(0) = this->VxNodes.GetValue(fR);
			vL(0) = this->VxNodes.GetValue(fL);
			vR(1) = this->VyNodes.GetValue(fR);
			vL(1) = this->VyNodes.GetValue(fL);
			double uf = 0.5 * (vR || n) + 0.5 * (vL || n);
			Vector3D gradP_P(this->gradPx.GetValue(fL), this->gradPy.GetValue(fL));
			Vector3D gradP_A(this->gradPx.GetValue(fR), this->gradPy.GetValue(fR));
			Vector3D P = fL->GetPoint();
			Vector3D A = fR->GetPoint();
			Vector3D Xi = A - P;
			Vector3D e_xi = Xi / Xi.abs();
			double aP = this->Helper.GetValue(fL);
			double aA = this->Helper.GetValue(fR);
			double areaP = this->triangulation->GetArea(fL);
			double areaA = this->triangulation->GetArea(fR);
			uf -= 0.5 * areaP / aP * (gradP_P || e_xi);
			uf -= 0.5 * areaA / aA * (gradP_A || e_xi);
			double pP = this->PNodes.GetValue(fL);
			double pA = this->PNodes.GetValue(fR);
			uf += 0.5 * (areaP / aP + areaA / aA) * (pP - pA) / Xi.abs();
			this->VxNodes.SetValue(e, uf * n(0));
			this->VyNodes.SetValue(e, uf * n(1));
			this->Helper.SetValue(e, uf);
		}
		e = ite.Next();
	}
}
void FVM::MakePressuresPositive()
{
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	double pMin = 0;
	{
		FaceIterator itf(qe);
		Face *f = itf.Next();
		while (f)
		{
			if (f != fb)
			{
				double p = this->PNodes.GetValue(f);
				if (p < pMin)
					pMin = p;
			}
			f = itf.Next();
		}
	}
	{
		FaceIterator itf(qe);
		Face *f = itf.Next();
		while (f)
		{
			if (f != fb)
			{
				double p = this->PNodes.GetValue(f);
				this->PNodes.SetValue(f, p - pMin);
			}
			f = itf.Next();
		}
	}
}
void FVM::populateVelocities(function<double(const Vector3D &P)> V[2])
{
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	if (!this->VxNodes.IsInitialized())
		this->VxNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::VERTICES);
	if (!this->VyNodes.IsInitialized())
		this->VyNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::VERTICES);
	this->VxNodes.SetValueAllNodes(V[0], false);
	this->VyNodes.SetValueAllNodes(V[1], false);
}
void FVM::populateVertexVelocities()
{
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	VertexIterator itv(qe);
	Vertex *v = itv.Next();
	while (v)
	{
		bool isAtBoundary = false;
		double weightSum = 0;
		double vx = 0;
		double vy = 0;
		Vector3D Pv = v->GetPoint();
		EdgeIterator ite(v);
		Edge *e = ite.Next();
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
	Face *fb = this->triangulation->GetBoundary();
	QuadEdge *qe = this->triangulation->GetMesh2D();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	int counter = 1;
	while (f)
	{
		if (f != fb)
		{
			Vector3D P = f->GetPoint();
			double Px = P(0);
			double Py = P(1);
			Node *n = this->VxNodes.GetNode(f);
			int i = n->index;
			double vx = this->VxNodes.GetValue(f);
			double vy = this->VyNodes.GetValue(f);
			double p = this->PNodes.GetValue(f);
			double dpdx = this->gradPx.GetValue(f);
			double dpdy = this->gradPy.GetValue(f);
			++counter;
		}
		f = itf.Next();
	}
	EdgeIterator ite(qe);
	Edge *e = ite.Next();
	while (e)
	{
		Vector3D P = e->GetPoint();
		double Px = P(0);
		double Py = P(1);
		Node *n = this->VxNodes.GetNode(e);
		int i = n->index;
		double vx = this->VxNodes.GetValue(e);
		double vy = this->VyNodes.GetValue(e);
		double p = this->PNodes.GetValue(e);
		double uf = this->Helper.GetValue(e);
		e = ite.Next();
		++counter;
	}
	VertexIterator itv(qe);
	Vertex *v = itv.Next();
	while (v)
	{
		Vector3D P = v->GetPoint();
		double Px = P(0);
		double Py = P(1);
		Node *n = this->VxNodes.GetNode(v);
		int i = n->index;
		double vx = this->VxNodes.GetValue(v);
		double vy = this->VyNodes.GetValue(v);
		double p = this->PNodes.GetValue(v);
		v = itv.Next();
		++counter;
	}
	//---------------------------
}
FVM::FVM(const Triangulation &T, CELL_GEOMETRY cell_geo_) : cell_geo(cell_geo_)
{
	this->triangulation = make_shared<Triangulation>(T);
}
FVM::~FVM()
{
}
void FVM::AddDiffusion(double conductivity)
{
	this->InitializeTNodes();
	if (this->cell_geo == CELL_GEOMETRY::VERTEX_BASED)
		this->AddDiffusion_vertexBased(conductivity);
	else
		this->constantConductivity = conductivity;
}
void FVM::AddDiffusionAndConvection(function<double(const Vector3D &P)> conductivity, double density, function<double(const Vector3D &P)> V[2])
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
void FVM::AddDiffusionAndConvection(double conductivity, double density, function<double(const Vector3D &P)> V[2])
{
	auto constant_conductivity = [conductivity](const Vector3D &P)
	{ return conductivity; };
	this->AddDiffusionAndConvection(constant_conductivity, density, V);
}
void FVM::SetThermalBoundaryConditions(function<double(const Vector3D &P)> BC)
{
	this->mode = MODE::CONSERVATION;
	this->InitializeTNodes();
	Face *fb = this->triangulation->GetBoundary();
	VertexIterator itv(fb);
	Vertex *v = itv.Next();
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
		Edge *e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetMidPoint();
			double T = BC(P);
			this->TNodes.SetConstantValue(e, T);
			e = ite.Next();
		}
	}
}
void FVM::SetVxBoundaryConditions(function<double(const Vector3D &P)> BC)
{
	Face *fb = this->triangulation->GetBoundary();
	VertexIterator itv(fb);
	Vertex *v = itv.Next();
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
		Edge *e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetMidPoint();
			double vx = BC(P);
			this->VxNodes.SetConstantValue(e, vx);
			e = ite.Next();
		}
	}
}
void FVM::SetVyBoundaryConditions(function<double(const Vector3D &P)> BC)
{
	Face *fb = this->triangulation->GetBoundary();
	VertexIterator itv(fb);
	Vertex *v = itv.Next();
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
		Edge *e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetMidPoint();
			double vy = BC(P);
			this->VyNodes.SetConstantValue(e, vy);
			e = ite.Next();
		}
	}
}
void FVM::SetTVD(TVD::MODE mode)
{
	this->TVD_mode = mode;
}
void FVM::SetFlowProblem(const LiquidProperties &liq)
{
	this->mode = MODE::SIMPLE;
	this->liquid = liq;
	QuadEdge *qe = this->triangulation->GetMesh2D();
	Face *fb = this->triangulation->GetBoundary();
	this->VxNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->VyNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->PNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS);
	this->Helper.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->gradPx.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS);
	this->gradPy.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS);
}
void FVM::SetFlowProblem(double Re)
{
	LiquidProperties liq;
	liq.density = 1.0;
	liq.viscosity = 1.0 / Re;
	this->SetFlowProblem(liq);
}
void FVM::SetMinconvergenceV(double value)
{
	this->MinConvergenceV = value;
}
void FVM::SetMinconvergenceP(double value)
{
	this->MinConvergenceP = value;
}
void FVM::SetMinconvergenceT(double value)
{
	this->MinConvergenceT = value;
}
void FVM::SetSolveMethod(NodeComposite::METHOD value, double tolerance)
{
	if (this->TNodes.IsInitialized())
		this->TNodes.SetSolveMethod(value, tolerance);
	if (this->VxNodes.IsInitialized())
		this->VxNodes.SetSolveMethod(value, tolerance);
	if (this->VyNodes.IsInitialized())
		this->VyNodes.SetSolveMethod(value, tolerance);
	if (this->PNodes.IsInitialized())
		this->PNodes.SetSolveMethod(value, tolerance);
}
void FVM::Solve(vector<Node *> &results)
{
	if (this->mode == MODE::SIMPLE)
	{
		this->constantConductivity = this->liquid.thermalConductionCoeff;
		this->constantDensity = this->liquid.density;
		int iter = -1;
		int max_iter = 50;
		while (iter < max_iter)
		{
			++iter;
			this->SolveSIMPLE_v();
			this->SolveSIMPLE_updateHelper();
			this->SolveSIMPLE_p();
			this->SolveSIMPLE_gradP();
			this->SolveSIMPLE_CorrectV();
			if (this->ConvergenceV[0] < this->MinConvergenceV && this->ConvergenceV[1] < this->MinConvergenceV && this->ConvergenceP < this->MinConvergenceP)
				break;
		}
		this->VxNodes.Populate();
		this->VyNodes.Populate();
		this->MakePressuresPositive();
		this->PNodes.Populate();
		this->PNodes.GetResults_Vertices(results);
	}
	else if (this->mode == MODE::CONSERVATION)
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
			while (iter < max_iter)
			{
				++iter;
				this->TNodes.InitializeEquations();
				double k = this->constantConductivity;
				auto conductivityFunction = [k](const Vector3D &P)
				{ return k; };
				this->AddDiffusion_cellBased(this->TNodes, conductivityFunction);
				if (this->vFunction[0] || this->vFunction[1])
				{
					this->populateVelocities(this->vFunction);
					this->AddConvection_cellBased(this->TNodes);
				}
				this->TNodes.StabilizeK();
				double convergence = this->TNodes.SolveAndUpdate();
				if (convergence < this->MinConvergenceT)
					break;
				this->TNodes.Populate();
			}
		}
		this->TNodes.GetResults_Vertices(results);
	}
}
void FVM::Solve(vector<Node *> &Vx, vector<Node *> &Vy, vector<Node *> &P)
{
	this->Solve(P);
	this->VxNodes.GetResults_Vertices(Vx);
	this->VyNodes.GetResults_Vertices(Vy);
}
QuadEdge *FVM::GetMesh2D()
{
	return this->triangulation->GetMesh2D();
}
Face *FVM::GetBoundary()
{
	if (!this->triangulation)
		return nullptr;
	return this->triangulation->GetBoundary();
}
double FVM::Get_TValue(const Vector &P)
{
	return this->TNodes.GetValue(P);
}
double FVM::GetTValue(GeoGraphObject *objPtr, bool &isValid)
{
	return this->TNodes.GetValue(objPtr, isValid);
}
double FVM::GetVxValue(GeoGraphObject *objPtr, bool &isValid)
{
	return this->VxNodes.GetValue(objPtr, isValid);
}
double FVM::GetVyValue(GeoGraphObject *objPtr, bool &isValid)
{
	return this->VyNodes.GetValue(objPtr, isValid);
}
double FVM::GetPValue(GeoGraphObject *objPtr, bool &isValid)
{
	return this->PNodes.GetValue(objPtr, isValid);
}
void FVM::SetAlphaP(double value)
{
	this->alpha_p = value;
}
void FVM::SetalphaV(double value)
{
	this->alpha_v = value;
}
bool tester_FVM(int &NumTests)
{
	cout << std::endl
		 << "tester_FVM_1 test ..." << flush;
	if (!tester_FVM_1(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_2 test ..." << flush;
	if (!tester_FVM_2(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_3 test ..." << flush;
	if (!tester_FVM_3(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_4 test ..." << flush;
	if (!tester_FVM_4(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_5 test ..." << flush;
	if (!tester_FVM_5(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_6 test ..." << flush;
	if (!tester_FVM_6(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_7 test ..." << flush;
	if (!tester_FVM_7(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_8 test ..." << flush;
	if (!tester_FVM_8(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_9 test ..." << flush;
	if (!tester_FVM_9(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_10 test ..." << flush;
	if (!tester_FVM_10(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_11 test ..." << flush;
	if (!tester_FVM_11(NumTests))
		return false;
	cout << "okay!" << endl
		 << "tester_FVM_12 test ..." << flush;
	if (!tester_FVM_12(NumTests))
		return false;
	cout << "okay!" << endl;
	++NumTests;
	return true;
}
template <class T>
ostream &operator<<(ostream &out, vector<T> &vec)
{
	out << "[";
	for (int i = 0; i < vec.size(); ++i)
	{
		out << vec[i];
		if (i < vec.size() - 1)
		{
			out << " , ";
		}
	}
	out << "]";
	return out;
}
bool tester_FVM_1(int &NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(2, 2, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	fvm.AddDiffusion(k);
	auto temperatureField = [cond](const Vector3D &P)
	{ return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node *> results;
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
bool tester_FVM_2(int &NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	fvm.AddDiffusion(k);
	auto temperatureField = [cond](const Vector3D &P)
	{ return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node *> results;
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
bool tester_FVM_3(int &NumTests)
{ // conduction-convection Rectangular 2 x 2
	double Lx = 10, Ly = 10;
	double k = 10.0;
	double rho = 1.0;
	double vx = 1.5, vy = 2.5, T0 = 1.0, T1 = 2.0;
	double kapa = k / rho;
	convectionConstVelocity conv(T0, T1, vx, vy, Lx, Ly, kapa);
	Triangulation T = Triangulation::OffDiagonalGrid(20, 20, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node *> results;
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
bool tester_FVM_4(int &NumTests)
{ // conduction-convection Rectangular 10 x 10, Baliga & Patankar 1985 Example 1 Case 1
	double Lx = sqrt(2.0), Ly = sqrt(2.0);
	double k = 1;
	double rho = 1.0;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	BaligaPatankarExample1 conv;
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node *> results;
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
bool tester_FVM_5(int &NumTests)
{ // conduction-convection Rectangular 10 x 10, Baliga & Patankar 1985 Example 1 Case 2
	double Lx = sqrt(2.0), Ly = sqrt(2.0);
	double rho = 1.0;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	BaligaPatankarExample1_case2 conv;
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	auto variableConductivity = [conv](const Vector3D &P)
	{ return conv.variableConductivity(P); };
	fvm.AddDiffusionAndConvection(variableConductivity, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node *> results;
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
bool tester_FVM_6(int &NumTests)
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
bool tester_FVM_6a(int &NumTests, vector<double> *verticalCenterlineResults)
{
	// conduction-convection Rectangular 10 x 10, Baliga & Patankar 1985 Example 2
	BaligaPatankarExample2 conv;
	conv.alpha = 0;
	conv.averageOnBoundary = false;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double rho = 1;
	double k = 0;
	double T00 = conv.T(0, 0, 0);
	vector<Node *> results;
	FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
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
		file << endl
			 << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
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
	if (max_abs_error_VerticalCenterline > 1.0)
		return false;
	if (max_error_percent > 25)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6b(int &NumTests, vector<double> *verticalCenterlineResults)
{
	BaligaPatankarExample2 conv;
	conv.alpha = 0.3;
	conv.averageOnBoundary = false;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node *> results;
	FVM fvm2(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	fvm2.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
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
		file << endl
			 << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
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
bool tester_FVM_6c(int &NumTests, vector<double> *verticalCenterlineResults)
{
	BaligaPatankarExample2 conv;
	conv.alpha = 0.5;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node *> results;
	FVM fvm3(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	fvm3.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
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
		file << endl
			 << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
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
bool tester_FVM_6d(int &NumTests, vector<double> *verticalCenterlineResults)
{
	BaligaPatankarExample2 conv;
	conv.alpha = 0.8;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node *> results;
	FVM fvm4(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	fvm4.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
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
		file << endl
			 << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
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
bool tester_FVM_6e(int &NumTests, vector<double> *verticalCenterlineResults)
{
	BaligaPatankarExample2 conv;
	conv.alpha = 1.0;
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, conv.Lx, conv.Ly);
	double max_abs_error_VerticalCenterline = 0;
	double max_error_percent = 0;
	double rho = 1;
	double k = 1.0e-7;
	vector<Node *> results;
	FVM fvm5(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
	function<double(const Vector3D &P)> V[2];
	V[0] = [conv](const Vector3D &P)
	{ return conv.vx(P); };
	V[1] = [conv](const Vector3D &P)
	{ return conv.vy(P); };
	fvm5.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D &P)
	{ return conv.T(P); };
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
		file << endl
			 << i << "," << x << "," << y << "," << conv.vx(x, y) << "," << conv.vy(x, y) << "," << conv.T(x, y) << "," << T << "," << abs_error << "," << error_percent << "%";
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
bool tester_FVM_7(int &NumTests)
{
	double Lx = 6, Ly = 6, T0 = 1.0, T1 = 10 / (exp(pi) - exp(-pi));
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(2, 2, Lx, Ly);
	FVM fvm(T);
	auto temperatureField = [cond](const Vector3D &P)
	{ return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	fvm.AddDiffusion(k);
	vector<Node *> results;
	fvm.SetSolveMethod(NodeComposite::METHOD::TRIANGULATION_CELL_BASED);
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
bool tester_FVM_8(int &NumTests)
{
	double Lx = 6, Ly = 6, T0 = 1.0, T1 = 10 / (exp(pi) - exp(-pi));
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(10, 10, Lx, Ly);
	FVM fvm(T);
	auto temperatureField = [cond](const Vector3D &P)
	{ return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	fvm.AddDiffusion(k);
	vector<Node *> results;
	fvm.SetSolveMethod(NodeComposite::METHOD::TRIANGULATION_CELL_BASED);
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
bool tester_FVM_9(int &NumTests)
{ // conduction-convection Rectangular 2 x 2
	double Lx = 1, Ly = 1;
	vector<double> k_cases{1e-5, 1e-4, 0.001, 0.01, 0.1, 1.0, 10, 100, 1000, 1e4, 1e5};
	double rho = 1.0;
	double theta = 15.8 / 180 * pi;
	double U = 1.0;
	double vx = U * cos(theta), vy = U * sin(theta), T0 = 1.0, T1 = 2.0;
	vector<int> N_cases{6, 7, 8};
	for (int c = 0; c < k_cases.size(); ++c)
	{
		double k = k_cases[c];
		double kapa = k / rho;
		convectionConstVelocity conv(T0, T1, vx, vy, Lx, Ly, kapa);
		vector<double> max_abs_error(N_cases.size(), 0);
		vector<double> max_error_percent(N_cases.size(), 0);
		for (int n = 0; n < N_cases.size(); ++n)
		{
			int N = N_cases[n];
			Triangulation T = Triangulation::OffDiagonalGrid(N, N, Lx, Ly);
			FVM fvm(T, FVM::CELL_GEOMETRY::VERTEX_BASED);
			function<double(const Vector3D &P)> V[2];
			V[0] = [conv](const Vector3D &P)
			{ return conv.vx(P); };
			V[1] = [conv](const Vector3D &P)
			{ return conv.vy(P); };
			fvm.AddDiffusionAndConvection(k, rho, V);
			auto temperatureField = [conv](const Vector3D &P)
			{ return conv.T(P); };
			fvm.SetThermalBoundaryConditions(temperatureField);
			vector<Node *> results;
			fvm.Solve(results);
			for (int i = 0; i < results.size(); ++i)
			{
				double T = results[i]->value;
				Vector3D P = results[i]->GetPoint();
				double T_actual = conv.T(P);
				double abs_error = fabs(T - T_actual);
				double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
				if (abs_error > max_abs_error[n])
					max_abs_error[n] = abs_error;
				if (error_percent > max_error_percent[n])
				{
					max_error_percent[n] = error_percent;
				}
			}
			if (n > 0 && max_abs_error[n] > max_abs_error[n - 1])
				return false;
			if (c == k_cases.size() - 1 && n == N_cases.size() - 1)
			{
				if (max_abs_error[n] > 1.0e-7 || max_error_percent[n] > 1.0e-5)
					return false;
			}
		}
		cout << endl
			 << "k = " << k << " , error % = " << max_error_percent << " , for mesh nxn, n = " << N_cases << endl;
	}
	++NumTests;
	return true;
}
bool tester_FVM_10(int &NumTests)
{ // Testing Lid Cavity results qualitatively
	//-----Geometry---------------------
	double Lx = 6.0;
	double Ly = 6.0;
	vector<Vector3D> input(5);
	input[1](0) = Lx;
	input[2](0) = Lx;
	input[2](1) = Ly;
	input[3](1) = Ly;
	input[4](0) = Lx / 2.0;
	input[4](1) = Ly / 2.0;
	Triangulation T = DelaunayLifting::Triangulate(input);
	T.Draw("D:/Bashar/MyStuff/A10/C616/tester_FVM_10.bmp"); // Temporary !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FVM fvm(T);
	fvm.SetTVD(TVD::MODE::UPSTREAM_DIFFERENCING);
	//-----Liquid properties---------------------
	LiquidProperties liq;
	liq.density = 1.0;
	liq.viscosity = 1.0;
	fvm.SetFlowProblem(liq);
	//-----Boundary conditions---------------------
	double VLid = 10;
	auto VxBC = [Ly, VLid](const Vector3D &P) -> double
	{ return P(1) < Ly - 0.00001 ? 0 : VLid; };
	fvm.SetVxBoundaryConditions(VxBC);
	auto VyBC = [](const Vector3D &P) -> double
	{ return 0; };
	fvm.SetVyBoundaryConditions(VyBC);
	//-----Solve---------------------
	fvm.SetMinconvergenceP(0.01);
	fvm.SetMinconvergenceV(0.01);
	vector<Node *> PNodes;
	fvm.SetSolveMethod(NodeComposite::METHOD::TRIANGULATION_CELL_BASED, 0.005);
	fvm.Solve(PNodes);
	//------Check-------------------
	QuadEdge *qe = fvm.GetMesh2D();
	FaceIterator itf(qe);
	Face *f = itf.Next();
	int counter = 0;
	double vx[]{0, 0, 0, 0};
	double vy[]{0, 0, 0, 0};
	double pres[]{0, 0, 0, 0};
	while (f)
	{
		if (f != fvm.GetBoundary())
		{
			Vector3D P = f->GetPoint();
			double Px = P(0);
			double Py = P(1);
			bool isValid = false;
			vx[counter] = fvm.GetVxValue(f, isValid);
			if (!isValid)
				return false;
			isValid = false;
			vy[counter] = fvm.GetVyValue(f, isValid);
			if (!isValid)
				return false;
			isValid = false;
			pres[counter] = fvm.GetPValue(f, isValid);
			if (!isValid)
				return false;
			++counter;
		}
		f = itf.Next();
	}
	double vx_excel[]{-0.735310619, -0.072300993, -0.542573676, +4.642079457};
	double vy_excel[]{-0.018760289, +0.000837618, -0.007909445, +0.022266843};
	double p_excel[]{-2.22359E-16, +2.26558847, -2.139272655, +0.006419896};

	if (fabs(vx[0] - vx_excel[0]) > 0.01)
		return false;
	if (fabs(vy[0] - vy_excel[0]) > 0.01)
		return false;
	if (fabs(vx[1] - vx_excel[1]) > 0.01)
		return false;
	if (fabs(vy[1] - vy_excel[1]) > 0.01)
		return false;
	if (fabs(pres[1] - pres[0] - p_excel[1] + p_excel[0]) > 0.01)
		return false;
	if (fabs(vx[2] - vx_excel[2]) > 0.01)
		return false;
	if (fabs(vy[2] - vy_excel[2]) > 0.01)
		return false;
	if (fabs(pres[2] - pres[0] - p_excel[2] + p_excel[0]) > 0.01)
		return false;
	if (fabs(vx[3] - vx_excel[3]) > 0.01)
		return false;
	if (fabs(vy[3] - vy_excel[3]) > 0.01)
		return false;
	if (fabs(pres[3] - pres[0] - p_excel[3] + p_excel[0]) > 0.01)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_11(int &NumTests)
{
	double Lx = 1.0;
	double Ly = 1.0;
	double density = 1;
	double viscosity = 1.0;
	double VLid = 1;
	int Nx = 6;
	int Ny = 6;
	Triangulation T = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	LiquidProperties liq;
	liq.density = density;
	liq.viscosity = viscosity;
	FVM fvm(T);
	fvm.SetFlowProblem(liq);
	auto VxBC = [Ly, VLid](const Vector3D &P) -> double
	{ return P(1) < Ly ? 0 : VLid; };
	fvm.SetVxBoundaryConditions(VxBC);
	auto VyBC = [](const Vector3D &P) -> double
	{ return 0; };
	fvm.SetVyBoundaryConditions(VyBC);
	fvm.SetAlphaP(0.1);
	fvm.SetalphaV(0.1);
	fvm.SetMinconvergenceP(0.1);
	fvm.SetMinconvergenceV(0.1);
	vector<Node *> Vx, Vy, P;
	fvm.SetSolveMethod(NodeComposite::METHOD::TRIANGULATION_CELL_BASED, 0.01);
	fvm.Solve(Vx, Vy, P);
	string fileName = "..\\Data\\FVM_Grid.output.txt";
	ofstream fout;
	fout.open(fileName.c_str(), fstream::app);
	fout << endl
		 << "Vx at Vertices:";
	for (int i = 0; i < Vx.size(); ++i)
	{
		fout << endl
			 << "Vx[" << i << "] = " << Vx[i]->value;
	}
	fout << endl
		 << "Vy at Vertices:";
	for (int i = 0; i < Vy.size(); ++i)
	{
		fout << endl
			 << "Vy[" << i << "] = " << Vy[i]->value;
	}
	fout << endl
		 << "P at Vertices:";
	for (int i = 0; i < P.size(); ++i)
	{
		fout << endl
			 << "P[" << i << "] = " << P[i]->value;
	}
	fout.close();
	int iPmin = 0, iPmax = 0, iUmin = 0, iUmax = 0, iVmin = 0, iVmax = 0;
	double Pmin = P[iPmin]->value, Pmax = P[iPmax]->value;
	double Umin = Vx[iUmin]->value, Umax = Vx[iUmax]->value;
	double Vmin = Vy[iVmin]->value, Vmax = Vy[iVmax]->value;
	for (int i = 0; i < P.size(); ++i)
	{
		double p = P[i]->value;
		double u = Vx[i]->value;
		double v = Vy[i]->value;
		if (p < Pmin)
		{
			Pmin = p;
			iPmin = i;
		}
		if (p > Pmax)
		{
			Pmax = p;
			iPmax = i;
		}
		if (u < Umin)
		{
			Umin = u;
			iUmin = i;
		}
		if (u > Umax)
		{
			Umax = u;
			iUmax = i;
		}
		if (v < Vmin)
		{
			Vmin = v;
			iVmin = i;
		}
		if (v > Vmax)
		{
			Vmax = v;
			iVmax = i;
		}
	}
	double dP_max = Pmax - Pmin;
	double dP_max_ANSYS = 22.73; // Value from ANSYS Fluent , SIMPLE, linear, 6x6 mesh
	double dP_max_error = fabs(dP_max - dP_max_ANSYS) / dP_max_ANSYS * 100;
	if (Nx >= 6 && Ny >= 6 && dP_max_error > 1)
		return false;
	double Umin_ANSYS = -0.1455;
	double Umin_error = fabs((Umin - Umin_ANSYS) / Umin_ANSYS) * 100.0;
	double Vmin_ANSYS = -0.2015;
	double Vmin_error = fabs((Vmin - Vmin_ANSYS) / Vmin_ANSYS) * 100.0;
	double Vmax_ANSYS = 0.2007;
	double Vmax_error = fabs((Vmax - Vmax_ANSYS) / Vmax_ANSYS) * 100.0;
	if (Nx >= 6 && Ny >= 6)
	{
		if (Umin_error > 11 || Vmin_error > 8 || Vmax_error > 3)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_FVM_12(int &NumTests)
{ // conduction-convection Rectangular 2 x 2
	double Lx = 1, Ly = 1;
	vector<double> k_cases{1e-5, 1e-4, 0.001, 0.01, 0.1, 1.0, 10, 100, 1000, 1e4, 1e5};
	double rho = 1.0;
	double theta = 15.8 / 180 * pi;
	double U = 1.0;
	double vx = U * cos(theta), vy = U * sin(theta), T0 = 1.0, T1 = 2.0;
	vector<int> N_cases{6, 7, 8};
	for (int c = 0; c < k_cases.size(); ++c)
	{
		double k = k_cases[c];
		double kapa = k / rho;
		convectionConstVelocity conv(T0, T1, vx, vy, Lx, Ly, kapa);
		vector<double> max_abs_error(N_cases.size(), 0);
		vector<double> max_error_percent(N_cases.size(), 0);
		for (int n = 0; n < N_cases.size(); ++n)
		{
			int N = N_cases[n];
			Triangulation T = Triangulation::OffDiagonalGrid(N, N, Lx, Ly);
			FVM fvm(T);
			fvm.SetTVD(TVD::MODE::VAN_LEER);
			function<double(const Vector3D &P)> V[2];
			V[0] = [conv](const Vector3D &P)
			{ return conv.vx(P); };
			V[1] = [conv](const Vector3D &P)
			{ return conv.vy(P); };
			fvm.AddDiffusionAndConvection(k, rho, V);
			auto temperatureField = [conv](const Vector3D &P)
			{ return conv.T(P); };
			fvm.SetThermalBoundaryConditions(temperatureField);
			vector<Node *> results;
			fvm.SetSolveMethod(NodeComposite::METHOD::TRIANGULATION_CELL_BASED);
			fvm.Solve(results);
			for (int i = 0; i < results.size(); ++i)
			{
				double T = results[i]->value;
				Vector3D P = results[i]->GetPoint();
				double T_actual = conv.T(P);
				double abs_error = fabs(T - T_actual);
				double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
				if (abs_error > max_abs_error[n])
					max_abs_error[n] = abs_error;
				if (error_percent > max_error_percent[n])
				{
					max_error_percent[n] = error_percent;
				}
			}
			// if (n > 0 && max_abs_error[n] > max_abs_error[n - 1])
			// 	return false;
			if (c == k_cases.size() - 1 && n == N_cases.size() - 1)
			{
				// if (max_abs_error[n] > 1.0e-7 || max_error_percent[n] > 1.0e-5)
				// 	return false;
			}
		}
		cout << endl
			 << "k = " << k << " , error % = " << max_error_percent << " , for mesh nxn, n = " << N_cases << endl;
	}
	++NumTests;
	return true;
}