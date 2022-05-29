#include <iostream>
#include "FVM_Grid.h"
#include "FVM.h"//for tester calsses
#include "MathUtils.h"
#include "Vector.h"
#include "Vector3D.h"
#include "SquareMatrix.h"
#include "QuadEdge.h"
#include "FaceIterator.h"
#include "VertexIterator.h"
#include "EdgeIterator.h"
#include "GeoGraphObject.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include "ShapeFunction.h"
double FVM_Grid::alpha_p = 0.3;
double FVM_Grid::alpha_v = 0.3;
double FVM_Grid::ConvergenceTolVx = 1.0e-6;
double FVM_Grid::ConvergenceTolVy = 1.0e-6;
double FVM_Grid::ConvergenceTolP = 1.0e-6;
bool FVM_Grid::change_alpha = false;
double FVM_Grid::Lx = 0;
double FVM_Grid::Ly = 0;
void BoundaryValues::copyBody(const BoundaryValues& rhs)
{
	for (int j = 0; j < 6; ++j) {
		this->boundaryValue[j] = NULL;
		if (rhs.boundaryValue[j])
			this->boundaryValue[j] = new double(*rhs.boundaryValue[j]);
	}
}
BoundaryValues::BoundaryValues()
{
	for (int j = 0; j < 6; ++j) {
		this->boundaryValue[j] = NULL;
	}
}
BoundaryValues::BoundaryValues(double v)
{
	for (int j = 0; j < 6; ++j) {
		this->boundaryValue[j] = new double(v);
	}
}
BoundaryValues::BoundaryValues(const BoundaryValues& rhs)
{
	this->copyBody(rhs);
}
BoundaryValues::~BoundaryValues()
{
	for (int j = 0; j < 6; ++j) {
		if (this->boundaryValue[j])
			delete this->boundaryValue[j];
		this->boundaryValue[j] = NULL;
	}
}
void BoundaryValues::operator=(const BoundaryValues& rhs)
{
	this->copyBody(rhs);
}
void BoundaryValues::Set(GRID_DIRECTION side, double value)
{
	if (!this->boundaryValue[side])
		this->boundaryValue[side] = new double;
	*this->boundaryValue[side] = value;
}
void FVM_Grid::AddDiffusionEquations()
{
	double k = this->thermalConductivity;
	Vector3D ez(0, 0, 1);
	int N = this->grid->NumCells();
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	FaceIterator itf(qe);
	Face* f = itf.Next();
	while (f)
	{
		if (f != fb)
		{
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			while (e)
			{
				double L = e->GetVector().abs();
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
				Vector3D delta = C1 - C;
				double deltaL = delta.abs();
				double Diff = k * L / deltaL;
				this->TNodes.AddToK(f, colPtr, Diff);
				aP += Diff;
				e = ite.Next();
			}
			this->TNodes.AddToK(f, f, -aP);
		}
		f = itf.Next();
	}
}
void FVM_Grid::AddDiffusionConvectionEquations()
{
	//Based on page 99 of Patankar's book
	double k = this->thermalConductivity;
	Vector3D ez(0, 0, 1);
	int N = this->grid->NumCells();
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
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
				double dL = e->GetVector().abs();
				Face* f1 = e->GetOtherFace(f);
				Vector3D C1;
				GeoGraphObject* gPtr;
				if (f1 != fb)
				{
					C1 = f1->GetCenteriod();
					gPtr = f1;
				}
				else
				{
					C1 = (e->GetOrig()->GetPoint() + e->GetDest()->GetPoint()) / 2.0;
					gPtr = e;
				}
				Vector3D delta = C1 - C;
				double deltaL = delta.abs();
				Vector3D edgeCenter = (e->GetOrig()->GetPoint() + e->GetDest()->GetPoint()) / 2.0;
				Vector3D et = e->GetVector() / dL;
				Vector3D en = ez && et;
				en(2) = 0;
				if ((delta || en) < 0)
					en = en * (-1);
				bool Horizontal = this->grid->IsHorizontal(e);
				double vn = Horizontal ? this->VyNodes.GetValue(e) : this->VxNodes.GetValue(e);
				if ((Horizontal && en(1) < 0) || (!Horizontal && en(0) < 0))
					vn *= -1;
				double Diff = k * dL / deltaL;
				double Flux = density * vn * dL;
				double Pe = fabs(Flux / Diff);
				double A0 = pow(1.0 - 0.1 * Pe, 5);
				double A = A0 > 0 ? A0 : 0;
				double a = Diff * A + (Flux < 0 ? -Flux : 0);//Flux < 0 means an incoming Flux
				this->TNodes.AddToK(f, gPtr, a);
				aP += a;
				e = ite.Next();
			}
			this->TNodes.AddToK(f, f, -aP);
			double area = this->grid->GetCellArea();
			if (this->Source_TVx)
			{
				double vx = this->VxNodes.GetNode(f)->value;
				b += this->Source_TVx * vx * area;
			}
			if (this->Source_TVy)
			{
				double vy = this->VyNodes.GetNode(f)->value;
				b += this->Source_TVy * vy * area;
			}
			this->TNodes.AddToC(f, -b);
		}
		f = itf.Next();
	}
}
void FVM_Grid::SolveSIMPLE_vx()
{
	this->VxNodes.InitializeEquations();
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e) {
		if (!this->grid->IsHorizontal(e) && e->GetRight() != fb && e->GetLeft() != fb)
		{
			Face* fR = e->GetRight();
			Face* fL = e->GetLeft();
			Vertex* vO = e->GetOrig();
			Vertex* vD = e->GetDest();
			Node* n[4];//east, west, north, south
			n[EAST] = this->VxNodes.GetNode(fR);
			n[WEST] = this->VxNodes.GetNode(fL);
			n[NORTH] = this->VyNodes.GetNode(vD);
			n[SOUTH] = this->VyNodes.GetNode(vO);
			GeoGraphObject* eNB[4];//east, west, north, south
			eNB[EAST] = e->GetNext(fR)->GetNext(fR);
			eNB[WEST] = e->GetNext(fL)->GetNext(fL);
			Edge* e_north = e->GetNext(vD);
			while (this->grid->IsHorizontal(e_north)) {
				e_north = e_north->GetNext(vD);
			}
			eNB[NORTH] = e_north;
			if (e_north == e)
				eNB[NORTH] = vD;
			Edge* e_south = e->GetNext(vO);
			while (this->grid->IsHorizontal(e_south)) {
				e_south = e_south->GetNext(vO);
			}
			eNB[SOUTH] = e_south;
			if (e_south == e)
				eNB[SOUTH] = vO;
			double vNormal[4];//east, west, north, south
			vNormal[EAST] = n[EAST]->value;
			vNormal[WEST] = -n[WEST]->value;
			vNormal[NORTH] = n[NORTH]->value;
			vNormal[SOUTH] = -n[SOUTH]->value;
			double aP = 0;
			double b = 0;
			double a[4];
			double mu = this->liquid.viscosity;
			double rho = this->liquid.density;
			double dL[4];
			double dy = dL[EAST] = dL[WEST] = e->GetVector().abs();//dy
			Vector3D dX = n[EAST]->GetPoint() - n[WEST]->GetPoint();
			dL[NORTH] = dL[SOUTH] = dX.abs();
			Vector3D C = e->GetMidPoint();
			for (int side = 0; side < 4; ++side)
			{
				a[side] = 0;
				Vector3D C1 = eNB[side]->GetPoint();
				Vector3D delta = C1 - C;
				double deltaL = delta.abs();
				double Diff = mu * dL[side] / deltaL;
				double Flux = rho * vNormal[side] * dL[side];
				double Pe = fabs(Flux / Diff);
				double A0 = pow(1.0 - 0.1 * Pe, 5);
				double A = A0 > 0 ? A0 : 0;
				a[side] = Diff * A + (Flux < 0 ? -Flux : 0);//Flux < 0 means an incoming Flux 
				this->VxNodes.AddToK(e, eNB[side], a[side]);
				aP += a[side];
			}
			this->VxNodes.AddToK(e, e, -aP);
			double area = this->grid->GetCellArea();
			if (this->Source_VxT)
			{
				double T = this->TNodes.GetNode(e)->value;
				b += this->Source_VxT * T * area;
			}
			double Pe = this->PNodes.GetValue(fR);
			double Pw = this->PNodes.GetValue(fL);
			double rhs = (Pe - Pw) * dy - b;
			this->VxNodes.AddToC(e, rhs);
			this->dValue[e] = dy / aP;
		}
		e = ite.Next();
	}
	this->ConvergedVx = this->VxNodes.SolveAndUpdate(FVM_Grid::alpha_v);
}
void FVM_Grid::SolveSIMPLE_vy()
{
	this->VyNodes.InitializeEquations(); 
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e) {
		if (this->grid->IsHorizontal(e) && e->GetRight() != fb && e->GetLeft() != fb)
		{
			Face* fR = e->GetRight();
			Face* fL = e->GetLeft();
			Vertex* vO = e->GetOrig();
			Vertex* vD = e->GetDest();
			Node* n[4];//east, west, north, south
			n[EAST] = this->VxNodes.GetNode(vD);
			n[WEST] = this->VxNodes.GetNode(vO);
			n[NORTH] = this->VyNodes.GetNode(fL);
			n[SOUTH] = this->VyNodes.GetNode(fR);
			GeoGraphObject* eNB[4];//east, west, north, south
			eNB[NORTH] = e->GetNext(fL)->GetNext(fL);
			eNB[SOUTH] = e->GetNext(fR)->GetNext(fR);
			Edge* e_East = e->GetNext(vD);
			while (!this->grid->IsHorizontal(e_East)) {
				e_East = e_East->GetNext(vD);
			}
			eNB[EAST] = e_East;
			if (e_East == e)
				eNB[EAST] = vD;
			Edge* e_West = e->GetNext(vO);
			while (!this->grid->IsHorizontal(e_West)) {
				e_West = e_West->GetNext(vO);
			}
			eNB[WEST] = e_West;
			if (e_West == e)
				eNB[WEST] = vO;
			double vNormal[4];//east, west, north, south
			vNormal[EAST] = n[EAST]->value;
			vNormal[WEST] = -n[WEST]->value;
			vNormal[NORTH] = n[NORTH]->value;
			vNormal[SOUTH] = -n[SOUTH]->value;
			double aP = 0;
			double b = 0;
			double a[4];
			double mu = this->liquid.viscosity;
			double rho = this->liquid.density;
			double dL[4];
			Vector3D dY = n[NORTH]->GetPoint() - n[SOUTH]->GetPoint();
			dL[EAST] = dL[WEST] = dY.abs();
			double dx = dL[NORTH] = dL[SOUTH] = e->GetVector().abs();//dx
			Vector3D C = e->GetMidPoint();
			for (int side = 0; side < 4; ++side) {
				a[side] = 0;
				Vector3D C1 = eNB[side]->GetPoint();
				Vector3D delta = C1 - C;
				double deltaL = delta.abs();
				double Diff = mu * dL[side] / deltaL;
				double Flux = rho * vNormal[side] * dL[side];
				double Pe = fabs(Flux / Diff);
				double A0 = pow(1.0 - 0.1 * Pe, 5);
				double A = A0 > 0 ? A0 : 0;
				a[side] = Diff * A + (Flux < 0 ? -Flux : 0);//Flux < 0 means an incoming Flux
				this->VyNodes.AddToK(e, eNB[side], a[side]);
				aP += a[side];
			}
			this->VyNodes.AddToK(e, e, -aP);
			double area = this->grid->GetCellArea();
			if (this->Source_VyT)
			{
				double T = this->TNodes.GetNode(e)->value;
				b += this->Source_VyT * T * area;
			}
			double Pn = this->PNodes.GetValue(fL);
			double Ps = this->PNodes.GetValue(fR);
			double rhs = (Pn - Ps) * dx - b;
			this->VyNodes.AddToC(e, rhs);
			this->dValue[e] = dx / aP;
		}
		e = ite.Next();
	}
	this->ConvergedVy = this->VyNodes.SolveAndUpdate(FVM_Grid::alpha_v);
}
void FVM_Grid::SolveSIMPLE_p()
{
	this->PNodes.InitializeEquations(); 
	QuadEdge* qe = this->grid->GetMesh2D();
	Vector3D ez(0, 0, 1);
	FaceIterator itf(qe);
	Face* fb = this->grid->GetBoundary();
	Face* f = itf.Next();
	Face* f0 = f;
	while (f) {
		if (f != fb) {
			Vector3D C = f->GetCenteriod();
			double aP = 0;
			double b = 0;
			EdgeIterator ite(f);
			Edge* e = ite.Next();
			while (e) {
				double dL = e->GetVector().abs();
				Face* f1 = e->GetOtherFace(f);
				Vector3D C1 = f1->GetCenteriod();
				Vector3D delta = C1 - C;
				Vector3D et = e->GetVector() / dL;
				Vector3D en = ez && et;
				en(2) = 0;
				if ((delta || en) < 0)
					en = en * (-1);
				bool Horizontal = this->grid->IsHorizontal(e);
				double vn = Horizontal ? this->VyNodes.GetValue(e) : this->VxNodes.GetValue(e);
				if ((Horizontal && en(1) < 0) || (!Horizontal && en(0) < 0))
					vn *= -1;
				double rho = this->liquid.density;
				if (f1 != fb)
				{
					double d = this->dValue[e];
					double a = rho * d * dL;
					aP += a;
					this->PNodes.AddToK(f, f1, a);
				}
				b += -rho * vn * dL;
				e = ite.Next();
			}
			this->PNodes.AddToK(f, f, -aP);
			this->PNodes.AddToC(f, -b);
		}
		f = itf.Next();
	}
	this->PNodes.SetValueInEquations(f0, 0);
	this->ConvergedP = this->PNodes.SolveAndAdd(FVM_Grid::alpha_p);
}
void FVM_Grid::SolveSIMPLE_CorrectV()
{
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e) {
		if (e->GetRight() != fb && e->GetLeft() != fb)
		{
			bool Horizontal = this->grid->IsHorizontal(e); 
			Node* nV = Horizontal ? this->VyNodes.GetNode(e) : this->VxNodes.GetNode(e);
			double P_Right = this->PNodes.GetSolution(e->GetRight());
			double P_Left = this->PNodes.GetSolution(e->GetLeft());
			double d = this->dValue[e];
			double correction = 1.0;//Correction to the sign of terms in (6.17-19) in Patankar's book.
			nV->value += correction * FVM_Grid::alpha_p * d * (Horizontal ? (P_Right - P_Left) : (P_Left - P_Right));
		}
		e = ite.Next();
	}
	this->VxNodes.ApplyBC();
	this->VyNodes.ApplyBC();
}
void FVM_Grid::SetThermalBoundaryConditions_internal(fxy BC, GRID_DIRECTION side, string type)
{
	Face* fb = this->grid->GetBoundary();
	EdgeIterator ite(fb);
	Edge* e = ite.Next();
	bool H = this->grid->IsHorizontal(e);
	while (e)
	{
		if (side == GRID_DIRECTION::ANY_BOUNDARY || this->grid->IsBoundary(e, side))
		{
			Vector3D P = (e->GetOrig()->GetPoint() + e->GetDest()->GetPoint()) / 2.0;
			double T = (*BC)(P(0), P(1));
			if(type == "Drichlit")
				this->TNodes.SetConstantValue(e, T);
			else if(type == "Neumann")
				this->TNodes.SetConstantNormalGradient(e, T);
		}
		e = ite.Next();
	}
}
double FVM_Grid::RBC_initial_u(double x, double y)
{
	return x * (Lx - x) * y * (Ly - y) * sin(pi * y/Ly - 0.5);
}
double FVM_Grid::RBC_initial_v(double x, double y)
{
	return x * (Lx - x) * y * (Ly - y) * sin(pi * x/Lx - 0.5);
}
double FVM_Grid::RBC_initial_T(double x, double y)
{
	return x * (Lx - x) * y * (Ly - y) * sin(0.5 - pi * y / Ly);
}
FVM_Grid::FVM_Grid(const bGrid& G)
{
	this->grid = new bGrid(G);
	this->thermalConductivity = 0;
	this->mode = MODE::NONE;
	this->Source_VxT = 0;
	this->Source_VyT = 0;
	this->Source_TVx = 0;
	this->Source_VyT = 0;
	FVM_Grid::Lx = G.GetLx();
	FVM_Grid::Ly = G.GetLy();
}
FVM_Grid::~FVM_Grid()
{
	if (this->grid)
		delete this->grid;
	this->grid = 0;
}
void FVM_Grid::SetThermalConductionProblem(double conductivity)
{
	this->mode = MODE::CONDUCTION;
	this->thermalConductivity = conductivity;
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	this->TNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->TNodes.InitializeEquations();
}
void FVM_Grid::SetThermalConductionConvectionProblem(double conductivity, double Density)
{
	//Assumes velocities are already set up
	this->mode = MODE::CONVECTION;
	this->thermalConductivity = conductivity;
	this->density = Density;
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	this->TNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->TNodes.InitializeEquations();
}
void FVM_Grid::SetThermalConductionConvectionProblem(double conductivity, double Density, velocityField Velocity)
{
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary(); 
	this->VxNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::EDGES);
	this->VyNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::EDGES);
	EdgeIterator ite(qe);
	Edge* e = ite.Next();
	while (e)
	{
		this->VxNodes.SetValue(e, Velocity[0]);
		this->VyNodes.SetValue(e, Velocity[1]);
		e = ite.Next();
	}
	this->SetThermalConductionConvectionProblem(conductivity, Density);
}
void FVM_Grid::SetFlowProblem(const LiquidProperties& liq)
{
	this->mode = MODE::FLOW; 
	this->liquid = liq;
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	unordered_map<Edge*, bool> isHorizontal, isVertical;
	this->grid->GetIsHorizontal(isHorizontal);
	this->grid->GetIsVertical(isVertical);
	this->VxNodes.SetEdgeCondition(qe, &isVertical);
	this->VyNodes.SetEdgeCondition(qe, &isHorizontal);
	this->VxNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::EDGES_AND_BOUNDARY);
	this->VyNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::EDGES_AND_BOUNDARY);
	this->PNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS);
	this->dValue.clear();
	this->dValue.rehash(qe->NumEdges());
}
void FVM_Grid::SetFlowProblem(double Re)
{
	LiquidProperties liq;
	liq.density = 1.0;
	liq.viscosity = 1.0 / Re;
	this->SetFlowProblem(liq);

}
void FVM_Grid::SetFlowForcedConvectionProblem(const LiquidProperties& liq)
{
	this->SetFlowProblem(liq);
	this->mode = MODE::FLOW_CONVECTION_FORCED;
}
void FVM_Grid::SetFlowForcedConvectionProblem(double Re, double Pe)
{
	LiquidProperties liq;
	liq.density = 1.0;
	liq.viscosity = 1.0 / Re;
	liq.heatCapacity = 1.0;
	liq.thermalConductionCoeff = 1.0 / Pe;
	this->SetFlowForcedConvectionProblem(liq);
}
void FVM_Grid::SetFlowNaturalConvectionProblem(double Ra, double Pr)
{
	//Applies Boussinesq approximation 
	this->mode = MODE::FLOW_CONVECTION_NATURAL;
	this->Source_VyT = Ra;
	this->Source_TVy = 1.0;
	LiquidProperties liq;
	this->liquid.density = 1.0 / Pr;
	this->liquid.viscosity = 1.0;
	this->liquid.heatCapacity = Pr;
	this->liquid.thermalConductionCoeff = 1.0;
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	unordered_map<Edge*, bool> isHorizontal, isVertical;
	this->grid->GetIsHorizontal(isHorizontal);
	this->grid->GetIsVertical(isVertical);
	this->VxNodes.SetEdgeCondition(qe, &isVertical);
	this->VyNodes.SetEdgeCondition(qe, &isHorizontal);
	this->VxNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::EDGES_AND_BOUNDARY);
	this->VyNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::EDGES_AND_BOUNDARY);
	this->PNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS);
	this->dValue.clear();
	this->dValue.rehash(qe->NumEdges());
	this->TNodes.Initialize(qe, fb, NODE_COMPOSITE_TYPE::CELLS_AND_BOUNDARY);
	this->TNodes.InitializeEquations();
	this->VxNodes.SetValueAllNodes(FVM_Grid::RBC_initial_u);
	this->VyNodes.SetValueAllNodes(FVM_Grid::RBC_initial_v);
	this->TNodes.SetValueAllNodes(FVM_Grid::RBC_initial_T);
}
void FVM_Grid::SetThermalBoundaryConditions(fxy BC, GRID_DIRECTION side)
{
	this->SetThermalBoundaryConditions_internal(BC, side, "Drichlit");
}
void FVM_Grid::SetThermalGradientBoundaryConditions(fxy BC, GRID_DIRECTION side)
{
	this->SetThermalBoundaryConditions_internal(BC, side, "Neumann");
}
void FVM_Grid::SetThermalBoundaryConditions(const BoundaryValues& BT)
{
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	EdgeIterator ite(fb);
	Edge* e = ite.Next();
	while (e) {
		bool H = this->grid->IsHorizontal(e);
		Face* fR = e->GetRight();
		Face* fL = e->GetLeft();
		Node* n = this->TNodes.GetNode(e);
		GRID_DIRECTION side = this->grid->WhichBoundary(e);
		if (side != GRID_DIRECTION::NONE && BT.boundaryValue[side])
		{
			this->TNodes.SetConstantValue(e, *BT.boundaryValue[side]);
		}
		e = ite.Next();
	}
}
void FVM_Grid::SetThermalGradientBoundaryConditions(const BoundaryValues& BT)
{
	QuadEdge* qe = this->grid->GetMesh2D();
	Face* fb = this->grid->GetBoundary();
	EdgeIterator ite(fb);
	Edge* e = ite.Next();
	while (e) {
		bool H = this->grid->IsHorizontal(e);
		Face* fR = e->GetRight();
		Face* fL = e->GetLeft();
		Node* n = this->TNodes.GetNode(e);
		GRID_DIRECTION side = this->grid->WhichBoundary(e);
		if (side != GRID_DIRECTION::NONE && BT.boundaryValue[side])
		{
			this->TNodes.SetConstantNormalGradient(e, *BT.boundaryValue[side]);
		}
		e = ite.Next();
	}
}
void FVM_Grid::SetVelocityBoundaryConditions(const BoundaryValues& Vx, const BoundaryValues& Vy)
{
	EdgeIterator ite(this->grid->GetBoundary());
	Edge* e = ite.Next();
	while (e) 
	{
		GRID_DIRECTION side = this->grid->WhichBoundary(e);
		if (side != GRID_DIRECTION::NONE)
		{
			if (this->grid->IsHorizontal(e) && Vy.boundaryValue[side])
				this->VyNodes.SetConstantValue(e, *Vy.boundaryValue[side]);
			if (!this->grid->IsHorizontal(e) && Vx.boundaryValue[side])
				this->VxNodes.SetConstantValue(e, *Vx.boundaryValue[side]);
		}
		e = ite.Next();
	}
	VertexIterator itv(this->grid->GetBoundary());
	Vertex* v = itv.Next();
	while (v) 
	{
		GRID_DIRECTION side = this->grid->WhichBoundary(v);
		if (side != GRID_DIRECTION::NONE)
		{
			if (Vx.boundaryValue[side])
				this->VxNodes.SetConstantValue(v, *Vx.boundaryValue[side]);
			if (Vy.boundaryValue[side])
				this->VyNodes.SetConstantValue(v, *Vy.boundaryValue[side]);
		}
		v = itv.Next();
	}
}
void FVM_Grid::SetVelocityGradientBoundaryConditions(const BoundaryValues& Vx, const BoundaryValues& Vy)
{
	VertexIterator itv(this->grid->GetBoundary());
	Vertex* v = itv.Next();
	while (v)
	{
		GRID_DIRECTION side = this->grid->WhichBoundary(v);
		if (side != GRID_DIRECTION::NONE)
		{
			if (Vx.boundaryValue[side])
				this->VxNodes.SetConstantNormalGradient(v, *Vx.boundaryValue[side]);
			if (Vy.boundaryValue[side])
				this->VyNodes.SetConstantNormalGradient(v, *Vy.boundaryValue[side]);
		}
		v = itv.Next();
	}
}
void FVM_Grid::SolveThermalProblem(vector<Node*>& cell_results)
{
	if (this->mode == MODE::CONDUCTION)
		this->AddDiffusionEquations();
	else
		this->AddDiffusionConvectionEquations();
	this->TNodes.SolveAndUpdate();
	this->TNodes.Populate();
	this->TNodes.GetResults_Cells(cell_results);
}
void FVM_Grid::SolveSIMPLE(vector<Node*>& Vx, vector<Node*>& Vy, vector<Node*>& P, int maxIter)
{
	string fileName = "..\\Data\\FVM_Grid.output.txt";
	ofstream fout;
	fout.open(fileName.c_str(), fstream::out);
	fout.close();
	int iter = -1;
	vector<double> ConvergedVxArr;
	vector<double> ConvergedVyArr;
	vector<double> ConvergedPArr;
	double minConvergedVx = 10000;
	double minConvergedVy = 10000;
	double minConvergedP = 10000;
	if (FVM_Grid::change_alpha)
	{
		ConvergedVxArr.resize(maxIter, 0);
		ConvergedVyArr.resize(maxIter, 0);
		ConvergedPArr.resize(maxIter, 0);
	}
	while (iter < maxIter)
	{
		++iter;
		this->SolveSIMPLE_vx();
		this->SolveSIMPLE_vy();
		this->SolveSIMPLE_p();
		this->SolveSIMPLE_CorrectV();
		this->VxNodes.Populate();
		this->VyNodes.Populate();
		this->PNodes.Populate();
		if (this->mode == MODE::FLOW_CONVECTION_NATURAL)
		{
			this->AddDiffusionConvectionEquations();
			this->TNodes.SolveAndUpdate();
			this->TNodes.Populate();
		}
		this->printEquations(fileName, iter);
		if (this->ConvergedVx < FVM_Grid::ConvergenceTolVx &&
			this->ConvergedVy < FVM_Grid::ConvergenceTolVy &&
			this->ConvergedP < FVM_Grid::ConvergenceTolP)
			break;
		if (this->ConvergedVx > 1e10 || this->ConvergedVy > 1e10 || this->ConvergedP > 1e10)
			break;
		if (FVM_Grid::change_alpha)
		{
			ConvergedVxArr[iter] = this->ConvergedVx;
			ConvergedVyArr[iter] = this->ConvergedVy;
			ConvergedPArr[iter] = this->ConvergedP;
			if (iter > 5)
			{
				minConvergedVx = minConvergedVx < ConvergedVxArr[iter - 4] ? minConvergedVx : ConvergedVxArr[iter - 4];
				minConvergedVy = minConvergedVy < ConvergedVyArr[iter - 4] ? minConvergedVy : ConvergedVyArr[iter - 4];
				minConvergedP = minConvergedP < ConvergedPArr[iter - 4] ? minConvergedP : ConvergedPArr[iter - 4];
				if (this->ConvergedVx > minConvergedVx || this->ConvergedVy > minConvergedVy)
				{
					if (FVM_Grid::alpha_v > 0.05)
						FVM_Grid::alpha_v *= 0.9;
					minConvergedVx *= 1.2;
					minConvergedVy *= 1.2;

				}
				if (this->ConvergedP > minConvergedP)
				{
					if (FVM_Grid::alpha_p > 0.05)
						FVM_Grid::alpha_p *= 0.9;
					minConvergedP *= 1.2;
				}
			}
		}
	}
	this->VxNodes.GetResults_Vertices(Vx);
	this->VyNodes.GetResults_Vertices(Vy);
	this->PNodes.GetResults_Vertices(P);
}
void FVM_Grid::SolveSIMPLE(vector<Node*>& Vx, vector<Node*>& Vy, vector<Node*>& P, vector<Node*>& T, int maxIter)
{
	this->SolveSIMPLE(Vx, Vy, P, maxIter);
	if (this->mode == MODE::FLOW_CONVECTION_FORCED)
	{
		this->SetThermalConductionConvectionProblem(this->liquid.thermalConductionCoeff,
			this->liquid.density * this->liquid.heatCapacity);
		this->TNodes.SolveAndUpdate();
		this->TNodes.Populate();
	}
	this->TNodes.GetResults_Vertices(T);
}
void FVM_Grid::printEquations(string fileName, int iterNum)
{
	/*ofstream fout;
	fout.open(fileName.c_str(), fstream::app);
	fout << endl << "===========================================================";
	fout << endl << "Iteration #" << iterNum;
	this->VxNodes.PrintEquations(fout, "Vx");
	this->VyNodes.PrintEquations(fout, "Vy");
	this->PNodes.PrintEquations(fout, "P");
	fout.close();*/
	cout << endl << "Iter #" << iterNum;
	cout << ", Convergence(Vx,Vy,P) = (" << this->ConvergedVx << ", ";
	cout << this->ConvergedVy << ", " << this->ConvergedP << ")";
	cout << endl << "		(aV,aP) = (" << FVM_Grid::alpha_v << "," << FVM_Grid::alpha_p << ")";
}
void FVM_Grid::printThermalEquations(string fileName)
{
	ofstream fout;
	fout.open(fileName.c_str(), fstream::app);
	fout << endl << "===========================================================";
	fout << endl << "Thermal Equations";
	this->TNodes.PrintEquations(fout, "T");
	fout.close();
}
bool tester_FVM_Grid(int& NumTests)
{
	if (!tester_FVM_Grid_1(NumTests))
		return false;
	if (!tester_FVM_Grid_2(NumTests))
		return false; 
	if (!tester_FVM_Grid_3(NumTests))
		return false;
	if (!tester_FVM_Grid_7(NumTests))
		return false;
	if (!tester_FVM_Grid_9(NumTests))
		return false;
	if (!tester_FVM_Grid_10(NumTests))
		return false;
	if (!tester_FVM_Grid_11(NumTests))
		return false;
	if (!tester_FVM_Grid_12(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_Grid_1(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	testerFVM1::Set(Lx, Ly, T0, T1);
	bGrid G(1, 1, Lx, Ly);
	FVM_Grid fvm(G);
	fvm.SetThermalConductionProblem(k);
	fvm.SetThermalBoundaryConditions(&testerFVM1::value);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = testerFVM1::value(x, y);
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
bool tester_FVM_Grid_2(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	testerFVM1::Set(Lx, Ly, T0, T1);
	bGrid G(10, 10, Lx, Ly);
	FVM_Grid fvm(G);
	fvm.SetThermalConductionProblem(k);
	fvm.SetThermalBoundaryConditions(&testerFVM1::value);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = testerFVM1::value(x, y);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 1.1)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_Grid_3(int& NumTests)
{
	//conduction-convection Rectangular 2 x 2
	int N = 10;
	double Lx = 10, Ly = 10;
	double k = 10.0;
	double rho = 1.0;
	double vx = 1.5, vy = 2.5, T0 = 1.0, T1 = 1;
	double kapa = k / rho;
	testerFVM3::Set(T0, T1, vx, vy, Lx, Ly, kapa);
	velocityField V;
	V[0] = testerFVM3::Vx;
	V[1] = testerFVM3::Vy;
	bGrid G(N, N, Lx, Ly);
	G.GetMesh2D()->print("..\\Data\\tester_FVM_Grid_3.Grid.txt");
	FVM_Grid fvm(G);
	fvm.SetThermalConductionConvectionProblem(k, rho, V);
	fvm.SetThermalBoundaryConditions(&testerFVM3::value);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		double x = results[i]->GetPoint()(0);
		double y = results[i]->GetPoint()(1);
		double T_actual = testerFVM3::value(x, y);
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (N >= 10 && max_error_percent > 0.04)
		return false;
	if (N >= 2 && max_error_percent > 0.5)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_Grid_7(int& NumTests)
{
	testerFVM7::Lx = 1.0;
	testerFVM7::Ly = 1.0;
	testerFVM7::density = 1;
	testerFVM7::viscosity = 1.0;
	testerFVM7::VLid = 1;
	testerFVM7::Nx = 6;
	testerFVM7::Ny = 6;
	int maxIter = 1000;
	BoundaryValues VBx(0), VBy(0);
	*VBx.boundaryValue[NORTH] = testerFVM7::VLid;
	LiquidProperties liq;
	liq.density = testerFVM7::density;
	liq.viscosity = testerFVM7::viscosity;
	bGrid G(testerFVM7::Nx, testerFVM7::Ny, testerFVM7::Lx, testerFVM7::Ly);
	G.GetMesh2D()->print("..\\Data\\tester_FVM_Grid_7.Grid.txt");
	FVM_Grid::alpha_p = 0.3;
	FVM_Grid::alpha_v = 0.3;
	FVM_Grid::ConvergenceTolVx = 1.0e-6;
	FVM_Grid::ConvergenceTolVy = 1.0e-6;
	FVM_Grid::ConvergenceTolP = 1.0e-6;
	FVM_Grid fvm(G);
	fvm.SetFlowProblem(liq);
	fvm.SetVelocityBoundaryConditions(VBx, VBy);
	vector<Node*> Vx,Vy,P;
	fvm.SolveSIMPLE(Vx, Vy, P, maxIter);
	string fileName = "..\\Data\\FVM_Grid.output.txt";
	ofstream fout;
	fout.open(fileName.c_str(), fstream::app);
	fout << endl << "Vx at Vertices:";
	for (int i = 0; i < Vx.size(); ++i)
	{
		fout << endl << "Vx[" << i << "] = " << Vx[i]->value;
	}
	fout << endl << "Vy at Vertices:";
	for (int i = 0; i < Vy.size(); ++i)
	{
		fout << endl << "Vy[" << i << "] = " << Vy[i]->value;
	}
	fout << endl << "P at Vertices:";
	for (int i = 0; i < P.size(); ++i)
	{
		fout << endl << "P[" << i << "] = " << P[i]->value;
	}
	fout.close();
	int iPmin = 0, iPmax = 0, iUmin = 0, iUmax = 0, iVmin = 0, iVmax = 0;
	double Pmin = P[iPmin]->value, Pmax = P[iPmax]->value;
	double Umin = Vx[iUmin]->value, Umax = Vx[iUmax]->value;
	double Vmin = Vy[iVmin]->value, Vmax = Vy[iVmax]->value;
	for (int i = 0; i < P.size(); ++i) {
		double p = P[i]->value;
		double u = Vx[i]->value;
		double v = Vy[i]->value;
		if (p < Pmin) { Pmin = p; iPmin = i; }
		if (p > Pmax) { Pmax = p; iPmax = i; }
		if (u < Umin) { Umin = u; iUmin = i; }
		if (u > Umax) { Umax = u; iUmax = i; }
		if (v < Vmin) { Vmin = v; iVmin = i; }
		if (v > Vmax) { Vmax = v; iVmax = i; }
	}
	double dP_max = Pmax - Pmin;
	double dP_max_ANSYS = 22.73;//Value from ANSYS Fluent , SIMPLE, linear, 6x6 mesh
	double dP_max_error = fabs(dP_max - dP_max_ANSYS) / dP_max_ANSYS * 100;
	if (testerFVM7::Nx >= 6 && testerFVM7::Ny >= 6 && dP_max_error > 9)
		return false;
	double Umin_ANSYS = -0.1455;
	double Umin_error = fabs((Umin - Umin_ANSYS) / Umin_ANSYS) * 100.0;
	double Vmin_ANSYS = -0.2015;
	double Vmin_error = fabs((Vmin - Vmin_ANSYS) / Vmin_ANSYS) * 100.0;
	double Vmax_ANSYS = 0.2007;
	double Vmax_error = fabs((Vmax - Vmax_ANSYS) / Vmax_ANSYS) * 100.0;
	if (testerFVM7::Nx >= 6 && testerFVM7::Ny >= 6)
	{
		if (Umin_error > 14 || Vmin_error > 12 || Vmax_error > 12)
			return false;
	}
	++NumTests;
	return true;
	/*Nx x Ny Number of iterations that the code converges with alpha_p = alpha_v = 0.3, with iteration tolerance of 1.0e6 for all 3 equations
	6x6     89
	*/
}
bool tester_FVM_Grid_9(int& NumTests)
{
	testerFVM9 test9;
	test9.Re = 100;
	test9.Pe = test9.Re;
	test9.Nx = 6;
	test9.Ny = 6;
	int maxIter = 1000;
	BoundaryValues VBx(0), VBy(0), TB(0);
	*VBx.boundaryValue[NORTH] = test9.VLid;
	*TB.boundaryValue[NORTH] = test9.TLid;
	bGrid G(test9.Nx, test9.Ny, test9.Lx, test9.Ly);
	G.GetMesh2D()->print("..\\Data\\tester_FVM_Grid_7.Grid.txt");
	FVM_Grid::alpha_p = 0.3;
	FVM_Grid::alpha_v = 0.3;
	FVM_Grid::ConvergenceTolVx = 1e-3;
	FVM_Grid::ConvergenceTolVy = 1e-3;
	FVM_Grid::ConvergenceTolP = 1e-3;
	FVM_Grid fvm(G);
	fvm.SetFlowForcedConvectionProblem(test9.Re, test9.Pe);
	fvm.SetVelocityBoundaryConditions(VBx, VBy);
	fvm.SetThermalBoundaryConditions(TB);
	vector<Node*> Vx, Vy, P, T;
	fvm.SolveSIMPLE(Vx, Vy, P, T, maxIter);
	string fileName = "..\\Data\\FVM_Grid.output.txt";
	fvm.printThermalEquations("..\\Data\\FVM_Grid.output.txt");
	ofstream fout;
	fout.open(fileName.c_str(), fstream::app);
	fout << endl << "Vx at Vertices:";
	for (int i = 0; i < Vx.size(); ++i)
	{
		fout << endl << "Vx[" << i << "] = " << Vx[i]->value;
	}
	fout << endl << "Vy at Vertices:";
	for (int i = 0; i < Vy.size(); ++i)
	{
		fout << endl << "Vy[" << i << "] = " << Vy[i]->value;
	}
	fout << endl << "P at Vertices:";
	for (int i = 0; i < P.size(); ++i)
	{
		fout << endl << "P[" << i << "] = " << P[i]->value;
	}
	fout << endl << "T at Vertices:";
	for (int i = 0; i < T.size(); ++i)
	{
		fout << endl << "T[" << i << "] = " << T[i]->value;
	}
	fout.close();
	string fileNameThermal = "..\\Data\\FVM_Grid.T.csv";
	fout.open(fileNameThermal.c_str(), fstream::out);
	fout << "x,y,T";
	for (int i = 0; i < T.size(); ++i)
	{
		fout << endl << T[i]->GetPoint()(0) << "," << T[i]->GetPoint()(1) << "," << T[i]->value;
	}
	fout.close();
	int iPmin = 0, iPmax = 0, iUmin = 0, iUmax = 0, iVmin = 0, iVmax = 0;
	double Pmin = P[iPmin]->value, Pmax = P[iPmax]->value;
	double Umin = Vx[iUmin]->value, Umax = Vx[iUmax]->value;
	double Vmin = Vy[iVmin]->value, Vmax = Vy[iVmax]->value;
	for (int i = 0; i < P.size(); ++i) {
		double p = P[i]->value;
		double u = Vx[i]->value;
		double v = Vy[i]->value;
		if (p < Pmin) { Pmin = p; iPmin = i; }
		if (p > Pmax) { Pmax = p; iPmax = i; }
		if (u < Umin) { Umin = u; iUmin = i; }
		if (u > Umax) { Umax = u; iUmax = i; }
		if (v < Vmin) { Vmin = v; iVmin = i; }
		if (v > Vmax) { Vmax = v; iVmax = i; }
	}
	++NumTests;
	return true;
}
bool tester_FVM_Grid_10(int& NumTests)//RBC Box
{
	double Ra = 3000;
	double Pr = 7;
	int Nx = 6;
	int Ny = 6;
	int maxIter = 1000;
	double Lx = 1.0;
	double Ly = 1.0;
	bGrid G(Nx, Ny, Lx, Ly);
	FVM_Grid::alpha_p = 0.3;
	FVM_Grid::alpha_v = 0.3;
	FVM_Grid::ConvergenceTolVx = 1e-3;
	FVM_Grid::ConvergenceTolVy = 1e-3;
	FVM_Grid::ConvergenceTolP = 1e-3;
	FVM_Grid fvm(G); 
	fvm.SetFlowNaturalConvectionProblem(Ra, Pr);
	BoundaryValues ZeroBC(0);
	fvm.SetVelocityBoundaryConditions(ZeroBC, ZeroBC);
	fvm.SetThermalBoundaryConditions(ZeroBC);
	vector<Node*> Vx, Vy, P, T;
	fvm.SolveSIMPLE(Vx, Vy, P, T, maxIter);
	ofstream fout; 
	string fileNameThermal = "..\\Data\\FVM_Grid.T.csv";
	fout.open(fileNameThermal.c_str(), fstream::out);
	fout << "x,y,T";
	for (int i = 0; i < T.size(); ++i)
	{
		fout << endl << T[i]->GetPoint()(0) << "," << T[i]->GetPoint()(1) << "," << T[i]->value;
	}
	fout.close();
	++NumTests;
	return true;
}
bool tester_FVM_Grid_11(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	int Nx = 6;
	int Ny = Nx;
	testerFVM1::Set(Lx, Ly, T0, T1);
	bGrid G(Nx, Ny, Lx, Ly);
	FVM_Grid fvm(G);
	fvm.SetThermalConductionProblem(k);
	fvm.SetThermalBoundaryConditions(&testerFVM1::value, WEST);
	fvm.SetThermalBoundaryConditions(&testerFVM1::value, SOUTH);
	fvm.SetThermalGradientBoundaryConditions(&testerFVM1::DvalueDy, NORTH);
	fvm.SetThermalGradientBoundaryConditions(&testerFVM1::DvalueDx, EAST);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
	fvm.printThermalEquations("..\\Data\\FVM_Grid.output.txt");
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		Vector3D P = results[i]->GetPoint();
		double x = P(0);
		double y = P(1);
		double T_actual = testerFVM1::value(x, y);
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
bool tester_FVM_Grid_12(int& NumTests)//test # 10 , infinitely extended RBC.
{
	double Ra = 3000;
	double Pr = 7;
	int Nx = 9;
	int Ny = 6;
	int maxIter = 1000;
	double Lx = 1.5;
	double Ly = 1.0;
	bGrid G(Nx, Ny, Lx, Ly);
	FVM_Grid::alpha_p = 0.1;
	FVM_Grid::alpha_v = 0.1;
	FVM_Grid::ConvergenceTolVx = 1e-3;
	FVM_Grid::ConvergenceTolVy = 1e-3;
	FVM_Grid::ConvergenceTolP = 1e-3;
	FVM_Grid fvm(G);
	fvm.SetFlowNaturalConvectionProblem(Ra, Pr);
	BoundaryValues nullBC;
	BoundaryValues BVx(0);
	BoundaryValues BVy, BdVy_dx;
	BVy.Set(NORTH, 0);
	BVy.Set(SOUTH, 0);
	BdVy_dx.Set(EAST, 0);
	BdVy_dx.Set(WEST, 0);
	BoundaryValues BT, dBT_dx;
	BT.Set(NORTH, 0);
	BT.Set(SOUTH, 0);
	dBT_dx.Set(EAST, 0);
	dBT_dx.Set(WEST, 0);
	fvm.SetVelocityBoundaryConditions(BVx, BVy);
	fvm.SetVelocityGradientBoundaryConditions(nullBC, BdVy_dx);
	fvm.SetThermalBoundaryConditions(BT);
	fvm.SetThermalGradientBoundaryConditions(dBT_dx);
	vector<Node*> Vx, Vy, P, T;
	fvm.SolveSIMPLE(Vx, Vy, P, T, maxIter);
	ofstream fout;
	string fileNameThermal = "..\\Data\\FVM_Grid.T.csv";
	fout.open(fileNameThermal.c_str(), fstream::out);
	fout << "x,y,T";
	for (int i = 0; i < T.size(); ++i)
	{
		fout << endl << T[i]->GetPoint()(0) << "," << T[i]->GetPoint()(1) << "," << T[i]->value;
	}
	fout.close();
	++NumTests;
	return true;
}
