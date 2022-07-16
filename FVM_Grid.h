#ifndef FVM_GridH
#define FVM_GridH
#include <string>
#include <fstream>
#include <unordered_map>
#include <functional>
#include "bGrid.h"
#include "Vector3D.h"
#include "Node.h"
#include "NodeComposite.h"
#include "TVD.h"
#include "LiquidProperties.h"
using namespace std;
class Face;
class Edge;
class Vertex;
class SquareMatrix;
class TVD;
typedef double(*fxy)(double x, double y);
typedef fxy velocityField[2];
class BoundaryValues
{
	void copyBody(const BoundaryValues& rhs);
public:
	double* boundaryValue[6];
	BoundaryValues();
	BoundaryValues(double v);
	BoundaryValues(const BoundaryValues& rhs);
	~BoundaryValues();
	void operator=(const BoundaryValues& rhs);
	void Set(GRID_DIRECTION side, double value);
};
class FVM_Grid
{
	enum MODE { NONE, CONDUCTION, CONVECTION, FLOW, FLOW_CONVECTION_FORCED, FLOW_CONVECTION_NATURAL} mode;
	bGrid* grid;
	NodeComposite TNodes;
	NodeComposite VxNodes;
	NodeComposite VyNodes;
	NodeComposite PNodes;
	LiquidProperties liquid;
	double thermalConductivity;
	double density;
	unordered_map<Edge*, double> dValue;
	double ConvergenceVx;
	double ConvergenceVy;
	double ConvergenceP;
	double ConvergenceT;
	double Source_VxT;//source term in Vx equation: Source_VxT * T; 
	double Source_VyT;//source term in Vy equation: Source_VyT * T;
	double Source_TVx;//source term in Energy equation: Source_TVx * Vx;
	double Source_TVy;//source term in Energy equation: Source_TVy * Vy;
	void AddDiffusionEquations();
	void AddDiffusionConvectionEquations();
	double GetTVD(double Flux, Face* f, Edge* e);
	void SolveSIMPLE_vx();
	void SolveSIMPLE_vy();
	void SolveSIMPLE_p();
	void SolveSIMPLE_CorrectV();
	void SetThermalBoundaryConditions_internal(function<double(const Vector3D& P)> BC, GRID_DIRECTION side, string type/*Drichlit or Neumann*/);
	double Lx;
	double Ly;
	TVD psi;
public:
	static double alpha_p;//based on Patankar's equation (6.24), page 128
	static double alpha_v;//relaxation value for V
	static double ConvergenceTolVx;
	static double ConvergenceTolVy;
	static double ConvergenceTolP;
	static double ConvergenceTolT;
	static bool change_alpha;
	FVM_Grid(const bGrid& G);
	~FVM_Grid();
	void SetThermalConductionProblem(double conductivity);
	void SetThermalConductionConvectionProblem(double conductivity, double Density);
	void SetThermalConductionConvectionProblem(double conductivity, double Density, function<double(const Vector3D& P)> Velocity[2]);
	void SetTVD(TVD::MODE TVD_Type = TVD::VAN_ALBADA);
	void SetFlowProblem(const LiquidProperties& liq);
	void SetFlowProblem(double Re);
	void SetFlowForcedConvectionProblem(const LiquidProperties& liq);
	void SetFlowForcedConvectionProblem(double Re, double Pe);
	void SetFlowNaturalConvectionProblem(double Ra, double Pr);
	void SetThermalBoundaryConditions(function<double(const Vector3D& P)> BC, GRID_DIRECTION side = ANY_BOUNDARY);
	void SetThermalGradientBoundaryConditions(function<double(const Vector3D& P)> BC, GRID_DIRECTION side = ANY_BOUNDARY);
	void SetThermalBoundaryConditions(const BoundaryValues& BT);
	void SetThermalGradientBoundaryConditions(const BoundaryValues& BT);
	void SetVelocityBoundaryConditions(const BoundaryValues& Vx, const BoundaryValues& Vy);
	void SetVelocityGradientBoundaryConditions(const BoundaryValues& Vx, const BoundaryValues& Vy);
	void SolveThermalProblem(vector<Node*>& results);
	void SolveSIMPLE(vector<Node*>& Vx, vector<Node*>& Vy, vector<Node*>& P, int maxIter);
	void SolveSIMPLE(vector<Node*>& Vx, vector<Node*>& Vy, vector<Node*>& P, vector<Node*>& T, int maxIter);
	void printEquations(string fileName, int iterNum);
	void printThermalEquations(string fileName);
};
bool tester_FVM_Grid(int& NumTests);
bool tester_FVM_Grid_1(int& NumTests);
bool tester_FVM_Grid_2(int& NumTests);
bool tester_FVM_Grid_3(int& NumTests);
bool tester_FVM_Grid_7(int& NumTests);//lid-driven cavity
bool tester_FVM_Grid_9(int& NumTests);//Convective lid-driven cavity
bool tester_FVM_Grid_10(int& NumTests);//RBC Box
bool tester_FVM_Grid_11(int& NumTests);//test # 1 with Neumann's boundary condition
bool tester_FVM_Grid_12(int& NumTests);//test # 10 , infinitely extended RBC.
bool tester_FVM_Grid_13(int& NumTests);//Conduction-convection using TVD (based on tester_FVM_Grid_3)
#endif
