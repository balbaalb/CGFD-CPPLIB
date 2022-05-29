#ifndef FVM_GridH
#define FVM_GridH
#include <string>
#include <fstream>
#include "bGrid.h"
#include "Vector3D.h"
#include <unordered_map>
#include "Node.h"
#include "NodeComposite.h"
using namespace std;
class Face;
class Edge;
class Vertex;
class SquareMatrix;
typedef double(*fxy)(double x, double y);
typedef fxy velocityField[2];
struct LiquidProperties
{
	double density, viscosity, thermalConductionCoeff, heatCapacity;
};
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
	double ConvergedVx;
	double ConvergedVy;
	double ConvergedP;
	double Source_VxT;//source term in Vx equation: Source_VxT * T; 
	double Source_VyT;//source term in Vy equation: Source_VyT * T;
	double Source_TVx;//source term in Energy equation: Source_TVx * Vx;
	double Source_TVy;//source term in Energy equation: Source_TVy * Vy;
	void AddDiffusionEquations();
	void AddDiffusionConvectionEquations();
	void SolveSIMPLE_vx();
	void SolveSIMPLE_vy();
	void SolveSIMPLE_p();
	void SolveSIMPLE_CorrectV();
	void SetThermalBoundaryConditions_internal(fxy BC, GRID_DIRECTION side,string type/*Drichlit or Neumann*/);
	static double Lx;
	static double Ly;
	static double RBC_initial_u(double x, double y);
	static double RBC_initial_v(double x, double y);
	static double RBC_initial_T(double x, double y);
public:
	static double alpha_p;//based on Patankar's equation (6.24), page 128
	static double alpha_v;//relaxation value for V
	static double ConvergenceTolVx;
	static double ConvergenceTolVy;
	static double ConvergenceTolP;
	static bool change_alpha;
	FVM_Grid(const bGrid& G);
	~FVM_Grid();
	void SetThermalConductionProblem(double conductivity);
	void SetThermalConductionConvectionProblem(double conductivity, double Density);
	void SetThermalConductionConvectionProblem(double conductivity, double Density, velocityField Velocity);
	void SetFlowProblem(const LiquidProperties& liq);
	void SetFlowProblem(double Re);
	void SetFlowForcedConvectionProblem(const LiquidProperties& liq);
	void SetFlowForcedConvectionProblem(double Re, double Pe);
	void SetFlowNaturalConvectionProblem(double Ra, double Pr);
	void SetThermalBoundaryConditions(fxy BC, GRID_DIRECTION side = ANY_BOUNDARY);
	void SetThermalGradientBoundaryConditions(fxy BC, GRID_DIRECTION side = ANY_BOUNDARY);
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
#endif
