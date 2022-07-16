#ifndef FVMH
#define FVMH
#include "Triangulation.h"
#include "Vector3D.h"
#include <unordered_map>
#include "Node.h"
#include "NodeComposite.h"
#include "LiquidProperties.h"
#include <functional>
using namespace std;
class Vector;
class SquareMatrix;
class FVM
{
	enum MODE { NONE, CONSERVATION, SIMPLE } mode {MODE::NONE};
	Triangulation* triangulation{0};
	double constantConductivity{0};
	double constantDensity{0};
	double ConvergenceV[2] {0, 0};
	double ConvergenceP{0};
	function<double(const Vector3D& P)> conductivityFunction {0};
	function<double(const Vector3D& P)> vFunction[2] {0, 0};
	NodeComposite TNodes;
	NodeComposite VxNodes;
	NodeComposite VyNodes;
	NodeComposite PNodes;
	NodeComposite aPu;
	NodeComposite aPv;
	LiquidProperties liquid;
	void InitializeTNodes();
	void AddDiffusion_vertexBased(double conductivity);
	void AddDiffusion_cellBased(NodeComposite& nodes, function<double(const Vector3D& P)> conductivit);
	void AddDiffusionAndConvection_vertexBased(function<double(const Vector3D& P)> conductivity, double density, function<double(const Vector3D& P)> V[2]);
	Vector3D GetGradient(NodeComposite& nodes, Face* f);
	void AddConvection_cellBased(NodeComposite& nodes);
	void SolveSIMPLE_v(int i);
	void SolveSIMPLE_p();
	void SolveSIMPLE_CorrectV();
	void populateVelocities(function<double(const Vector3D& P)> V[2]);
	void populateVertexVelocities();
	void Debugger();
public: enum CELL_GEOMETRY {VERTEX_BASED, CELL_BASED};
private: CELL_GEOMETRY cell_geo;
public:
	FVM(const Triangulation& T, CELL_GEOMETRY cell_geo_ = CELL_GEOMETRY::CELL_BASED);
	~FVM();
	void AddDiffusion(double conductivity);
	void AddDiffusionAndConvection(function<double(const Vector3D& P)> conductivity, double density, function<double(const Vector3D& P)> V[2]);
	void AddDiffusionAndConvection(double conductivity, double density, function<double(const Vector3D& P)> V[2]);
	void SetThermalBoundaryConditions(function<double(const Vector3D& P)> BC);
	void SetVxBoundaryConditions(function<double(const Vector3D& P)> BC);
	void SetVyBoundaryConditions(function<double(const Vector3D& P)> BC);
	void SetFlowProblem(const LiquidProperties& liq);
	void SetFlowProblem(double Re);
	void Solve(vector<Node*>& results);
	double Get_TValue(const Vector& P);
	double GetTValue(GeoGraphObject* objPtr, bool& isValid);
	double GetVxValue(GeoGraphObject* objPtr, bool& isValid);
	double GetVyValue(GeoGraphObject* objPtr, bool& isValid);
	double GetPValue(GeoGraphObject* objPtr, bool& isValid);

};
bool tester_FVM(int& NumTests); //tester_FVM_6 has errors !!!!!!!!!!!!!!!!!!!!!!!!
bool tester_FVM_1(int& NumTests);//conduction Rectangular 2 x 2
bool tester_FVM_2(int& NumTests);//conduction Rectangular 20 x 20
bool tester_FVM_3(int& NumTests);//conduction-convection Rectangular 20 x 20
bool tester_FVM_4(int& NumTests);//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1980 Example 1 Case 1
bool tester_FVM_5(int& NumTests);//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1980 Example 1 Case 2
bool tester_FVM_6(int& NumTests);//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1980 Example 2
bool tester_FVM_6a(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6b(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6c(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6d(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6e(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_7(int& NumTests);// tester_FVM_1 using the cell_based FVM.
bool tester_FVM_8(int& NumTests);// tester_FVM_2 using the cell_based FVM.
bool tester_FVM_9(int& NumTests);// tester_FVM_3 using the cell_based FVM and TVD.
#endif