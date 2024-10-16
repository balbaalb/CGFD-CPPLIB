#ifndef FVMH
#define FVMH
#include <functional>
#include <unordered_map>
#include <memory>
#include "Triangulation.h"
#include "Vector3D.h"
#include "Node.h"
#include "NodeComposite.h"
#include "LiquidProperties.h"
#include "TVD.h"
using namespace std;
class Vector;
class SquareMatrix;
class FVM
{
	enum MODE { NONE, CONSERVATION, SIMPLE } mode {MODE::NONE};
	shared_ptr<Triangulation> triangulation{nullptr};
	double constantConductivity{0};
	double constantDensity{0};
	double ConvergenceV[2] {0, 0};
	double ConvergenceP{0};
	double MinConvergenceV{ 0.001 };
	double MinConvergenceP{ 0.001 };
	double MinConvergenceT{ 0.001 };
	double alpha_p{ 0.5 }, alpha_v{ 0.5 };
	function<double(const Vector3D& P)> conductivityFunction {0};
	function<double(const Vector3D& P)> vFunction[2] {0, 0};
	NodeComposite TNodes;
	NodeComposite VxNodes;//on edges keeps ufx
	NodeComposite VyNodes;//on edges keeps ufy
	NodeComposite PNodes;
	NodeComposite Helper;/*on faces keeps ap of momentum equations, 
						   on edges keeps uf* (uf with u_faces and grad p terms but without p terms).
						   See eq 11.90 of Versteeg & Malalasekera (2007)
						   Normal defined as n = e x k */
	NodeComposite gradPx;
	NodeComposite gradPy;
	LiquidProperties liquid;
	TVD::MODE TVD_mode {TVD::MODE::VAN_ALBADA};
	void InitializeTNodes();
	void AddDiffusion_vertexBased(double conductivity);
	void AddDiffusion_cellBased(NodeComposite& nodes, function<double(const Vector3D& P)> conductivit);
	void AddDiffusionAndConvection_vertexBased(function<double(const Vector3D& P)> conductivity, double density, function<double(const Vector3D& P)> V[2]);
	Vector3D GetGradient(NodeComposite& nodes, Face* f);
	Vector3D GetGradient_p(Face* f);
	void AddConvection_cellBased(NodeComposite& nodes);
	void SolveSIMPLE_v();
	void SolveSIMPLE_updateHelper();
	void SolveSIMPLE_p();
	void SolveSIMPLE_gradP();
	void SolveSIMPLE_CorrectV();
	void MakePressuresPositive();
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
	void SetTVD(TVD::MODE mode);
	void SetFlowProblem(const LiquidProperties& liq);
	void SetFlowProblem(double Re);
	void SetMinconvergenceV(double value);
	void SetMinconvergenceP(double value);
	void SetMinconvergenceT(double value);
	void SetSolveMethod(NodeComposite::METHOD value, double tolerance = 1.0e-6);
	void Solve(vector<Node*>& results);
	void Solve(vector<Node*>& Vx, vector<Node*>& Vy, vector<Node*>& P);
	QuadEdge* GetMesh2D();
	Face* GetBoundary();
	double Get_TValue(const Vector& P);
	double GetTValue(GeoGraphObject* objPtr, bool& isValid);
	double GetVxValue(GeoGraphObject* objPtr, bool& isValid);
	double GetVyValue(GeoGraphObject* objPtr, bool& isValid);
	double GetPValue(GeoGraphObject* objPtr, bool& isValid);
	void SetAlphaP(double value);
	void SetalphaV(double value);
};
#endif