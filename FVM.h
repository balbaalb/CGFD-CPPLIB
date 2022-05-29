#ifndef FVMH
#define FVMH
#include "Triangulation.h"
#include "Vector3D.h"
#include <unordered_map>
#include "Node.h"
#include "NodeComposite.h"
using namespace std;
class Vector;
class SquareMatrix;
class FVM
{
	Triangulation* triangulation;
	static double constant_conductivity;
	static double get_constant_conductivity(double x, double y);
	NodeComposite TNodes;
public:
	FVM(const Triangulation& T);
	~FVM();
	void AddDiffusion(double conductivity);
	void AddDiffusionAndConvection(fxy conductivity, double density, velocityField V);
	void AddDiffusionAndConvection(double conductivity, double density, velocityField V);
	void SetBoundaryConditions(fxy BC);
	void Solve(vector<Node*>& vertex_results);
	double Get_TValue(const Vector& P);
};
bool tester_FVM(int& NumTests); //tester_FVM_6 has errors !!!!!!!!!!!!!!!!!!!!!!!!
class testerFVM1
{
	static double Lx, Ly, T0, T1, A;
public:
	static void Set(double Lx_, double Ly_, double T0_, double T1_);
	static double value(double x, double y);
	static double DvalueDx(double x, double y);
	static double DvalueDy(double x, double y);
};
bool tester_FVM_1(int& NumTests);//conduction Rectangular 2 x 2
bool tester_FVM_2(int& NumTests);//conduction Rectangular 20 x 20
class testerFVM3
{
	static double T0, T1, vx, vy, Lx, Ly;
	static double kapa;//conductivity / density
public:
	static void Set(double T0_, double T1_, double vx_, double vy_, double Lx_, double Ly_, double kinematic_diffusivity);
	static double Vx(double x, double y);
	static double Vy(double x, double y);
	static double value(double x, double y);
};
bool tester_FVM_3(int& NumTests);//conduction-convection Rectangular 20 x 20
class testerFVM4
{
public:
	static double Vx(double x, double y);
	static double Vy(double x, double y);
	static double value(double x, double y);
	static double Pe, x0, y0;
	static double variableConductivity(double x, double y);
	static double value_variableConductivity(double x, double y);
};
bool tester_FVM_4(int& NumTests);//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1980 Example 1 Case 1
bool tester_FVM_5(int& NumTests);//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1980 Example 1 Case 2
class testerFVM6
{
public:
	static double Lx, Ly, alpha /* yc / Ly*/, v0, T0, T1;
	static bool averageOnBoundary;
	static double Vx(double x, double y);
	static double Vy(double x, double y);
	static double value(double x, double y);
};
bool tester_FVM_6(int& NumTests);//conduction-convection Rectangular 10 x 10, Baliga & Patankar 1980 Example 2
bool tester_FVM_6a(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6b(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6c(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6d(int& NumTests, vector<double>* verticalCenterlineResults = 0);
bool tester_FVM_6e(int& NumTests, vector<double>* verticalCenterlineResults = 0);
class testerFVM7//Lid-driven cavity
{
public:
	static double Lx, Ly, VLid, density, viscosity;
	static int Nx, Ny;
};
class testerFVM8//Duct flow 1D
{
public:
	static double Lx, Ly, density, viscosity, InletV;
	static int Nx, Ny;
};
class testerFVM9//Convective Lid Driven Cavity Flow, forced convection
{
public:
	double Re, Pe;
	int Nx, Ny;
	double RetOne();
	__declspec(property(get = RetOne)) double Lx;
	__declspec(property(get = RetOne)) double Ly;
	__declspec(property(get = RetOne)) double VLid;
	__declspec(property(get = RetOne)) double TLid;
};
#endif