#ifndef FVMH
#define FVMH
#include "Triangulation.h"
#include "Vector3D.h"
#include <unordered_map>
#include "Node.h"
#include "NodeComposite.h"
#include <functional>
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
	void AddDiffusionAndConvection(function<double(const Vector3D& P)> conductivity, double density, function<double(const Vector3D& P)> V[2]);
	void AddDiffusionAndConvection(double conductivity, double density, function<double(const Vector3D& P)> V[2]);
	void SetBoundaryConditions(function<double(const Vector3D& P)> BC);
	void Solve(vector<Node*>& vertex_results);
	double Get_TValue(const Vector& P);
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
#endif