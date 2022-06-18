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
#include "ConductionConvectionProblem.h"
double FVM::constant_conductivity = 0;
double FVM::get_constant_conductivity(double x, double y)
{
	return FVM::constant_conductivity;
}
FVM::FVM(const Triangulation& T)
{
	this->triangulation = new Triangulation(T);
	QuadEdge* qe = this->triangulation->GetMesh2D();
	Face* fb = this->triangulation->GetBoundary();
	this->TNodes.Initialize(qe, fb, VERTICES);
	this->TNodes.InitializeEquations();
}
FVM::~FVM()
{
	if (this->triangulation)
		delete this->triangulation;
	this->triangulation = 0;
}
void FVM::AddDiffusion(double conductivity)
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
void FVM::AddDiffusionAndConvection(function<double(const Vector3D& P)> conductivity, double density, function<double(const Vector3D& P)> V[2])
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
void FVM::AddDiffusionAndConvection(double conductivity, double density, function<double(const Vector3D& P)> V[2])
{
	auto constant_conductivity = [conductivity](const Vector3D& P) {return conductivity; };
	this->AddDiffusionAndConvection(constant_conductivity, density, V);
}
void FVM::SetBoundaryConditions(function<double(const Vector3D& P)> BC)
{
	int N = this->triangulation->NumVertices(); 
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
}
void FVM::Solve(vector<Node*>& vertex_results)
{
	this->TNodes.SolveAndUpdate();
	this->TNodes.Populate();
	this->TNodes.GetResults_Vertices(vertex_results);
}
double FVM::Get_TValue(const Vector& P)
{
	return this->TNodes.GetValue(P);
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
	++NumTests;
	return true;
}
bool tester_FVM_1(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	Triangulation T = Triangulation::OffDiagonalGrid(2, 2, Lx, Ly);
	FVM fvm(T);
	fvm.AddDiffusion(k);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetBoundaryConditions(temperatureField);
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
	FVM fvm(T);
	fvm.AddDiffusion(k);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetBoundaryConditions(temperatureField);
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
	FVM fvm(T);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k,rho,V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetBoundaryConditions(temperatureField);
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
	FVM fvm(T);
	BaligaPatankarExample1 conv;
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetBoundaryConditions(temperatureField);
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
	FVM fvm(T);
	BaligaPatankarExample1_case2 conv;
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	auto variableConductivity = [conv](const Vector3D& P) {return conv.variableConductivity(P); };
	fvm.AddDiffusionAndConvection(variableConductivity, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetBoundaryConditions(temperatureField);
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
	FVM fvm(T);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetBoundaryConditions(temperatureField);
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
	FVM fvm2(T);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm2.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm2.SetBoundaryConditions(temperatureField);
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
	FVM fvm3(T);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm3.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm3.SetBoundaryConditions(temperatureField);
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
	FVM fvm4(T);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm4.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm4.SetBoundaryConditions(temperatureField);
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
	FVM fvm5(T);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm5.AddDiffusionAndConvection(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm5.SetBoundaryConditions(temperatureField);
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

