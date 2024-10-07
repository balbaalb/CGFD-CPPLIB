#ifndef Test_FVM_GridH
#define Test_FVM_GridH
#include "../src/FVM_Grid.h"
#include "../src/ConductionConvectionProblem.h"
#include "../src/QuadEdge.h"
bool tester_FVM_Grid_1(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	bGrid G(1, 1, Lx, Ly);
	FVM_Grid fvm(G);
	fvm.SetThermalConductionProblem(k);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
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
	if (max_error_percent > 25)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_Grid_2(int& NumTests)
{
	double Lx = 2, Ly = 2, T0 = 1.0, T1 = 10;
	double k = 1.0;
	ConductionExample cond(Lx, Ly, T0, T1);
	bGrid G(10, 10, Lx, Ly);
	FVM_Grid fvm(G);
	fvm.SetThermalConductionProblem(k);
	auto temperatureField = [cond](const Vector3D& P) {return cond.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
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
	convectionConstVelocity conv(T0, T1, vx, vy, Lx, Ly, kapa);
	bGrid G(N, N, Lx, Ly);
	G.GetMesh2D()->print("..\\Data\\tester_FVM_Grid_3.Grid.txt");
	FVM_Grid fvm(G);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.SetThermalConductionConvectionProblem(k, rho, V);
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		double T_actual = conv.T(results[i]->GetPoint());
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
	double Lx = 1.0;
	double Ly = 1.0;
	double density = 1;
	double viscosity = 1.0;
	double VLid = 1;
	int Nx = 6;
	int Ny = 6;
	int maxIter = 1000;
	BoundaryValues VBx(0), VBy(0);
	*VBx.boundaryValue[NORTH] = VLid;
	LiquidProperties liq;
	liq.density = density;
	liq.viscosity = viscosity;
	bGrid G(Nx, Ny, Lx, Ly);
	G.GetMesh2D()->print("..\\Data\\tester_FVM_Grid_7.Grid.txt");
	FVM_Grid::alpha_p = 0.3;
	FVM_Grid::alpha_v = 0.3;
	FVM_Grid fvm(G);
	fvm.SetFlowProblem(liq);
	fvm.SetVelocityBoundaryConditions(VBx, VBy);
	fvm.SetMinconvergenceP(0.1);
	fvm.SetMinconvergenceV(0.1);
	fvm.SetMinconvergenceT(0.1);
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
	if (Nx >= 6 && Ny >= 6 && dP_max_error > 9)
		return false;
	double Umin_ANSYS = -0.1455;
	double Umin_error = fabs((Umin - Umin_ANSYS) / Umin_ANSYS) * 100.0;
	double Vmin_ANSYS = -0.2015;
	double Vmin_error = fabs((Vmin - Vmin_ANSYS) / Vmin_ANSYS) * 100.0;
	double Vmax_ANSYS = 0.2007;
	double Vmax_error = fabs((Vmax - Vmax_ANSYS) / Vmax_ANSYS) * 100.0;
	if (Nx >= 6 && Ny >= 6)
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
	//testerFVM9 test9;
	double Re = 100;
	double Pe = Re;
	int Nx = 6;
	int Ny = 6;
	double Lx = 1.0;
	double Ly = 1.0;
	int maxIter = 1000;
	BoundaryValues VBx(0), VBy(0), TB(0);
	*VBx.boundaryValue[NORTH] = 1.0;
	*TB.boundaryValue[NORTH] = 1.0;
	bGrid G(Nx, Ny, Lx, Ly);
	G.GetMesh2D()->print("..\\Data\\tester_FVM_Grid_7.Grid.txt");
	FVM_Grid::alpha_p = 0.3;
	FVM_Grid::alpha_v = 0.3;
	FVM_Grid fvm(G);
	fvm.SetFlowForcedConvectionProblem(Re, Pe);
	fvm.SetVelocityBoundaryConditions(VBx, VBy);
	fvm.SetThermalBoundaryConditions(TB);
	fvm.SetMinconvergenceP(0.1);
	fvm.SetMinconvergenceV(0.1);
	fvm.SetMinconvergenceT(0.1);
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
	FVM_Grid fvm(G); 
	fvm.SetFlowNaturalConvectionProblem(Ra, Pr);
	BoundaryValues ZeroBC(0);
	fvm.SetVelocityBoundaryConditions(ZeroBC, ZeroBC);
	fvm.SetThermalBoundaryConditions(ZeroBC);
	fvm.SetMinconvergenceP(0.1);
	fvm.SetMinconvergenceV(0.1);
	fvm.SetMinconvergenceT(0.1);
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
	ConductionExample cond(Lx, Ly, T0, T1);
	bGrid G(Nx, Ny, Lx, Ly);
	FVM_Grid fvm(G);
	fvm.SetThermalConductionProblem(k);
	auto lambdaBC = [&cond](const Vector3D& P) {return cond.T(P); };
	auto lambdaBCx = [&cond](const Vector3D& P) {return cond.dT_dx(P); };
	auto lambdaBCy = [&cond](const Vector3D& P) {return cond.dT_dy(P); };
	fvm.SetThermalBoundaryConditions(lambdaBC, WEST);
	fvm.SetThermalBoundaryConditions(lambdaBC, SOUTH);
	fvm.SetThermalGradientBoundaryConditions(lambdaBCy, NORTH);
	fvm.SetThermalGradientBoundaryConditions(lambdaBCx, EAST);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
	fvm.printThermalEquations("..\\Data\\FVM_Grid.output.txt");
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
	fvm.SetMinconvergenceP(0.1);
	fvm.SetMinconvergenceV(0.1);
	fvm.SetMinconvergenceT(0.1);
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
bool tester_FVM_Grid_13(int& NumTests)
{
	//conduction-convection Rectangular 2 x 2 + TVD
	int N = 10;
	double Lx = 10, Ly = 10;
	double k = 10.0;
	double rho = 1.0;
	double vx = 1.5, vy = 2.5, T0 = 1.0, T1 = 1;
	double kapa = k / rho;
	convectionConstVelocity conv(T0, T1, vx, vy, Lx, Ly, kapa);
	bGrid G(N, N, Lx, Ly);
	G.GetMesh2D()->print("..\\Data\\tester_FVM_Grid_3.Grid.txt");
	FVM_Grid fvm(G);
	function<double(const Vector3D& P)> V[2];
	V[0] = [conv](const Vector3D& P) {return conv.vx(P); };
	V[1] = [conv](const Vector3D& P) {return conv.vy(P); };
	fvm.SetThermalConductionConvectionProblem(k, rho, V);
	fvm.SetTVD(); //Change from tester_FVM_Grid_3
	auto temperatureField = [conv](const Vector3D& P) {return conv.T(P); };
	fvm.SetThermalBoundaryConditions(temperatureField);
	vector<Node*> results;
	fvm.SolveThermalProblem(results);
	double max_abs_error = 0;
	double max_error_percent = 0;
	for (int i = 0; i < results.size(); ++i)
	{
		double T = results[i]->value;
		double T_actual = conv.T(results[i]->GetPoint());
		double abs_error = fabs(T - T_actual);
		double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
		if (abs_error > max_abs_error)
			max_abs_error = abs_error;
		if (error_percent > max_error_percent)
			max_error_percent = error_percent;
	}
	if (max_error_percent > 1.0e-10)//TVD Effect!
		return false;
	++NumTests;
	return true;
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
	if (!tester_FVM_Grid_13(NumTests))
		return false;
	++NumTests;
	return true;
}
#endif