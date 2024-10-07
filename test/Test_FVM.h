#ifndef Test_FVMH
#define Test_FVMH
#include "../src/FVM.h"
#include "../src/ConductionConvectionProblem.h"
#include "../src/DelaunayLifting.h"
#include "../src/FaceIterator.h"
#include "../src/Face.h"
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
bool tester_FVM_6a(int &NumTests)
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
bool tester_FVM_6b(int &NumTests)
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
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.24)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6c(int &NumTests)
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
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.0001)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6d(int &NumTests)
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
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.23)
		return false;
	++NumTests;
	return true;
}
bool tester_FVM_6e(int &NumTests)
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
		}
	}
	file.close();
	if (max_abs_error_VerticalCenterline > 0.3)
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
		vector<double> max_error_percent(N_cases.size(), 0);
		for (int n = 0; n < N_cases.size(); ++n)
		{
			int N = N_cases[n];
			Triangulation T = Triangulation::OffDiagonalGrid(N, N, Lx, Ly);
			FVM fvm(T);
			fvm.SetTVD(TVD::MODE::NO_TVD);
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
				double error_percent = T_actual != 0 ? fabs(T - T_actual) / T_actual * 100 : fabs(T - T_actual) / T0 * 100;
				if (error_percent > max_error_percent[n])
				{
					max_error_percent[n] = error_percent;
				}
			}
			if (c == k_cases.size() - 1 && n == N_cases.size() - 1)
			{
				if (k > 10 && k < 0.1 && max_error_percent[n] > 1.0e-5)
					return false;
				if (fabs(k - 0.1) < 1.0e-10 && max_error_percent[n] > 18)
					return false;
				if (fabs(k - 1) < 1.0e-10 && max_error_percent[n] > 29)
					return false;
				if (fabs(k - 10) < 1.0e-10 && max_error_percent[n] > 5)
					return false;
			}
		}
		cout << endl
			 << "k = " << k << " , error % = " << max_error_percent << " , for mesh nxn, n = " << N_cases << endl;
	}
	++NumTests;
	return true;
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
#endif