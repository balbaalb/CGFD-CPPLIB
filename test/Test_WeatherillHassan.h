#ifndef Test_WeatherillHassanH
#define Test_WeatherillHassanH
#include "../src/WeatherillHassan.h"
#include "../src/QuadEdge.h"
bool test_WeatherillHassan_1(int& NumTests)
{
	vector<Vector3D> input;
	input.resize(8);
	input[0](0) = 0; input[0](1) = 0;
	input[1](0) = 1.0; input[1](1) = 0;
	input[2](0) = 2.0; input[2](1) = 0;
	input[3](0) = 2.0; input[3](1) = 1.0;
	input[4](0) = 2.0; input[4](1) = 2.0;
	input[5](0) = 1.0; input[5](1) = 2.0;
	input[6](0) = 0; input[6](1) = 2.0;
	input[7](0) = 0; input[7](1) = 1.0;
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation T = WeatherillHassan::Tessellate(input,100);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	++NumTests;
	return true;
}
bool test_WeatherillHassan_2(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 4;
	int Ny = 4;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2*(Nx + Ny));
	for (int i = 0; i < Nx; ++i){
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i){
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i){
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i){
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation T = WeatherillHassan::Tessellate(input, 100);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;
	++NumTests;
	return true;
}
bool test_WeatherillHassan_3(int& NumTests)
{
	vector<Vector3D> input;
	int Nx = 4;
	int Ny = 4;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
	for (int i = 0; i < Nx; ++i) {
		input[i](0) = hx * i; input[i](1) = 0;
	}
	for (int i = 0; i < Ny; ++i) {
		int n = Nx + i;
		input[n](0) = Lx; input[n](1) = i * hy;
	}
	for (int i = 0; i < Nx; ++i) {
		int n = Nx + Ny + i;
		input[n](0) = Lx - i * hx; input[n](1) = Ly;
	}
	for (int i = 0; i < Ny; ++i) {
		int n = 2 * Nx + Ny + i;
		input[n](0) = 0; input[n](1) = Ly - i * hy;
	}
	//Making the boundary non-convex:----
	int n1 = Nx / 2;
	input[n1](1) = hy;
	int n2 = Nx + Ny / 2;
	input[n2](0) = Lx + hx;
	//-----------------------------------
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation T = WeatherillHassan::Tessellate(input, 100);
	T.PrintTriangulation();
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	QuadEdge* mesh = T.GetMesh2D();
	if (!mesh->TestIntegrity())
		return false;/**/
	++NumTests;
	return true;
}
bool test_WeatherillHassan_4(int& NumTests)
{
	int N = 4;
	vector<Vector3D> input;
	Vector3D A, B(1, 0), C(2, 1);
	Vector3D P = A;
	Vector3D dP = B - A;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	P = B;
	dP = C - B;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	P = C;
	dP = A - C;
	dP = dP / double(N);
	for (int i = 0; i < N; ++i) {
		input.push_back(P);
		P = P + dP;
	}
	WeatherillHassan::alpha = 1.0;
	WeatherillHassan::beta = 0.1;
	Triangulation Triang = WeatherillHassan::Tessellate(input, 100);
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	if (fabs(alpha_max - 135) > 0.01)
		return false;
	int Nv = Triang.NumVertices();
	if (Nv != 12)
		return false;
	++NumTests;
	return true;
}
bool test_WeatherillHassan(int& NumTests)
{
	if (!test_WeatherillHassan_1(NumTests))
		return false;
	if (!test_WeatherillHassan_2(NumTests))
		return false;
	if (!test_WeatherillHassan_3(NumTests))
		return false;
	if (!test_WeatherillHassan_4(NumTests))
		return false;
	NumTests += 1;
	return true;
}
#endif