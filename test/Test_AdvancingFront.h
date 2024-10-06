#ifndef Test_AdvancingFrontH
#define Test_AdvancingFrontH
#include "../src/AdvancingFront.h"
#include "../src/QuadEdge.h"
bool tester_AdvancingFront_1(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 1;
	int Ny = 1;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = AdvancingFront::Tessellate(input);
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
bool tester_AdvancingFront_2(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 2;
	int Ny = 2;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = AdvancingFront::Tessellate(input);
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
bool tester_AdvancingFront_3(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 3;
	int Ny = 3;
	double Lx = 3;
	double Ly = 3;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = AdvancingFront::Tessellate(input);
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
bool tester_AdvancingFront_4(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 4;
	int Ny = 4;
	double Lx = 4;
	double Ly = 4;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = AdvancingFront::Tessellate(input);
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
bool tester_AdvancingFront_5(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 5;
	int Ny = 5;
	double Lx = 2;
	double Ly = 2;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = AdvancingFront::Tessellate(input); 
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
bool tester_AdvancingFront_6(int& NumTests)
{
	vector<Vector3D> input;
	vector<Vector3D> output;
	int Nx = 6;
	int Ny = 6;
	double Lx = 6;
	double Ly = 6;
	double hx = Lx / double(Nx);
	double hy = Ly / double(Ny);
	input.resize(2 * (Nx + Ny));
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
	Triangulation T = AdvancingFront::Tessellate(input); 
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
bool tester_AdvancingFront_7(int& NumTests)
{
	//Like test 5 but  input[2:3] are shifted upward
	vector<Vector3D> input;
	int Nx = 5;
	int Ny = 5;
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
	input.erase(input.begin() + 19);
	input.erase(input.begin() + 18);
	input.erase(input.begin() + 17);
	input.erase(input.begin() + 16);
	input[13](1) *= 2.0;
	input.erase(input.begin() + 11);
	input.erase(input.begin() + 9);
	input.erase(input.begin() + 8);
	input.erase(input.begin() + 7);
	input.erase(input.begin() + 6);
	input[2](1) = 0.5;
	input[3](1) = 0.5;
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	ofstream file;
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront_9(int& NumTests)
{
	//Like test 5 but  input[2:3] are shifted upward
	vector<Vector3D> input;
	double R1 = 1, R2 = 2;
	int Nr = 5;
	int Nt = 5;
	double dr = (R2 - R1) / double(Nr);
	double dt = (pi / 2.0) / double(Nt);
	int n = -1;
	Vector3D buffer;
	for (int i = 0; i < Nr; ++i) {
		buffer(0) = R1 + double(i) * dr; buffer(1) = 0;
		input.push_back(buffer);
	}
	for (int i = 0; i < 2 * Nt; ++i) {
		double theta = double(i) * dt / 2.0;
		buffer(0) = R2 * cos(theta); buffer(1) = R2 * sin(theta);
		input.push_back(buffer);
	}
	for (int i = 0; i < Nr; ++i) {
		buffer(0) = 0; buffer(1) = R2 - double(i) * dr;
		input.push_back(buffer);
	}
	for (int i = 0; i < Nt; ++i) {
		double theta = pi / 2.0 - double(i) * dt;
		buffer(0) = R1 * cos(theta); buffer(1) = R1 * sin(theta);
		input.push_back(buffer);
	}
	Triangulation T = AdvancingFront::Tessellate(input);
	T.PrintTriangulation();
	ofstream file;
	if (!T.TestIntegrity())
		return false;
	if (!T.TestDelaunay())
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront_10(int& NumTests)
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
	Triangulation Triang = AdvancingFront::Tessellate(input);
	Triang.PrintTriangulation();
	if (!Triang.TestIntegrity())
		return false;
	if (!Triang.TestDelaunay())
		return false;
	double alpha_max = Triang.GetMaxAngle() / pi * 180.0;
	if (fabs(alpha_max - 135) > 0.01)
		return false;
	int Nv = Triang.NumVertices();
	if (Nv != 13)
		return false;
	++NumTests;
	return true;
}
bool tester_AdvancingFront(int& NumTests)
{
	if (!tester_AdvancingFront_1(NumTests))
		return false;
	if (!tester_AdvancingFront_2(NumTests))
		return false;
	if (!tester_AdvancingFront_3(NumTests))
		return false;
	if (!tester_AdvancingFront_4(NumTests))
		return false;
	if (!tester_AdvancingFront_5(NumTests))
		return false;
	if (!tester_AdvancingFront_6(NumTests))
		return false;
	if (!tester_AdvancingFront_7(NumTests))
		return false;
	if (!tester_AdvancingFront_9(NumTests))
		return false;
	if (!tester_AdvancingFront_10(NumTests))
		return false;
	++NumTests;
	return true;
}
#endif