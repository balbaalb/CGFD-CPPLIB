#include "ConductionConvectionProblem.h"
double ConductionConvectionProblem::T(const Vector3D& P) const
{
	return T(P(0), P(1), P(2));
}

double ConductionConvectionProblem::vx(const Vector3D& P) const
{
	return vx(P(0), P(1), P(2));
}

double ConductionConvectionProblem::vy(const Vector3D& P) const
{
	return vy(P(0), P(1), P(2));
}

ConductionExample::ConductionExample(double Lx_, double Ly_, double T0_, double T1_) :
	Lx(Lx_), Ly(Ly_), T0(T0_), T1(T1_)
{
	A = pi * Ly / Lx;
}
double ConductionExample::T(double x, double y, double z) const
{
	return T0 + T1 * sin(pi * x / Lx) * (exp(A * y / Ly) - exp(-A * y / Ly));
}

double ConductionExample::T(const Vector3D& P) const
{
	return T(P(0), P(1), P(2));
}

double ConductionExample::vx(double x, double y, double z) const
{
	return 0.0;
}

double ConductionExample::vx(const Vector3D& P) const
{
	return 0.0;
}

double ConductionExample::vy(double x, double y, double z) const
{
	return 0.0;
}

double ConductionExample::vy(const Vector3D& P) const
{
	return 0.0;
}

double ConductionExample::dT_dx(double x, double y, double z) const
{
	return T1 * pi / Lx * cos(pi * x / Lx) * (exp(A * y / Ly) - exp(-A * y / Ly));
}

double ConductionExample::dT_dx(const Vector3D& P) const
{
	return dT_dx(P(0), P(1), P(2));
}

double ConductionExample::dT_dy(double x, double y, double z) const
{
	return T1 * sin(pi * x / Lx) * A / Ly * (exp(A * y / Ly) + exp(-A * y / Ly));
}

double ConductionExample::dT_dy(const Vector3D& P) const
{
	return dT_dy(P(0), P(1), P(2));
}

convectionConstVelocity::convectionConstVelocity(double T0_, double T1_, double vx_, double vy_, double Lx_, double Ly_, double kinematic_diffusivity) :
	T0(T0_), T1(T1_), constVx(vx_), constVy(vy_), Lx(Lx_), Ly(Ly_), kapa(kinematic_diffusivity)
{

}

double convectionConstVelocity::T(double x, double y, double z) const
{
	double argument = (x * constVx + y * constVy) / kapa;
	double maxArgument = (Lx * constVx + Ly * constVy) / kapa;
	double T = T0 + (T1 - T0) * (exp(argument) - 1.0) / (exp(maxArgument) - 1.0);
	return T;
}

double convectionConstVelocity::T(const Vector3D& P) const
{
	return T(P(0), P(1), P(2));
}

double convectionConstVelocity::vx(double x, double y, double z) const
{
	return constVx;
}

double convectionConstVelocity::vx(const Vector3D& P) const
{
	return constVx;
}

double convectionConstVelocity::vy(double x, double y, double z) const
{
	return constVy;
}

double convectionConstVelocity::vy(const Vector3D& P) const
{
	return constVy;
}

BaligaPatankarExample1::BaligaPatankarExample1(double Pe_) : Pe(Pe_), x0(sqrt(2.0) / 2.0), y0(sqrt(2.0) / 2.0)
{
}

BaligaPatankarExample1::BaligaPatankarExample1(double Pe_, double x0_, double y0_) : Pe(Pe_), x0(x0_), y0(y0_)
{
}

double BaligaPatankarExample1::T(double x, double y, double z) const
{
	x += x0;
	y += y0;
	double T = 1.0 - log(x * x + y * y) / 2.0 / log(3.0);
	return T;
}

double BaligaPatankarExample1::T(const Vector3D& P) const
{
	return T(P(0), P(1), P(2));
}

double BaligaPatankarExample1::vx(double x, double y, double z) const
{
	x += x0;
	y += y0;
	return 2.0 * y;
}

double BaligaPatankarExample1::vx(const Vector3D& P) const
{
	return vx(P(0), P(1), P(2));
}

double BaligaPatankarExample1::vy(double x, double y, double z) const
{
	x += x0;
	y += y0;
	return -2.0 * x;
}

double BaligaPatankarExample1::vy(const Vector3D& P) const
{
	return vy(P(0), P(1), P(2));
}

BaligaPatankarExample1_case2::BaligaPatankarExample1_case2(double Pe_) : BaligaPatankarExample1(Pe_)
{
}

double BaligaPatankarExample1_case2::T(double x, double y, double z) const
{
	x += x0;
	y += y0;
	double T = (9.0 - (x * x + y * y)) / 8.0;
	return T;
}

double BaligaPatankarExample1_case2::T(const Vector3D& P) const
{
	return T(P(0), P(1), P(2));
}

double BaligaPatankarExample1_case2::variableConductivity(double x, double y, double z) const
{
	x += x0;
	y += y0;
	return 1.0 / Pe / (x * x + y * y);
}

double BaligaPatankarExample1_case2::variableConductivity(const Vector3D& P) const
{
	return variableConductivity(P(0), P(1), P(2));
}

double BaligaPatankarExample2::T(double x, double y, double z) const
{
	double f = alpha - y / Ly + x / Lx * (1.0 - 2.0 * alpha);
	if (f < -1e-10)/*above*/
		return T1;
	if (f > 1e-10)/*below*/
		return T0;
	if (averageOnBoundary)
		return (T0 + T1) / 2.0;
	return T0;
}

double BaligaPatankarExample2::T(const Vector3D& P) const
{
	return T(P(0), P(1), P(2));
}

double BaligaPatankarExample2::vx(double x, double y, double z) const
{
	double dx = Lx;
	double dy = Ly * (1.0 - 2.0 * alpha);
	return v0 * dx / sqrt(dx * dx + dy * dy);
}

double BaligaPatankarExample2::vx(const Vector3D& P) const
{
	return vx(P(0), P(1), P(2));
}

double BaligaPatankarExample2::vy(double x, double y, double z) const
{
	double dx = Lx;
	double dy = Ly * (1.0 - 2.0 * alpha);
	return v0 * dy / sqrt(dx * dx + dy * dy);
}

double BaligaPatankarExample2::vy(const Vector3D& P) const
{
	return vy(P(0), P(1), P(2));
}



