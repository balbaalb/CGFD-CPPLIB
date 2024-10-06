#include "ConductionConvectionProblem.h"
double ConductionConvectionProblem::T(const Vector3D &P) const
{
	return T(P(0), P(1), P(2));
}

double ConductionConvectionProblem::vx(const Vector3D &P) const
{
	return vx(P(0), P(1), P(2));
}

double ConductionConvectionProblem::vy(const Vector3D &P) const
{
	return vy(P(0), P(1), P(2));
}

ConductionExample::ConductionExample(double Lx_, double Ly_, double T0_, double T1_) : Lx(Lx_), Ly(Ly_), T0(T0_), T1(T1_)
{
	A = pi * Ly / Lx;
}
double ConductionExample::T(double x, double y, double z) const
{
	return T0 + T1 * sin(pi * x / Lx) * (exp(A * y / Ly) - exp(-A * y / Ly));
}

double ConductionExample::T(const Vector3D &P) const
{
	return T(P(0), P(1), P(2));
}

double ConductionExample::vx(double x, double y, double z) const
{
	return 0.0;
}

double ConductionExample::vx(const Vector3D &P) const
{
	return 0.0;
}

double ConductionExample::vy(double x, double y, double z) const
{
	return 0.0;
}

double ConductionExample::vy(const Vector3D &P) const
{
	return 0.0;
}

double ConductionExample::dT_dx(double x, double y, double z) const
{
	return T1 * pi / Lx * cos(pi * x / Lx) * (exp(A * y / Ly) - exp(-A * y / Ly));
}

double ConductionExample::dT_dx(const Vector3D &P) const
{
	return dT_dx(P(0), P(1), P(2));
}

double ConductionExample::dT_dy(double x, double y, double z) const
{
	return T1 * sin(pi * x / Lx) * A / Ly * (exp(A * y / Ly) + exp(-A * y / Ly));
}

double ConductionExample::dT_dy(const Vector3D &P) const
{
	return dT_dy(P(0), P(1), P(2));
}

convectionConstVelocity::convectionConstVelocity(double T0_, double T1_, double vx_, double vy_, double Lx_, double Ly_, double kinematic_diffusivity, double y0_) : 
	T0(T0_), T1(T1_), constVx(vx_), constVy(vy_), Lx(Lx_), Ly(Ly_), kappa(kinematic_diffusivity), y0(y0_)
{
	Amax = fmax(Amax , fmax(Lx * constVx / kappa , fmax(Ly * constVy / kappa , (Lx * constVx + Ly * constVy) / kappa)));
}

double convectionConstVelocity::T(double x, double y, double z) const
{
	if (kappa > 0)
	{
		double argument = (x * constVx + y * constVy) / kappa;
		double maxArgument = (fabs(Lx * constVx) + fabs(Ly * constVy)) / kappa;
		double T = T0 + (T1 - T0) * exp(argument - maxArgument);
		return T;
	}
	else if (fabs(constVy) > 0)
	{
		double m = constVx / constVy;
		double f = m * x + y0 - y;
		return (f >= 0 ? T0 : T1);
	}
	return (y <= y0 ? T0 : T1);
}

double convectionConstVelocity::T(const Vector3D &P) const
{
	return T(P(0), P(1), P(2));
}

double convectionConstVelocity::vx(double x, double y, double z) const
{
	return constVx;
}

double convectionConstVelocity::vx(const Vector3D &P) const
{
	return constVx;
}

double convectionConstVelocity::vy(double x, double y, double z) const
{
	return constVy;
}

double convectionConstVelocity::vy(const Vector3D &P) const
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

double BaligaPatankarExample1::T(const Vector3D &P) const
{
	return T(P(0), P(1), P(2));
}

double BaligaPatankarExample1::vx(double x, double y, double z) const
{
	x += x0;
	y += y0;
	return 2.0 * y;
}

double BaligaPatankarExample1::vx(const Vector3D &P) const
{
	return vx(P(0), P(1), P(2));
}

double BaligaPatankarExample1::vy(double x, double y, double z) const
{
	x += x0;
	y += y0;
	return -2.0 * x;
}

double BaligaPatankarExample1::vy(const Vector3D &P) const
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

double BaligaPatankarExample1_case2::T(const Vector3D &P) const
{
	return T(P(0), P(1), P(2));
}

double BaligaPatankarExample1_case2::variableConductivity(double x, double y, double z) const
{
	x += x0;
	y += y0;
	return 1.0 / Pe / (x * x + y * y);
}

double BaligaPatankarExample1_case2::variableConductivity(const Vector3D &P) const
{
	return variableConductivity(P(0), P(1), P(2));
}

double BaligaPatankarExample2::T(double x, double y, double z) const
{
	double f = alpha - y / Ly + x / Lx * (1.0 - 2.0 * alpha);
	if (f < -1e-10) /*above*/
		return T1;
	if (f > 1e-10) /*below*/
		return T0;
	if (averageOnBoundary)
		return (T0 + T1) / 2.0;
	return T0;
}

double BaligaPatankarExample2::T(const Vector3D &P) const
{
	return T(P(0), P(1), P(2));
}

double BaligaPatankarExample2::vx(double x, double y, double z) const
{
	double dx = Lx;
	double dy = Ly * (1.0 - 2.0 * alpha);
	return v0 * dx / sqrt(dx * dx + dy * dy);
}

double BaligaPatankarExample2::vx(const Vector3D &P) const
{
	return vx(P(0), P(1), P(2));
}

double BaligaPatankarExample2::vy(double x, double y, double z) const
{
	double dx = Lx;
	double dy = Ly * (1.0 - 2.0 * alpha);
	return v0 * dy / sqrt(dx * dx + dy * dy);
}

double BaligaPatankarExample2::vy(const Vector3D &P) const
{
	return vy(P(0), P(1), P(2));
}

bool tester_ConductionConvectionProblem(int& NumTests)
{
	if(!tester_convectionConstVelocity(NumTests))
		return false;
	++NumTests;
	return true;
}

bool tester_convectionConstVelocity(int& NumTests)
{
	double T_00 {1.0}; // T(0,0)
	double T1 {10} , vx {1.25} , vy {2.5} , Lx {1.0} , Ly {1.5} , kappa {3.0};
	double Amax = (Lx * vx + Ly * vy) / kappa;
	double T0 = (T_00 - T1 * exp(-Amax)) / (1 - exp(-Amax));
	convectionConstVelocity problem(T0, T1, vx, vy, Lx, Ly, kappa);
	Vector3D P0(0 , 0), P1(Lx , Ly), P10(Lx , 0) , P01(0 , Ly);
	Vector3D Pm = (P0 + P1) / 2;
	if(problem.vx(P0) != 1.25 || problem.vx(P1) != 1.25 || problem.vx(Pm) != 1.25 || problem.vx(P10) != 1.25 || problem.vx(P01) != 1.25)
		return false;
	if(problem.vy(P0) != 2.5 || problem.vy(P1) != 2.5 || problem.vy(Pm) != 2.5 || problem.vy(P10) != 2.5 || problem.vy(P01) != 2.5)
		return false;
	if(fabs(problem.T(P0) - T_00) > 1.0e-10)
		return false;
	if(fabs(problem.T(P1) - T1) > 1.0e-10)
		return false;
	double T_10 = T0 + (T1 - T0) * exp(Lx * vx / kappa - Amax);
	if(fabs(problem.T(P10) - T_10) > 1.0e-10)
		return false;
	double T_01 = T0 + (T1 - T0) * exp(Ly * vy / kappa - Amax);
	if(fabs(problem.T(P01) - T_01) > 1.0e-10)
		return false;
	double Tm = T0 + (T1 - T0) * exp(-0.5 * Amax);
	if(fabs(problem.T(Pm) - Tm) > 1.0e-10)
		return false;
	++NumTests;
	return true;
}
