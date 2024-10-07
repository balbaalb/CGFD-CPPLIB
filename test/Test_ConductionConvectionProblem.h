#ifndef Test_ConductionConvectionProblemH
#define Test_ConductionConvectionProblemH
#include "../src/ConductionConvectionProblem.h"
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
bool tester_ConductionConvectionProblem(int& NumTests)
{
	if(!tester_convectionConstVelocity(NumTests))
		return false;
	++NumTests;
	return true;
}
#endif