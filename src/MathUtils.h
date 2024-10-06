#ifndef MathUtilsH
#define MathUtilsH
#define _USE_MATH_DEFINES
#include <math.h>
#define pi M_PI
void SolveQuad(double a, double b, double c, double& Root1Real, double& Root1imag, double& Root2Real, double& Root2imag);
bool IsEqual(double a, double b);
bool IsNotEqual(double a, double b);
//=======================================================================================
bool tester_MathUtils(int& NumTests);
#endif
