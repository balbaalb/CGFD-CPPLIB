#ifndef MathUtilsH
#define MathUtilsH
#define _USE_MATH_DEFINES
#include <math.h>
#define pi M_PI
typedef double(*fxy)(double x, double y);
typedef double(*fxyz)(double x, double y, double z);
typedef fxy velocityField[2];
void SolveQuad(double a, double b, double c, double& Root1Real, double& Root1imag, double& Root2Real, double& Root2imag);
void swap(double& a, double& b);
bool IsEqual(double a, double b);
bool IsNotEqual(double a, double b);
//=======================================================================================
bool tester_MathUtils(int& NumTests);
#endif
