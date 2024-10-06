#ifndef Test_MathUtilsH
#define Test_MathUtilsH
#include "../src/MathUtils.h"
bool tester_MathUtils(int& NumTests)
{
	double u = 1550.126;
	if (!IsEqual(u, 1550.126) || IsNotEqual(u, 1550.126))
		return false; 
	double a, b, c, Root1Real, Root1imag, Root2Real, Root2imag;
	a = b = c = Root1Real = Root1imag = Root2Real = Root2imag = 0;
	a = 1; b = -8; c = 15;
	SolveQuad(a, b, c, Root1Real, Root1imag, Root2Real, Root2imag);
	if (IsNotEqual(Root1Real, 3) || IsNotEqual(Root1imag, 0) || IsNotEqual(Root2Real, 5) || IsNotEqual(Root2imag, 0))
		return false;
	a = 1; b = 0; c = 25;
	SolveQuad(a, b, c, Root1Real, Root1imag, Root2Real, Root2imag);
	if (IsNotEqual(Root1Real, 0) || IsNotEqual(Root1imag, -5) || IsNotEqual(Root2Real, 0) || IsNotEqual(Root2imag,5))
		return false;
	a = 3; b = 2; c = 5;
	SolveQuad(a, b, c, Root1Real, Root1imag, Root2Real, Root2imag);
	if (IsNotEqual(Root1Real, -1.0 / 3.0) || IsNotEqual(Root1imag, -sqrt(14.0) / 3.0) || IsNotEqual(Root2Real,-1.0 / 3.0)
		|| IsNotEqual(Root2imag,sqrt(14.0) / 3.0))
		return false;
	a = 0; b = 2; c = 5;
	SolveQuad(a, b, c, Root1Real, Root1imag, Root2Real, Root2imag);
	if (IsNotEqual(Root1Real, -5.0 / 2.0) || IsNotEqual(Root1imag, 0) || IsNotEqual(Root2Real, -5.0 / 2.0) || IsNotEqual(Root2imag,0))
		return false;
	NumTests += 1;
	return true;
}
#endif