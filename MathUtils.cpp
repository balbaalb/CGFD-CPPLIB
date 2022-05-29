#include <cmath>
#include <MathUtils.h>
void SolveQuad(double a,double b, double c, double& Root1Real,double& Root1imag, double& Root2Real, double& Root2imag)
//solving a*x^2+b*x+c = 0.
{
	if (a == 0 && b == 0)
		throw "MathUtils::SolveQuad() 0*x^2+0*x+c = 0!";
	Root1Real = 0;
	Root1imag = 0;
	Root2Real = 0;
	Root2imag = 0;
	if (fabs(a) > 1.e-20)
	{
		double delta;
		double delta2 = b*b - 4.0*a*c;
		if (delta2 >= 0.)
		{
			delta = sqrt(delta2);
			Root1Real = (-b + delta) / 2. / a;
			Root2Real = (-b - delta) / 2. / a;
		}
		else
		{
			delta = sqrt(-delta2);
			Root1Real = -b / 2. / a;
			Root1imag = delta / 2. / a;
			Root2Real = Root1Real;
			Root2imag = -Root1imag;
		}
	}
	else if (fabs(b) > 1.e-20)
	{
		Root1Real = -c / b;
		Root2Real = Root1Real;
	}
	if ((Root1Real > Root2Real) || (Root1Real == Root2Real && Root1imag > Root2imag))
	{
		swap(Root1Real, Root2Real);
		swap(Root1imag, Root2imag);
	}
}
void swap(double& a, double& b)
{
	double c = a;
	a = b;
	b = c;
}
bool IsEqual(double a, double b)
{
	return (fabs(a - b) < 1.0e-10);
}
bool IsNotEqual(double a, double b)
{
	return (fabs(a - b) > 1.0e-10);
}
//=======================================================================================
bool tester_MathUtils(int& NumTests)
{
	double u = 1550.126;
	if (!IsEqual(u, 1550.126) || IsNotEqual(u, 1550.126))
		return false; 
	double v = 100, w = -1000;
	swap(v, w);
	if (IsNotEqual(v, -1000) || IsNotEqual(w, 100))
		return false;
	if (!IsEqual(v, -1000) || !IsEqual(w, 100))
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