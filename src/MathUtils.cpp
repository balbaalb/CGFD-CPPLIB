#include <cmath>
#include <algorithm>
#include "MathUtils.h"
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
		std::swap(Root1Real, Root2Real);
		std::swap(Root1imag, Root2imag);
	}
}
bool IsEqual(double a, double b)
{
	return (fabs(a - b) < 1.0e-10);
}
bool IsNotEqual(double a, double b)
{
	return (fabs(a - b) > 1.0e-10);
}