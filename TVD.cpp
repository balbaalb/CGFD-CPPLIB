#include "TVD.h"
#include <cmath>
TVD::TVD(TVD::MODE mode_) : mode(mode_)
{
}
double TVD::operator()(double r) const
{
	switch (mode)
	{
	case TVD::NO_TVD:
		return 0;
	case TVD::QUICK:
		return fmax(0, fmin(2.0, fmin(2.0 * r, (3.0 + r) / 4.0)));
	case TVD::VAN_LEER:
		return (r + fabs(r)) / (1.0 + r);
	case TVD::VAN_ALBADA:
		return (r + r * r) / (1.0 + r * r);
	default:
		return 0;
	}
}
