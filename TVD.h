#pragma once
class TVD
{
public:
	enum MODE { NO_TVD, QUICK, VAN_LEER, VAN_ALBADA } mode;
	TVD(TVD::MODE mode_ = TVD::MODE::NO_TVD);
	double operator()(double r) const;
};

