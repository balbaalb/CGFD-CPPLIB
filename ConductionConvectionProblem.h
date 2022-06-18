#pragma once
#include "MathUtils.h"
#include "Vector3D.h"
#include <functional>
typedef function<double(const Vector3D& P)> VelocityField[2];
class ConductionConvectionProblem
{
public:
	virtual double T(double x, double y, double z = 0) const = 0;
	virtual double T(const Vector3D& P) const = 0;
	virtual double vx(double x, double y, double z = 0) const = 0;
	virtual double vx(const Vector3D& P) const = 0;
	virtual double vy(double x, double y, double z = 0) const = 0;
	virtual double vy(const Vector3D& P) const = 0;
};

class ConductionExample : public ConductionConvectionProblem
{
	double Lx{1}, Ly{1}, T0{0}, T1{1}, A{pi};
public:
	ConductionExample() = default;
	ConductionExample(double Lx_, double Ly_, double T0_, double T1_);
	virtual double T(double x, double y, double z = 0) const;
	virtual double T(const Vector3D& P) const;
	virtual double vx(double x, double y, double z = 0) const;
	virtual double vx(const Vector3D& P) const;
	virtual double vy(double x, double y, double z = 0) const;
	virtual double vy(const Vector3D& P) const;
	double dT_dx(double x, double y, double z = 0) const;
	double dT_dx(const Vector3D& P) const;
	double dT_dy(double x, double y, double z = 0) const;
	double dT_dy(const Vector3D& P) const;
};

class convectionConstVelocity : public ConductionConvectionProblem
{
	double T0{0}, T1{1}, constVx{0}, constVy{0}, Lx{1}, Ly{1};
	double kapa{1};//conductivity / density
public:
	convectionConstVelocity() = default;
	convectionConstVelocity(double T0_, double T1_, double vx_, double vy_, double Lx_, double Ly_, double kinematic_diffusivity);
	virtual double T(double x, double y, double z) const;
	virtual double T(const Vector3D& P) const;
	virtual double vx(double x, double y, double z = 0) const;
	virtual double vx(const Vector3D& P) const;
	virtual double vy(double x, double y, double z = 0) const;
	virtual double vy(const Vector3D& P) const;
};

class BaligaPatankarExample1 : public ConductionConvectionProblem
{
protected:
	double Pe{1}, x0{ sqrt(2.0) / 2.0}, y0{ sqrt(2.0) / 2.0 };
public:
	BaligaPatankarExample1() = default;
	BaligaPatankarExample1(double Pe_);
	BaligaPatankarExample1(double Pe_, double x0_, double y0_);
	virtual double T(double x, double y, double z = 0) const;
	virtual double T(const Vector3D& P) const;
	virtual double vx(double x, double y, double z = 0) const;
	virtual double vx(const Vector3D& P) const;
	virtual double vy(double x, double y, double z = 0) const;
	virtual double vy(const Vector3D& P) const;
};

class BaligaPatankarExample1_case2 : public BaligaPatankarExample1
{
public:
	BaligaPatankarExample1_case2() = default;
	BaligaPatankarExample1_case2(double Pe_);
	virtual double T(double x, double y, double z = 0) const;
	virtual double T(const Vector3D& P) const;
	double variableConductivity(double x, double y, double z = 0) const;
	double variableConductivity(const Vector3D& P) const;
};

class BaligaPatankarExample2 : public ConductionConvectionProblem
{
public:
	double Lx{10}, Ly{10}, alpha{0}, v0{1}, T0{0}, T1{1};
	bool averageOnBoundary{true};
	virtual double T(double x, double y, double z = 0) const;
	virtual double T(const Vector3D& P) const;
	virtual double vx(double x, double y, double z = 0) const;
	virtual double vx(const Vector3D& P) const;
	virtual double vy(double x, double y, double z = 0) const;
	virtual double vy(const Vector3D& P) const;
};

