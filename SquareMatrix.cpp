#include <math.h>
#include "SquareMatrix.h"
SquareMatrix::SquareMatrix(int N) : bMatrix(N,N)
{

}
SquareMatrix::SquareMatrix(const SquareMatrix& rhs)
{
	this->bMatrix::CopyBody(rhs);
}
SquareMatrix& SquareMatrix::operator=(const SquareMatrix& rhs)
{
	this->bMatrix::DeleteBody();
	this->bMatrix::CopyBody(rhs);
	return *this;
}
int SquareMatrix::GetDim() const
{
	return this->GetDim1();
}
void SquareMatrix::SetGaussSeidelTolerance(double value)
{
	this->GaussSeidelTolerance = value;
}
double SquareMatrix::MakeUpperTriangular(Vector* b)
{
	bool doPivotting = true;
	double detCorrection = 1.0;
	int N = this->GetDim();
	Vector scale(N);
	for (unsigned i = 0; i < N; ++i)
	{
		scale(i) = fabs((*this)(i, 1));
		for (int j = 1; j < N; ++j)
		{
			if (fabs((*this)(i, j)) > scale(i))
			{
				scale(i) = fabs((*this)(i, j));
			}
		}
	}
	//Triangulation
	for (unsigned i = 0; i < N - 1; ++i)
	{
		//Finding the best Pivot: find a row bewteen rows i to N that has the maximum value of |A(ii,i)|/scale[ii], i<=ii<=N.
		//IF row p has this maximum swap rows p and i
		if (scale(i) > 1e-15)
		{
			if (doPivotting)
			{
				double scaledPivot = fabs((*this)(i, i)) / scale(i);
				int p = i;
				for (unsigned ii = i + 1; ii < N; ++ii)
				{
					if (fabs((*this)(ii, i)) / scale(ii) > scaledPivot)
					{
						scaledPivot = fabs((*this)(ii, i)) / scale(ii);
						p = ii;
					}
				}
				if (p != i)
				{
					this->SwitchRows(i, p);
					detCorrection *= -1;
					if (b)
					{
						b->SwitchRows(i, p);
					}
					scale.SwitchRows(i, p);
				}
			}
			//End of Finding the best Pivot
			for (unsigned ii = i + 1; ii < N; ++ii)
			{
				if (fabs((*this)(i, i)) > 1e-20)
				{
					double Coeff = -(*this)(ii, i) / (*this)(i, i);
					for (unsigned j = 0; j < N; ++j)
					{
						(*this)(ii, j) += Coeff * (*this)(i, j);
					}
					if (b)
					{
						(*b)(ii) += Coeff * (*b)(i);
					}
				}
			}
		}
	}
	return detCorrection;
}
Vector SquareMatrix::Solve_GaussEliniation(const Vector& B) const
{
	SquareMatrix K = (*this);
	Vector F = B;
	Vector X(this->GetDim());
	for (unsigned i = 0; i < K.GetDim(); ++i)
	{
		bool acceptable = false;
		for (unsigned j = 0; j < K.GetDim(); ++j)
		{
			if (fabs(K(i, j)) > 1.0e-10)
			{
				acceptable = true;
				break;
			}
		}
		if (!acceptable)
		{
			K(i, i) = 1;
			F(i) = 0;
		}
	}
	K.MakeUpperTriangular(&F);
	for (int i = K.GetDim() - 1; i >= 0; --i)
	{
		double RHS = F(i);
		for (int j = K.GetDim() - 1; j > i; --j)
		{
			RHS -= K(i, j) * X(j);
		}
		X(i) = RHS / K(i, i);
	}
	return X;
}
Vector SquareMatrix::Solve_GaussSeidel(const Vector& B) const
{
	double alpha = 0.5;
	Vector convergence(this->GetDim());
	int counter = 0;
	Vector X(this->GetDim());
	bool can_use_GaussSeidel = true;
	do
	{
		Vector X_old = X;
		for (int i = 0; i < this->GetDim(); ++i)
		{
			double sum = B(i);
			for (int j = 0; j < this->GetDim(); ++j)
			{
				sum -= (i != j) ? (*this)(i, j) * X(j) : 0;
			}
			X(i) = sum / (*this)(i, i);
		}
		X = X * alpha + X_old * (1.0 - alpha);
		convergence = X - X_old;
		++counter;
	} while (counter < 1000 && (counter < 1 || convergence.abs() > this->GaussSeidelTolerance));
	return X;
}
bool SquareMatrix::CanUse_GaussSeidel()
{
	for (int i = 0; i < this->GetDim(); ++i)
	{
		bool acceptable = false;
		double sum = 0;
		for (int j = 0; j < this->GetDim(); ++j)
		{
			sum += (i != j) ? fabs((*this)(i, j)) : 0;
			if (fabs((*this)(i, j)) > 1.0e-10)
			{
				acceptable = true;
			}
		}
		if (!acceptable || fabs((*this)(i, i)) <= sum)
			return false;
	}
	return true;
}
Vector SquareMatrix::Solve(const Vector& B, METHOD method) const
{
	if (B.GetDim1() != this->GetDim())
		throw "SquareMatrix::Solve()";
	if (method == GAUSS_ELIMINATION)
		return this->Solve_GaussEliniation(B);
	else
		return this->Solve_GaussSeidel(B);
}
bool tester_SquareMatrix(int& NumTests)
{
	SquareMatrix A(2);
	A(0, 0) = 11.5;	A(0, 1) = -2.35;
	A(1, 0) = 17;
	SquareMatrix A2(A), A3;
	A3 = A2;
	if (A3.GetDim1() != 2 || A3.GetDim1() != 2)
		return false;
	if (A3 != A)
		return false;
	if (!tester_SquareMatrix_1(NumTests))
		return false;
	if (!tester_SquareMatrix_2(NumTests))
		return false;
	if (!tester_SquareMatrix_3(NumTests))
		return false;
	if (!tester_SquareMatrix_4(NumTests))
		return false;
	if (!tester_SquareMatrix_5(NumTests))
		return false;
	++NumTests;
	return true;
}
bool tester_SquareMatrix_1(int& NumTests)
{
	SquareMatrix A(5);
	A(1, 1) = 4; 	A(1, 2) = 5; 	A(1, 3) = 12; 	A(1, 4) = 13.5;
	A(1, 1) = 0.2; 	A(1, 2) = 4;	A(2, 3) = 0.5; 	A(2, 4) = -17;
	A(3, 1) = 6;  	A(3, 2) = 14;	A(3, 3) = 13; 	A(3, 4) = -18;
	A(4, 1) = 5; 	A(4, 2) = 7; 	A(4, 3) = -9;	A(4, 4) = -9;
	Vector b(5);
	b(1) = 11;
	b(2) = -3;
	b(3) = 0;
	b(4) = 5;
	Vector X = A.Solve(b);
	Vector bb = A * X;
	if (bb != b)
		return false;
	for (int i = 0; i < 5; ++i){
		if (fabs(bb(i) - b(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_2(int& NumTests)
{
	//solution of A*xx=yy : Needs correct pivoting
	SquareMatrix A(4);
	A(1, 1) = 1.0; A(1, 2) = 0.0,	A(1, 3) = 0.0;
	A(2, 1) = 1.0; A(2, 2) = 0.0;	A(2, 3) = 1.0;
	A(3, 1) = 1.0; A(3, 2) = 1.0;	A(3, 3) = 0.0;
	Vector x(4); 
	x(1) = 1.0;
	x(2) = 0.5;
	x(3) = -2.0;
	Vector yy = A * x;
	Vector xx = A.Solve(yy);
	if (xx != x)
		return false;
	for (int i = 0; i < 4; ++i){
		if (fabs(xx(i) - x(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_3(int& NumTests)
{
	
	SquareMatrix A6(7);
	A6(1, 1) = 1; A6(1, 2) = 0; A6(1, 3) = 0; A6(1, 4) = 0; A6(1, 5) = 0; A6(1, 6) = 0;
	A6(2, 1) = 1; A6(2, 2) = 1; A6(2, 3) = 1; A6(2, 4) = 0; A6(2, 5) = 0; A6(2, 6) = 0;
	A6(3, 1) = 1; A6(3, 2) = 0; A6(3, 3) = 0; A6(3, 4) = 0; A6(3, 5) = 1; A6(3, 6) = 0;
	A6(4, 1) = 1; A6(4, 2) = 1; A6(4, 3) = 1; A6(4, 4) = 1; A6(4, 5) = 1; A6(4, 6) = 1;
	A6(5, 1) = 1; A6(5, 2) = 1; A6(5, 3) = 1; A6(5, 4) = 2; A6(5, 5) = 4; A6(5, 6) = 4;
	A6(6, 1) = 1; A6(6, 2) = 2; A6(6, 3) = 8; A6(6, 4) = 4; A6(6, 5) = 1; A6(6, 6) = 2;
	Vector x(7);
	x(1) = 0.5;
	x(2) = 1.6;
	x(3) = -1;
	x(4) = 2;
	x(5) = 3;
	x(6) = 1;
	Vector yy = A6 * x;
	Vector xx = A6.Solve(yy);
	if (xx != x)
		return false;
	for (int i = 0; i < 7; ++i){
		if (fabs(xx(i) - x(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_4(int& NumTests)
{
	
	SquareMatrix A(4);
	A(1, 1) = 5.0; A(1, 2) = 0.0, A(1, 3) = 0.0;
	A(2, 1) = -1.0; A(2, 2) = 0.0;	A(2, 3) = 16.0;
	A(3, 1) = 4.0; A(3, 2) = -6.0;	A(3, 3) = 0.0;
	Vector x(4);
	x(1) = 1.0;
	x(2) = 0.5;
	x(3) = -2.0;
	Vector yy = A * x;
	Vector xx = A.Solve(yy);
	if (xx != x)
		return false;
	for (int i = 0; i < 4; ++i){
		if (fabs(xx(i) - x(i)) > 1.0e-10)
			return false;
	}
	++NumTests;
	return true;
}
bool tester_SquareMatrix_5(int& NumTests)
{

	SquareMatrix A(4);
	Vector B(4), X(4);
	A(0, 0) = 10;	A(0, 1) = 2;	A(0, 2) = -1;	A(0, 3) = 5;	X(0) = 5;	B(0) = 64.5;
	A(1, 0) = 5;	A(1, 1) = 12;	A(1, 2) = 2;	A(1, 3) = 4;	X(1) = -11; B(1) = -61;
	A(2, 0) = -3;	A(2, 1) = 5;	A(2, 2) = 18; 	A(2, 3) = 8;	X(2) = 6;	B(2) = 106;
	A(3, 0) = -2;	A(3, 1) = 4;	A(3, 2) = 3;	A(3, 3) = 20;	X(3) = 8.5;	B(3) = 134;
	if (!A.CanUse_GaussSeidel())
		return false;
	Vector xx = A.Solve(B, SquareMatrix::METHOD::GAUSS_SEIDEL);
	for (int i = 0; i < 4; ++i) {
		if (fabs(xx(i) - X(i)) > 1.0e-5)
			return false;
	}
	A(2, 2) = 15.99;
	if (A.CanUse_GaussSeidel())
		return false;
	++NumTests;
	return true;
}