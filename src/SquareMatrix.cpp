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
double SquareMatrix::GetGaussSeidelTolerance() const
{
	return this->GaussSeidelTolerance;
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
bool SquareMatrix::CanUse_GaussSeidel(bool fixRoundOffError)
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
		if (!acceptable)
			return false;
		if (fabs((*this)(i, i)) < sum)
		{
			double aP = fabs((*this)(i, i));
			double delta = sum - aP;
			if (delta < 0.0001 * aP && fixRoundOffError)
			{
				aP = sum * (1.0001);
			}
			else 
			{
				return false;
			}
		}
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