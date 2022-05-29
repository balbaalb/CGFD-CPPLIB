#include "bMatrix.h"
#include "MathUtils.h"
//=======================================================================================
double bMatrix::ErrorTolerance = 1.0e-10;
void bMatrix::ConstructionBody(int M, int N)
{
	this->dim1 = M;
	this->dim2 = N;
	this->mat.resize(M * N, 0);
}
void bMatrix::CopyBody(const bMatrix& A)
{
	this->mat.resize(A.mat.size(), 0);
	for (int i = 0; i < A.mat.size(); ++i)
	{
		this->mat[i] = A.mat[i];
	}
	this->dim1 = A.dim1;
	this->dim2 = A.dim2;
}
void bMatrix::DeleteBody()
{
	this->mat.resize(0);
	this->dim1 = 0;
	this->dim2 = 0;
}
bMatrix::bMatrix(int M, int N)
{
	this->ConstructionBody(M, N);
}
bMatrix::bMatrix(const bMatrix& A)
{
	this->CopyBody(A);
}
void bMatrix::operator=(const bMatrix& A)
{
	this->DeleteBody();
	this->CopyBody(A);
}
int bMatrix::GetDim1() const
{
	return this->dim1;
}
int bMatrix::GetDim2() const
{
	return this->dim2;
}
double& bMatrix::Set(int i, int j)
{
	if (i < 0 || j < 0 || i >= this->dim1 || j >= this->dim2)
	{
		throw "bMatrix::operator(): Index i or j is out of bounds";
	}
	return this->mat[i * this->dim2 + j];
}
const double& bMatrix::Get(int i, int j) const
{
	if (i < 0 || j < 0 || i >= this->dim1 || j >= this->dim2)
	{
		throw "bMatrix::operator(): Index i or j is out of bounds";
	}
	return this->mat[i * this->dim2 + j];
}
double& bMatrix::operator()(int i, int j)
{
	return this->Set(i, j);
}
const double& bMatrix::operator()(int i, int j) const
{
	return this->Get(i, j);
}
bMatrix bMatrix::Add(const bMatrix& B) const
{
	if (B.dim1 != this->dim1 || B.dim2 != this->dim2)
		throw "bMatrix::Add(): Input size is incorrect.";
	bMatrix C(this->dim1, this->dim2);
	for (int i = 0; i < this->dim1; ++i)
	{
		for (int j = 0; j < this->dim2; ++j)
		{
			C(i, j) = (*this)(i, j) + B(i, j);
		}
	}
	return C;
}
bMatrix bMatrix::Multiply(double a) const
{
	bMatrix C(this->dim1, this->dim2);
	for (int i = 0; i < this->dim1; ++i)
	{
		for (int j = 0; j < this->dim2; ++j)
		{
			C(i, j) = (*this)(i, j) * a;
		}
	}
	return C;
}
bMatrix bMatrix::Multiply(const bMatrix& B) const
{
	if (B.dim1 != this->dim2)
		throw "bMatrix::Multiply(Matrix&): Input size is incorrect.";
	bMatrix C(this->dim1, B.dim2);
	for (int i = 0; i < this->dim1; ++i)
	{
		for (int j = 0; j < B.dim2; ++j)
		{
			for (int k = 0; k < this->dim2; ++k)
			{
				C(i, j) += (*this)(i, k) * B(k, j);
			}
		}
	}
	return C;
}
bMatrix bMatrix::Subtract(const bMatrix& B) const
{
	return (this->Add(B.Multiply(-1)));
}
bMatrix bMatrix::Divide(double a) const
{
	if (a == 0)
		throw "bMatrix::Divide(double a) : Division over zero!";
	return this->Multiply(1.0 / a);
}
bool bMatrix::IsEqual(const bMatrix& B) const
{
	if (B.dim1 != this->dim1 || B.dim2 != this->dim2)
		return false;
	for (int i = 0; i < this->dim1; ++i)
	{
		for (int j = 0; j < this->dim2; ++j)
		{
			if (fabs(B(i, j) - (*this)(i, j)) > bMatrix::ErrorTolerance)
				return false;
		}
	}
	return true;
}
bMatrix bMatrix::SubMatrix(int rowNum, int colNum, int drow, int dcol) const
{
	if (rowNum < 0 || rowNum >= this->dim1 || colNum < 0 || colNum >= this->dim2 ||
		drow < 0 || drow >= this->dim1 || dcol < 0 || dcol >= this->dim2)
		throw "bMatrix::SubMatrix(): Inputs are incorrect";
	bMatrix C(drow, dcol);
	for (int i = rowNum; i < rowNum + drow; ++i)
	{
		for (int j = colNum; j < colNum + dcol; ++j)
		{
			C(i, j) = (*this)(i, j);
		}
	}
	return C;
}
bMatrix bMatrix::insertSubMatrix(int rowNum, int colNum, const bMatrix& SubMatrix) const
{
	int drow = SubMatrix.dim1;
	int dcol = SubMatrix.dim2;
	bMatrix C = this->insertSubMatrix(rowNum, colNum, drow, dcol);
	for (int i = rowNum; i < rowNum + drow; ++i)
	{
		for (int j = colNum; j < colNum + dcol; ++j)
		{
			C(i, j) = SubMatrix(i - rowNum, j - colNum);
		}
	}
	return C;
}
bMatrix bMatrix::insertSubMatrix(int rowNum, int colNum, int drow, int dcol) const
{
	if (rowNum < 0 || rowNum >= this->dim1 || colNum < 0 || colNum >= this->dim2)
		throw "bMatrix::insertSubMatrix(): Inputs are incorrect";
	bMatrix C(this->dim1 + drow, this->dim2 + dcol);
	for (int i = 0; i < this->dim1; ++i)
	{
		for (int j = 0; j < this->dim2; ++j)
		{
			int ii = (i < rowNum) ? i : i + drow;
			int jj = (j < colNum) ? j : j + dcol;
			C(ii, jj) = (*this)(i, j);
		}
	}
	return C;
}
bMatrix bMatrix::operator+(const bMatrix& B) const
{
	return this->Add(B);
}
bMatrix bMatrix::operator*(double a) const
{
	return this->Multiply(a);
}
bMatrix bMatrix::operator*(const bMatrix& B) const
{
	return this->Multiply(B);
}
bMatrix bMatrix::operator-(const bMatrix& B) const
{
	return this->Subtract(B);
}
bMatrix bMatrix::operator/(double a) const
{
	return this->Divide(a);
}
bool bMatrix::operator==(const bMatrix& B) const
{
	return this->IsEqual(B);
}
bool bMatrix::operator!=(const bMatrix& B) const
{
	return !this->IsEqual(B);
}
bMatrix bMatrix::ChangeDim1(int newDim1) const
{
	if (newDim1 <= 0)
		throw "bMatrix::ChangeDim1(): Inputs are incorrect";
	bMatrix C(newDim1, this->dim2);
	for (int i = 0; i < C.dim1; ++i)
	{
		if (i < this->dim1)
		{
			for (int j = 0; j < C.dim2; ++j)
			{
				C(i, j) = (*this)(i, j);
			}
		}
	}
	return C;
}
bMatrix bMatrix::ChangeDim2(int newDim2) const
{
	if (newDim2 <= 0)
		throw "bMatrix::ChangeDim2(): Inputs are incorrect";
	bMatrix C(this->dim1, newDim2);
	for (int i = 0; i < C.dim1; ++i)
	{
		for (int j = 0; j < C.dim2; ++j)
		{
			if (j < this->dim2)
			{
				C(i, j) = (*this)(i, j);
			}
		}
	}
	return C;
}
bMatrix bMatrix::Transpose() const
{
	bMatrix C(this->dim2, this->dim1);
	for (int i = 0; i < this->dim1; ++i)
	{
		for (int j = 0; j < this->dim2; ++j)
		{
			C(j, i) = (*this)(i, j);
		}
	}
	return C;
}
void bMatrix::SwitchRows(int i, int j)
{
	if (i < 0 || j < 0 || i >= this->GetDim1() || j > this->GetDim1())
		throw "Matrix::SwitchRows(i,j): Inputs are incorrect";
	for (int k = 0; k < this->GetDim2(); ++k)
	{
		swap((*this)(i, k), (*this)(j, k));
	}
}
double bMatrix::Min(int* iMin, int* jMin) const
{
	double MinValue = (*this)(0, 0);
	if (iMin) *iMin = 0;
	if (jMin) *jMin = 0;
	for (int i = 0; i < this->dim1; ++i) {
		for (int j = 0; j < this->dim2; ++j) {
			if ((*this)(i, j) < MinValue){
				MinValue = (*this)(i, j);
				if (iMin) *iMin = i;
				if (jMin) *jMin = j;
			}
		}
	}
	return MinValue;
}
double bMatrix::Max(int* iMax, int* jMax) const
{
	double MaxValue = (*this)(0, 0);
	if (iMax) *iMax = 0;
	if (jMax) *jMax = 0;
	for (int i = 0; i < this->dim1; ++i) {
		for (int j = 0; j < this->dim2; ++j) {
			if ((*this)(i, j) > MaxValue) {
				MaxValue = (*this)(i, j);
				if (iMax) *iMax = i;
				if (jMax) *jMax = j;
			}
		}
	}
	return MaxValue;
}
void bMatrix::SetRowZero(int i)
{
	for (int j = 0; j < this->dim2; ++j)
		(*this)(i, j) = 0;
}
void bMatrix::print(ofstream& f) const
{
	if (f.is_open())
	{
		for (int i = 0; i < this->dim1; ++i) {
			f << endl << "[";
			for (int j = 0; j < this->dim2; ++j) {
				f << " " << (*this)(i, j) << " ";
			}
			f << "]";
		}
	}
}
void bMatrix::SetErrorTolerance(double errTol)
{
	bMatrix::ErrorTolerance = errTol;
}
double bMatrix::GetErrorTolerance()
{
	return bMatrix::ErrorTolerance;
}
bool tester_bMatrix(int& NumTests)
{
	bMatrix A(3, 5);
	A(0, 0) = 12; A(0, 1) = 13; A(0, 2) = -1; A(0, 3) = -5.5; A(0, 4) = -13;
	A(1, 0) = 9; A(1, 1) = 7; A(1, 2) = 6.5; A(1, 3) = -7.5; A(1, 4) = 5.3;
	A(2, 0) = 0; A(2, 1) = 15; A(2, 2) = -18; A(2, 3) = 19; A(2, 4) = -11.5;
	bMatrix AA(A);
	bMatrix AAA(10,10);
	AAA = AA;
	if (AAA.GetDim1() != 3)
		return false;
	if (AAA.GetDim2() != 5)
		return false;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			if (fabs(AAA(i, j) - A(i, j)) > 1.0e-15)
				return false;
		}
	}
	if (fabs(AAA(2, 3) - 19) > 1.0e-10)
		return false;
	bool isEqual = (AAA == A);
	if (!isEqual)
		return false;
	bool notEqual = (AAA != A);
	if (notEqual)
		return false;
	bMatrix B = A * 2.56;
	if (fabs(B(2, 3) - 19 * 2.56) > 1.0e-10)
		return false;
	bMatrix C = B + A;
	if (fabs(C(2, 3) - 19 * 3.56) > 1.0e-10)
		return false;
	bMatrix D = C - B;
	if (fabs(D(2, 3) - 19) > 1.0e-10)
		return false;
	B = B / 6.16;
	if (fabs(B(2, 3) - 19 * 2.56 / 6.16) > 1.0e-10)
		return false;
	bMatrix E(5, 2);
	E(0, 0) = 1.5; E(0, 1) = 14.5;
	E(1, 0) = -6; E(1, 1) = 8.5;
	E(2, 0) = 5.3; E(2, 1) = 1.8;
	E(3, 0) = 8.5; E(3, 1) = -9.7;
	E(4, 0) = 4.3; E(4, 1) = 11.5;
	bMatrix F = A * E;
	if (F.GetDim1() != 3)
		return false;
	if (F.GetDim2() != 2)
		return false;
	if (fabs(F(0, 0) + 167.95) > 1.0e-10)
		return false;
	if (fabs(F(0, 1) - 186.55) > 1.0e-10)
		return false;
	if (fabs(F(1, 0) + 35.01) > 1.0e-10)
		return false;
	if (fabs(F(1, 1) - 335.4) > 1.0e-10)
		return false;
	if (fabs(F(2, 0) + 73.35) > 1.0e-10)
		return false;
	if (fabs(F(2, 1) + 221.45) > 1.0e-10)
		return false;
	AAA = AAA.ChangeDim1(4);
	AAA(3, 0) = 10; AAA(3, 1) = 1.5; AAA(3, 2) = 9; AAA(3, 3) = 36; AAA(3, 4) = 12.6;
	AAA = AAA.ChangeDim2(6);
	if (AAA.GetDim1() != 4)
		return false;
	if (AAA.GetDim2() != 6)
		return false;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			if (fabs(AAA(i, j) - A(i, j)) > 1.0e-15)
				return false;
		}
	}
	if (fabs(AAA(3, 0) - 10) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 1) - 1.5) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 2) - 9) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 3) - 36) > 1.0e-15)
		return false;
	if (fabs(AAA(3, 4) - 12.6) > 1.0e-15)
		return false;
	for (int i = 0; i < 4; ++i)
	{
		if (fabs(AAA(i, 5)) > 1.0e-15)
			return false;
	}
	AAA(0, 5) = -1;
	AAA(1, 5) = -2;
	AAA(2, 5) = -3;
	AAA(3, 5) = -4;
	AAA = AAA.insertSubMatrix(1, 4, F);
	if (AAA.GetDim1() != 7)
		return false;
	if (AAA.GetDim2() != 8)
		return false;
	for (int j = 0; j < 4; ++j)
	{
		if (fabs(AAA(0, j) - A(0, j)) > 1.0e-15)
			return false;
	}
	if (fabs(AAA(0, 6) - A(0, 4)) > 1.0e-15)
		return false;
	if (fabs(AAA(0, 7) + 1) > 1.0e-15)
		return false;
	for (int i = 1; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (fabs(AAA(i, j)) > 1.0e-15)
				return false;
		}
	}
	for (int i = 1; i < 4; ++i)
	{
		for (int j = 4; j < 6; ++j)
		{
			if (fabs(AAA(i, j) - F(i - 1, j - 4) )> 1.0e-15)
				return false;
		}
	}
	for (int i = 1; i < 4; ++i)
	{
		for (int j = 6; j < 8; ++j)
		{
			if (fabs(AAA(i, j))> 1.0e-15)
				return false;
		}
	}
	for (int i = 4; i < 6; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (fabs(AAA(i, j) - A(i - 3,j))> 1.0e-15)
				return false;
		}
	}
	if (fabs(AAA(6, 0) - 10)> 1.0e-15)
			return false;
	if (fabs(AAA(6, 1) - 1.5)> 1.0e-15)
		return false;
	if (fabs(AAA(6, 2) - 9)> 1.0e-15)
		return false;
	if (fabs(AAA(6, 3) - 36)> 1.0e-15)
		return false;
	for (int i = 4; i < 7; ++i)
	{
		for (int j = 4; j < 6; ++j)
		{
			if (fabs(AAA(i, j))> 1.0e-15)
				return false;
		}
	}
	for (int i = 4; i < 6; ++i)
	{
		if (fabs(AAA(i, 6) - A(i - 3, 4))> 1.0e-15)
			return false;
	}
	if (fabs(AAA(6, 6) - 12.6)> 1.0e-15)
		return false;
	if (fabs(AAA(4, 7) + 2)> 1.0e-15)
		return false;
	if (fabs(AAA(5, 7) + 3)> 1.0e-15)
		return false;
	if (fabs(AAA(6, 7) + 4)> 1.0e-15)
		return false;
	bMatrix AT = A.Transpose();
	if (AT.GetDim1() != 5)
		return false;
	if (AT.GetDim2() != 3)
		return false;
	if (fabs(AT(2, 1) - 6.5)> 1.0e-15)
		return false;
	if (fabs(AT(4,2) - A(2,4))> 1.0e-15)
		return false;
	if (fabs(AT(1,1) -7)> 1.0e-15)
		return false;
	if (fabs(AT(3,2) - A(2,3))> 1.0e-15)
		return false;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 5; ++j)
		{
			if (fabs(AT(j, i) - A(i, j))> 1.0e-15)
				return false;
		}
	}
	bMatrix A2 = A.ChangeDim1(2);
	A2 = A2.ChangeDim2(4);
	if (A2.GetDim1() != 2)
		return false;
	if (A2.GetDim2() != 4)
		return false;
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (fabs(A2(i, j) - A(i, j)) > 1.0e-15)
				return false;
		}
	}
	A(0, 0) = 12; A(0, 1) = 13; A(0, 2) = -1; A(0, 3) = -5.5; A(0, 4) = -13;
	A(1, 0) = 9; A(1, 1) = 7; A(1, 2) = 6.5; A(1, 3) = -7.5; A(1, 4) = 5.3;
	A(2, 0) = 0; A(2, 1) = 15; A(2, 2) = -18; A(2, 3) = 19; A(2, 4) = -11.5;
	bMatrix A_orig = A;
	bMatrix A_swaped(3, 5);
	A_swaped(2, 0) = 12; A_swaped(2, 1) = 13; A_swaped(2, 2) = -1; A_swaped(2, 3) = -5.5; A_swaped(2, 4) = -13;
	A_swaped(1, 0) = 9; A_swaped(1, 1) = 7; A_swaped(1, 2) = 6.5; A_swaped(1, 3) = -7.5; A_swaped(1, 4) = 5.3;
	A_swaped(0, 0) = 0; A_swaped(0, 1) = 15; A_swaped(0, 2) = -18; A_swaped(0, 3) = 19; A_swaped(0, 4) = -11.5;
	A.SwitchRows(0, 2);
	if (A == A_orig)
		return false;
	if (A != A_swaped)
		return false;
	A.SwitchRows(0, 2);
	if (A != A_orig)
		return false;
	if (A == A_swaped)
		return false;
	NumTests += 1;
	return true;
}