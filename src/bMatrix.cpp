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
bMatrix& bMatrix::operator=(const bMatrix& A)
{
	this->DeleteBody();
	this->CopyBody(A);
	return *this;
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