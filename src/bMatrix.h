#ifndef MatrixH
#define MatrixH
#include <vector>
#include <fstream>
using namespace std;
class bMatrix
{
	int dim1;//number of rows
	int dim2;//number of columns
	static double ErrorTolerance;
	vector<double> mat;
protected:
	void ConstructionBody(int M, int N);
	void CopyBody(const bMatrix& A);
	void DeleteBody();
	double& Set(int i, int j);
	const double& Get(int i, int j) const;
	bMatrix Add(const bMatrix& B) const;
	bMatrix Multiply(double a) const;
	bMatrix Multiply(const bMatrix& B) const;
	bMatrix Subtract(const bMatrix& B) const;
	bMatrix Divide(double a) const;
	bool IsEqual(const bMatrix& B) const;
public:
	bMatrix(int M = 1, int N = 1);
	bMatrix(const bMatrix& A);
	bMatrix& operator=(const bMatrix& A);
	double& operator()(int i, int j);
	const double& operator()(int i, int j) const;
	int GetDim1() const;
	int GetDim2() const;
	bMatrix SubMatrix(int rowNum, int colNum, int drow, int dcol) const;
	bMatrix insertSubMatrix(int rowNum, int colNum, const bMatrix& SubMatrix) const;
	bMatrix insertSubMatrix(int rowNum, int colNum, int drow, int dcol) const;
	bMatrix operator+(const bMatrix& B) const;
	bMatrix operator*(double a) const;
	bMatrix operator*(const bMatrix& B) const;
	bMatrix operator-(const bMatrix& B) const;
	bMatrix operator/(double a) const;
	bool operator==(const bMatrix& B) const;
	bool operator!=(const bMatrix& B) const;
	bMatrix ChangeDim1(int newDim1) const;
	bMatrix ChangeDim2(int newDim2) const;
	bMatrix Transpose() const;
	void SwitchRows(int i, int j);
	void print(ofstream& f) const;
	double Min(int* iMin = 0, int* jMin = 0) const;
	double Max(int* iMax = 0, int* jMax = 0) const;
	void SetRowZero(int i);
	static void SetErrorTolerance(double errTol);
	static double GetErrorTolerance();
};
#endif
