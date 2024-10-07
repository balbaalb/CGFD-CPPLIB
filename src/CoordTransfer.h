#ifndef CoordTransferH
#define CoordTransferH
#include "Vector3D.h"
#include "bMatrix.h"
class CoordTransfer//tested only for 2D conditions
{
	bMatrix* M;
	bMatrix* RevM;
	Vector3D X0;
	void CopyBody(const CoordTransfer& rhs);
	void DeleteBody();
	void UpdateReverseMatrix();
public:
	CoordTransfer();
	CoordTransfer(const CoordTransfer& rhs);
	~CoordTransfer();
	CoordTransfer& operator=(const CoordTransfer& rhs);
	void Translate(double x, double y = 0, double z = 0);//Needs Testing
	void Translate(const Vector3D& C);//Needs Testing
	void Rotate(double theta);//Needs Testing
	void Rotate(double dx, double dy);//Needs Testing
	void SetRotationMatrix(const bMatrix& A);
	void OriginalToNew(double x, double y, double& x1, double& y1) const;
	void NewToOriginal(double x1, double y1, double& x, double& y) const;
	Vector3D OriginalToNew(const Vector3D& p0) const;
	Vector3D NewToOriginal(const Vector3D& p1) const;
};
#endif