#ifndef Test_TriangulationQMH
#define Test_TriangulationQMH
#include <sstream>
#include "../src/TriangulationQM.h"
#include "../src/QuadEdge.h"
#include "../src/DelaunayLifting.h"
#include "../src/VertexIterator.h"
#include "../src/Vertex.h"
bool tester_TriangulationQM_1(int& NumTests)
{//testing shape factor
	Triangulation T = Triangulation::OffDiagonalGrid(3, 6, 1.35, 2.7);
	double Alpha = T.GetShapeFactor();
	if (fabs(Alpha - sqrt(3.0) / 2.0) > 1.0e-5)
		return false;
	double maxAngle = T.GetMaxAngle();
	if (fabs(maxAngle - pi / 2.0) > 1.0e-10)
		return false;
	double minAngle = T.GetMinAngle();
	if (fabs(minAngle - pi / 4.0) > 1.0e-10)
		return false;
	int maxVertexDegree = T.GetMaxVertexDegree();
	if (maxVertexDegree != 6)
		return false;
	int minVertexDegree = T.GetMinVertexDegree();
	if (minVertexDegree != 2)
		return false;
	double Lx = 1.78;
	double Ly = 4.56;
	int Nx = 3;
	int Ny = 5;
	Triangulation T2i = Triangulation::OffDiagonalGrid(Nx, Ny, Lx, Ly);
	TriangulationQM T2(T2i);
	double hx = 1.78 / 3.0;/*0.593333...*/
	double hy = 4.56 / 5.0;/*0.912*/
	double hxy = sqrt(hx * hx + hy * hy);/*1.088*/
	double weight = 3.0 * 6.0 * hx + 4.0 * 5.0 * hy + 3.0 * 5.0 * hxy;
	if (fabs(T2.GetMaxEdgeLength() - hxy) > 1.0e-10)
		return false;
	if (fabs(T2.GetMinEdgeLength() - hx) > 1.0e-10)
		return false;
	if (fabs(T2.GetWeight() - weight) > 1.0e-10)
		return false;
	int ench = T2.GetEncrochness();
	if (ench != 2)
		return false;
	double skinnyness = T2.GetSkinnyness();
	if (fabs(skinnyness - hxy / 2.0 / hx) > 1.0e-10)
		return false;
	T2.UpdateArea();
	if (fabs(T2.GetArea() - 1.78 * 4.56) > 1.0e-10)
		return false;
	double hmin = 0;
	double hmax = T2.MaxVoronoiDistance(&hmin);
	if (fabs(hmax - hxy / 2.0) > 1.0e-10)
		return false;
	if (fabs(hmin - hxy / 2.0) > 1.0e-10)
		return false;
	if (fabs(T2.GetPointDistributionMeasure() - 0.5) > 0.05)
		return false;
	if (fabs(T2.GetPointDistributionRatio()) > 1e-10)
		return false;
	double Chi = sqrt(3.0) / 2.0 * hxy / hx;
	if (fabs(T2.GetRegulatoryMeasure() - Chi + 1.0) > 1.0e-10)
		return false;
	double nu = T2.GetCellVolumeDeviation();
	if (fabs(nu - 3) > 1.0e-10)
		return false;
	double tau = T2.Get2ndMomentTraceMaxDev();
	if (fabs(tau - 0.60651) > 0.01)
		return false;
	double dmin = 0;
	double dmax = T2.MaxGeoAreaDistance(&dmin);
	if (fabs(dmax - 0.64) > 0.1 || fabs(dmin - 0.456) > 0.1)
		return false;
	double nuPrime = T2.GetCellVolumeDeviation_GeoAreaBased();
	if (fabs(nuPrime - 5.0) > 0.01)
		return false;
	double alpha = T2.GetMaxAreaMeasure();
	if (fabs(alpha) > 1.0e-10)
		return false;
	double R = hxy / 2.0;
	double beta = atan(hy / hx);
	double r = hx * tan(beta / 2.0) / (1.0 + tan(beta / 2.0));
	double q = R / 2.0 / r - 1.0;
	if (fabs(q - T2.GetCircleRatioMeasure()) > 1.0e-10)
		return false;
	if (fabs(T2.GetCircumRadDeviationMeasure()) > 1.0e-4)
		return false;
	Triangulation T3i = Triangulation::OffDiagonalGrid(3, 1, Lx, Ly);
	TriangulationQM T3(T3i);
	double sf3 = T3.GetStretchFactor();
	double sf_expected = (Lx + Ly) / sqrt(Lx * Lx + Ly * Ly);
	if (fabs(sf3 - sf_expected) > 1.0e-5)
		return false;
	Triangulation T4i = Triangulation::OffDiagonalGrid(1, 2, 1.0, 2.0);
	TriangulationQM T4(T4i);
	double sf4 = T4.GetStretchFactor();
	if (fabs(sf4 - sqrt(2.0)) > 1.0e-5)
		return false;
	double sf2 = T2.GetStretchFactor();
	sf_expected = (3.0 * hx + 2.0 * hy) / sqrt(9.0 * hx * hx + 4.0 * hy * hy);
	if (fabs(sf2 - sf_expected) > 1.0e-5)
		return false;
	Face* f2 = T2.GetMesh2D()->GetFace(2);
	if (f2 == T2.GetBoundary())
		f2 = T2.GetMesh2D()->GetFace(1);
	double Lmin, Lmed, Lmax;
	T2.GetSizeRelations(f2, Lmin, Lmed, Lmax);
	if (fabs(Lmin - hx) > 1.0e-10 || fabs(Lmed - hy) > 1.0e-10 || fabs(Lmax - hxy) > 1.0e-10)
		return false;
	int h = T2.NumBoundaryPoints();
	if (h != 2 * (Nx + Ny))
		return false;
	double thetaAvg = T2.EdgeOrientationAverage();
	if (fabs(thetaAvg - 0.87407694750795) > 1.0e-5)
		return false;
	double sigma = T2.EdgeOrientationStdDeviation();
	if (fabs(sigma - 0.674769284104802) > 1.0e-5)
		return false;
	double skewness = T2.EdgeOrientationSkewness();
	if (fabs(skewness + 0.0986904482846126) > 1.0e-5)
		return false;
	double Kurtosis = T2.EdgeOrientationKurtosis();
	if (fabs(Kurtosis - 0.287217569357575) > 1.0e-5)
		return false;
	int chi1 = T2.EdgeOrientationCoverage();
	if (chi1 != 17 * 17 + 14 * 14 + 19 * 19 + 50)
		return false;
	int chi2 = T2.EdgeOrientationCoverage_10boxes();
	if (chi2 != 17 * 17 + 14 * 14 + 19 * 19 + 7)
		return false;
	++NumTests;
	return true;
}
bool tester_TriangulationQM_2(int& NumTests)
{
	Vector3D A, B(2), C(2, 2), D(0, 2), E(0.5, 0.5), F(1.5, 0.5), G(1.5, 1.5), H(0.5, 1.5);
	vector<Vector3D> input = { A,B,C,D,E,F,G,H };
	Triangulation T0 = DelaunayLifting::Triangulate(input);
	stringstream T0_ss;
	T0_ss << T0;
	Triangulation T;
	T0_ss >> T;
	T.Write("..\\Data\\test8a.bqe");
	vector<Vector3D> OrigPoints;
	T.GetPoints(OrigPoints);
	for (int i = 0; i < input.size(); ++i)
	{
		bool found_point_i = false;
		for (int j = 0; j < OrigPoints.size(); ++j)
		{
			if (input[i] == OrigPoints[j])
			{
				found_point_i = true;
				OrigPoints.erase(OrigPoints.begin() + j);
				break;
			}
		}
		if (!found_point_i)
			return false;
	}
	QuadEdge* qe0 = T.GetMesh2D();
	VertexIterator itv(qe0);
	Vertex* v = itv.Next();
	while (v)
	{
		Vector3D P = v->GetPoint();
		if (P == A || P == B || P == C || P == D)
		{
			if (!T.IsOnBoundary(v))
				return false;
		}
		else
		{
			if (T.IsOnBoundary(v))
				return false;
		}
		v = itv.Next();
	}
	TriangulationQM::ShiftPointsRandomly(&T);
	T.Write("..\\Data\\test8b.bqe");
	T = T0;
	TriangulationQM::ShiftPointsRandomly(&T);
	T.Write("..\\Data\\test8c.bqe");
	T = T0;
	TriangulationQM::ShiftPointsRandomly(&T);
	T.Write("..\\Data\\test8d.bqe");
	Triangulation T2 = Triangulation::OffDiagonalGrid(5, 4, 6.25, 7.11);
	T2.Write("..\\Data\\test9.bqe");
	TriangulationQM::ShiftPointsRandomly(&T2);
	T2.Write("..\\Data\\test9a.bqe");
	TriangulationQM::ShiftPointsRandomly(&T2);
	T2.Write("..\\Data\\test9b.bqe");
	TriangulationQM::ShiftPointsRandomly(&T2);
	T2.Write("..\\Data\\test9c.bqe");
	TriangulationQM::ShiftPointsRandomly(&T2);
	T2.Write("..\\Data\\test9d.bqe");
	++NumTests;
	return true;
}
bool tester_TriangulationQM(int& NumTests)
{
	if (!tester_TriangulationQM_1(NumTests))
		return false;
	if (!tester_TriangulationQM_2(NumTests))
		return false;
	NumTests += 1;
	return true;
}
#endif