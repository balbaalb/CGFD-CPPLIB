#ifndef TriangulatioQMH
#define TriangulatioQMH
#include <functional>
#include "Triangulation.h"
class TriangulationQM : public Triangulation
{
public:
	TriangulationQM(Triangulation& T);
	double GetMaxShewchuk_T2_1(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 2
	double GetMaxShewchuk_T2_2(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 2
	double GetMaxShewchuk_T2_3(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 2
	double GetMaxShewchuk_T2_4(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 2
	double GetMaxShewchuk_T4_1(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T4_2(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T4_3(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T4_4(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T4_5(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T4_6(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T4_7(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T4_8(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T6_1(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_2(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_3(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_4(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_5(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_6(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_7(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_7a(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 4
	double GetMaxShewchuk_T6_8(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_9(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_10(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_11(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double GetMaxShewchuk_T6_12(double& minVal, double& piVal, double& avgVal);//Ref Shewchuk 2002 Table 6
	double MaxVoronoiDistance(Vertex* v);
	double MaxVoronoiDistance(double* min_of_max = 0);//Max distance of a point in a Voronoi cell from the Voronoi site.
	double GetPointDistributionMeasure();//Ref: D2R33 Nguyen et al. 2009
	double GetPointDistributionRatio();//Ref: D2R33 Nguyen et al. 2009
	double GetRegulatoryMeasure();//Ref: D2R33 Nguyen et al. 2009
	double GetVoronoiArea(Vertex* v);
	double GetCellVolumeDeviation();//Ref: D2R33 Nguyen et al. 2009
	double Get2ndMomentTrace(Vertex* v);//based on the geometric (not Voronoi) region around each vertex
	double Get2ndMomentDeviatoricDet(Vertex* v);//based on the geometric (not Voronoi) region around each vertex
	double Get2ndMomentTraceMaxDev();//Ref: D2R33 Nguyen et al. 2009, but based on geoArea
	double MaxGeoAreaDistance(Vertex* v);
	double MaxGeoAreaDistance(double* min_of_max = 0);//Max distance of a point in a Geometric Area around vertex from the Vertex.
	double GetPointDistributionMeasure_GeoAreaBased();//Ref: D2R33 Nguyen et al. 2009, but based on geoArea
	double GetPointDistributionRatio_GeoAreaBased();//Ref: D2R33 Nguyen et al. 2009, but based on geoArea
	double GetGeoArea(Vertex* v);
	double GetCellVolumeDeviation_GeoAreaBased();//Ref: D2R33 Nguyen et al. 2009, but based on geoArea
	double GetMaxArea();
	double GetMinArea();
	double GetMaxAreaMeasure();//Ref: D2R33 Nguyen 2009
	double GetCircleRatioMeasure();//Ref: D2R33 Nguyen 2009
	double GetCircumRadDeviationMeasure(); //Ref: D2R33 Nguyen 2009
	double GetStretchFactor(Vertex* v, shared_ptr<QuadEdgeIndex> ind = 0);//Ref: Cui, S., Kanj, I. A., & Xia, G. 2011 Comp. Geo. 44, pp 104-109  
	double GetStretchFactor();//Ref: D2R61 Cui, S., Kanj, I. A., & Xia, G. 2011 Comp. Geo. 44, pp 104-109  
	double GetMinCircumRadius();
	double GetMaxCircumRadius();
	double GetMinInRadius();
	double GetMaxInRadius();
	void GetSizeRelations(Face* f, double& Lmin, double& Lmed, double& Lmax);
	double GetMinContainmentRadius(Face* f);
	double EdgeOrientationAverage();
	double EdgeOrientationStdDeviation();
	double EdgeOrientationSkewness();
	double EdgeOrientationKurtosis();
	double EdgeOrientationCoverage();//0 is the edge orientations are distributted uniformlly from 0 to pi.
	double EdgeOrientationCoverage_10boxes();//0 is the edge orientations are distributted uniformlly from 0 to pi.
	double GetAvgVelDotEdge(function<double(const Vector3D& P)> V[2], double& minVal, double& maxVal, double& piVal);
	static bool ShiftPointsRandomly(Triangulation* T);//for now, only works for square domains
};
bool tester_TriangulationQM(int& NumTests);
bool tester_TriangulationQM_1(int& NumTests);//testing shape factor
bool tester_TriangulationQM_2(int& NumTests);//testing ShiftPoints() and some other functions
#endif
