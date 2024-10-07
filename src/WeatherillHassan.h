#ifndef WeatherillHassanH
#define WeatherillHassanH
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "Triangulation.h"
#include "Vector3D.h"
class Vertex;
using namespace std;
class WeatherillHassan
{
	double scale0;
	unordered_set<Vertex*> boundary;
	unordered_map<Vertex*, double> LenScale;
	int numInputPoints;
	unique_ptr<Triangulation> triang;
	WeatherillHassan();
	~WeatherillHassan();
	void UpdateBoundary();
	void UpdateLengthScales();
	int InsertVertices(vector<Vector3D>& points);
	void Build(const vector<Vector3D>& BoundaryPoints, int maxIter);
public:
	enum LENGTH_SCALE_DEFINITION
	{
		MAX_ADJ_EDGE_LENGTH,
		AVERAGE_ADJ_EDGE_LENGTHS
	};
	static LENGTH_SCALE_DEFINITION length_scale_definition;
	static double alpha;
	static double beta;
	static Triangulation Tessellate(const vector<Vector3D>& BoundaryPoints, int maxIter);
	/*
	Issues: 
	
	1. Read Weatherill & Hassan paper and see when the scale is calculated, is it based only on boundary edges or
	is it based on the all edges (internal and boundary) in the first trianguilation that will include only boundary points.
	The current code does the latter but the results seems to be the same as the boundary edges are the smaller edges compared
	to boundary edges.
	
	2. It seems in the Weatherill & Hassan paper, the scales of initial points are calculated and for the later points that are added,
	their scale is calculated based only on initial points and not all points inserted so far. They suggest to calculate this scale based on 
	interpolation of initial nodes scales. This effort is not done here and the it is assumed that boundary points are distributed
	uniformly and thus the min scale of boundary points is taken as a constant length scale for all later nodes being inserted.
	*/
};
#endif