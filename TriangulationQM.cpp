#include <fstream>
#include <algorithm> 
#include <queue>
#include "Triangulation.h"
#include "Vector3D.h"
#include "Circle.h"
#include "Triangle.h"
#include "Vertex.h"
#include "Face.h"
#include "Edge.h"
#include "EdgeIterator.h"
#include "VertexIterator.h"
#include "FaceIterator.h"
#include "QuadEdge.h"
#include "bData.h"
#include "DelaunayLifting.h"
#include "bPolygon.h"
#include "MathUtils.h"
#include "QuadEdgeIndex.h"
#include "TriangulationQM.h"
TriangulationQM::TriangulationQM(Triangulation& T)
{
	this->Triangulation::CopyBody(T);
}
double TriangulationQM::GetMaxShewchuk_T2_1(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double r = this->GetIncircleRadius(f);
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Val = Lmax * Lmed * (Lmin + 4.0 * r) / 4.0 / A;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T2_2(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Val = 3.0 * Lmin * Lmed * Lmax / 4.0 / A;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T2_3(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetCircumCircleRadius(f);
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double amax = 2.0 * A / Lmin;
				double amed = 2.0 * A / Lmed;
				double Val = std::fmax(R, std::fmax(amax, sqrt(fabs(Lmax * Lmax - amed * amed))));
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T2_4(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double Val = 1.0 / Lmax / Lmax;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_1(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double rmc = this->GetMinContainmentRadius(f);
				double Val = 1.0 / rmc / rmc;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_2(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double rmc = this->GetMinContainmentRadius(f);
				double A = this->GetArea(f);
				double Val = A / rmc / rmc;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_3(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double r = this->GetIncircleRadius(f);
				double A = this->GetArea(f);
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double Val = A / Lmax / Lmed / (Lmin + 4.0 * r);
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_4(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Val = A / Lmin / Lmed / Lmax;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_5(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double r = this->GetIncircleRadius(f);
				double A = this->GetArea(f);
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double Val = A / pow(Lmax * Lmed * (Lmin + 4.0 * r), 2.0 / 3.0);
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_6(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Val = A / pow(Lmin * Lmed * Lmax, 2.0 / 3.0);
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_7(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Lrms = sqrt((Lmin * Lmin + Lmed * Lmed + Lmax * Lmax) / 3.0);
				double Val = A / (3.0 * Lrms * Lrms + sqrt(fabs(9.0 * Lrms * Lrms * Lrms * Lrms - 48.0 * A * A)));
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T4_8(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Lrms = sqrt((Lmin * Lmin + Lmed * Lmed + Lmax * Lmax) / 3.0);
				double Val = A / Lrms / Lrms;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_1(double& minVal, double& piVal, double& avgVal)
{
	double maxVal = this->GetMaxShewchuk_T4_5(minVal, piVal, avgVal);
	minVal *= 4.0 * pow((12.0 + 7.0 * sqrt(3.0)), 1.0 / 3.0) / 3.0;
	maxVal *= 4.0 * pow((12.0 + 7.0 * sqrt(3.0)), 1.0 / 3.0) / 3.0;
	avgVal *= 4.0 * pow((12.0 + 7.0 * sqrt(3.0)), 1.0 / 3.0) / 3.0;
	if (this->Mesh2D)
	{
		piVal = 1;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double r = this->GetIncircleRadius(f);
				double A = this->GetArea(f);
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double Val = A / pow(Lmax * Lmed * (Lmin + 4.0 * r), 2.0 / 3.0);
				piVal *= Val * 4.0 * pow((12.0 + 7.0 * sqrt(3.0)), 1.0 / 3.0) / 3.0;
			}
			f = itf.Next();
		}
	}
	return maxVal;
}
double TriangulationQM::GetMaxShewchuk_T6_2(double& minVal, double& piVal, double& avgVal)
{
	double maxVal = this->GetMaxShewchuk_T4_6(minVal, piVal, avgVal);
	minVal *= 4.0 / sqrt(3.0);
	maxVal *= 4.0 / sqrt(3.0);
	avgVal *= 4.0 / sqrt(3.0);
	if (this->Mesh2D)
	{
		piVal = 1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Val = A / pow(Lmin * Lmed * Lmax, 2.0 / 3.0);
				piVal *= Val * 4.0 / sqrt(3.0);
			}
			f = itf.Next();
		}
	}
	return maxVal;
}
double TriangulationQM::GetMaxShewchuk_T6_3(double& minVal, double& piVal, double& avgVal)
{
	double maxVal = this->GetMaxShewchuk_T4_7(minVal, piVal, avgVal);
	minVal *= 4.0 * sqrt(3.0);
	maxVal *= 4.0 * sqrt(3.0);
	avgVal *= 4.0 * sqrt(3.0);
	if (this->Mesh2D)
	{
		piVal = 1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Lrms = sqrt((Lmin * Lmin + Lmed * Lmed + Lmax * Lmax) / 3.0);
				double Val = A / (3.0 * Lrms * Lrms + sqrt(fabs(9.0 * Lrms * Lrms * Lrms * Lrms - 48.0 * A * A)));
				piVal *= Val * 4.0 * sqrt(3.0);
			}
			f = itf.Next();
		}
	}
	return maxVal;
}
double TriangulationQM::GetMaxShewchuk_T6_4(double& minVal, double& piVal, double& avgVal)
{
	double maxVal = this->GetMaxShewchuk_T4_8(minVal, piVal, avgVal);
	minVal *= 4.0 / sqrt(3.0);
	maxVal *= 4.0 / sqrt(3.0);
	avgVal *= 4.0 / sqrt(3.0);
	if (this->Mesh2D)
	{
		piVal = 1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Lrms = sqrt((Lmin * Lmin + Lmed * Lmed + Lmax * Lmax) / 3.0);
				double Val = A / Lrms / Lrms;
				piVal *= Val * 4.0 / sqrt(3.0);
			}
			f = itf.Next();
		}
	}
	return maxVal;
}
double TriangulationQM::GetMaxShewchuk_T6_5(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				Triangle t = this->GetTriangle(f);
				double thetaMin = t.GetMinAngle();
				double Val = 3.0 / pi * thetaMin;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_6(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				Triangle t = this->GetTriangle(f);
				double thetaMin = t.GetMinAngle();
				double Val = 2.0 / sqrt(3.0) * sin(thetaMin);
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_7(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Val = 4.0 / sqrt(3.0) * A / Lmax / Lmax;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_7a(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double A = this->GetArea(f);
				double Val = 4.0 / sqrt(3.0) * A / Lmin / Lmin;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_8(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double r = this->GetIncircleRadius(f);
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double Val = 2.0 * sqrt(3) * r / Lmax;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_9(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double r = this->GetIncircleRadius(f);
				double R = this->GetCircumCircleRadius(f);
				double Val = 2.0 * r / R;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_10(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetCircumCircleRadius(f);
				double A = this->GetArea(f);
				double Val = 4.0 * sqrt(3.0) / 9.0 * A / R / R;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_11(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double Val = Lmin / Lmax;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::GetMaxShewchuk_T6_12(double& minVal, double& piVal, double& avgVal)
{
	if (this->Mesh2D)
	{
		double maxVal = -1;
		minVal = -1;
		piVal = 1;
		avgVal = 0;
		int counter = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetCircumCircleRadius(f);
				double Lmin, Lmed, Lmax;
				this->GetSizeRelations(f, Lmin, Lmed, Lmax);
				double Val = 0.5 * Lmax / R;
				if (minVal < 0 || Val < minVal)
					minVal = Val;
				if (Val > maxVal)
					maxVal = Val;
				piVal *= Val;
				avgVal += Val;
				++counter;
			}
			f = itf.Next();
		}
		avgVal = counter ? avgVal / double(counter) : -1;
		return maxVal;
	}
	return -1;
}
double TriangulationQM::MaxVoronoiDistance(Vertex* v)
{
	if (v)
	{
		Vector3D P = v->GetPoint();
		double hk = 0;
		FaceIterator itf(v);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				VertexIterator itvv(f);
				Vertex* a = itvv.Next();
				Vertex* b = itvv.Next();
				Vertex* c = itvv.Next();
				if (b == v)
					std::swap(a, b);
				if (c == v)
					std::swap(a, c);
				Vector3D A = a->GetPoint();
				Vector3D B = b->GetPoint();
				Vector3D C = c->GetPoint();
				Triangle t(A, B, C);
				double theta = t.GetMaxAngle();
				if (theta <= pi / 2.0 + 1.0e-10)
				{
					Circle Cir = this->GetCircumCircle(f);
					Vector3D Center = Cir.GetCenter();
					Vector3D A_Center = Center - A;
					double h = A_Center.abs();
					if (h > hk)
						hk = h;
				}
				else {
					Vector3D MidAB = (A + B) / 2.0;
					Vector3D A_MidAB = MidAB - A;
					double hB = A_MidAB.abs();
					Vector3D MidAC = (A + C) / 2.0;
					Vector3D A_MidAC = MidAC - A;
					double hC = A_MidAC.abs();
					double BAC = Vector3D::GetAngle(B, A, C);
					if (BAC >= pi / 2.0)
					{
						if (std::fmax(hB, hC) > hk)
							hk = std::fmax(hB, hC);
					}
					else
					{
						if (std::fmin(hB, hC) > hk)
							hk = std::fmin(hB, hC);
					}
				}
			}
			f = itf.Next();
		}
		return hk;
	}
	return -1;
}
double TriangulationQM::MaxVoronoiDistance(double* min_of_max)
{
	if (this->Mesh2D)
	{
		int Nv = this->Mesh2D->NumVertices();
		double* max_h_Arr = new double[Nv];
		int counter = 0;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			double hk = this->MaxVoronoiDistance(v);
			max_h_Arr[counter] = hk;
			++counter;
			v = itv.Next();
		}
		double h = *std::max_element(max_h_Arr, max_h_Arr + Nv);
		if (min_of_max)
		{
			(*min_of_max) = *std::min_element(max_h_Arr, max_h_Arr + Nv);
		}
		delete[] max_h_Arr;
		return h;
	}
	return -1;
}
double TriangulationQM::GetPointDistributionMeasure()
{
	if (this->Mesh2D)
	{
		double hmax = this->MaxVoronoiDistance();
		int Nv = this->Mesh2D->NumVertices();
		this->UpdateArea();
		double A = this->Area;
		double hPrime = sqrt(sqrt(12.0) / 9.0 * A / Nv);
		return (hmax / hPrime - 1.0);
	}
	return -1;
}
double TriangulationQM::GetPointDistributionRatio()
{
	if (this->Mesh2D)
	{
		double hmin;
		double hmax = this->MaxVoronoiDistance(&hmin);
		return (hmax / hmin - 1.0);
	}
	return -1;
}
double TriangulationQM::GetRegulatoryMeasure()
{
	if (this->Mesh2D)
	{
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		double chiMax = -1;
		while (v)
		{
			double gamma = this->GetMinDistance(v);//O(n)
			double h = this->MaxVoronoiDistance(v);//O(n)
			double chi = sqrt(3.0) * h / gamma;
			if (chi > chiMax)
				chiMax = chi;
			v = itv.Next();
		}
		return chiMax - 1.0;
	}
	return -1;
}
double TriangulationQM::GetVoronoiArea(Vertex* v)
{
	if (v)
	{
		Vector3D P = v->GetPoint();
		double VoronoiArea = 0;
		FaceIterator itf(v);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				VertexIterator itvv(f);
				Vertex* a = itvv.Next();
				Vertex* b = itvv.Next();
				Vertex* c = itvv.Next();
				if (b == v)
					std::swap(a, b);
				if (c == v)
					std::swap(a, c);
				Vector3D A = a->GetPoint();
				Vector3D B = b->GetPoint();
				Vector3D C = c->GetPoint();
				Vector3D MidAB = (A + B) / 2.0;
				Vector3D A_MidAB = MidAB - A;
				Vector3D MidAC = (A + C) / 2.0;
				Vector3D A_MidAC = MidAC - A;
				double hB = A_MidAB.abs();
				double hC = A_MidAC.abs();
				Triangle t(A, B, C);
				Circle Cir = this->GetCircumCircle(f);
				Vector3D Center = Cir.GetCenter();
				Vector3D A_Center = Center - A;
				double h = A_Center.abs();
				double theta = t.GetMaxAngle();
				double BAC = Vector3D::GetAngle(B, A, C);
				if (theta <= pi / 2.0 + 1.0e-10)
				{
					Triangle t1(A, MidAB, Center);
					Triangle t2(A, MidAC, Center);
					VoronoiArea += t1.GetArea();
					VoronoiArea += t2.GetArea();
				}
				else {

					if (BAC >= pi / 2.0)
					{
						double dA = 0;
						Triangle t1(A, MidAB, Center);
						dA += t1.GetArea();
						Triangle t2(A, MidAC, Center);
						dA += t2.GetArea();
						Vector3D MidBC = (B + C) / 2.0;
						Vector3D Center_MidBC = MidBC - Center;
						double h = Center_MidBC.abs();
						double alpha = Vector3D::GetAngle(MidBC, Center, MidAB);
						double beta = Vector3D::GetAngle(MidBC, Center, MidAC);
						dA -= 0.5 * h * h * tan(alpha);
						dA -= 0.5 * h * h * tan(beta);
						if (dA < 0)
							dA = 0;
						VoronoiArea += dA;
					}
					else
					{
						double smin = std::fmin(hB, hC);
						VoronoiArea += 0.5 * smin * smin * tan(BAC);
					}
				}
			}
			f = itf.Next();
		}
		return VoronoiArea;
	}
	return -1;
}
double TriangulationQM::GetCellVolumeDeviation()
{
	if (this->Mesh2D)
	{
		double Amin = -1;
		double Amax = -1;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			double A = this->GetVoronoiArea(v);//O(n)
			if (Amin < 0 || A < Amin)
				Amin = A;
			if (A > Amax)
				Amax = A;
			v = itv.Next();
		}
		return (Amax / Amin - 1.0);
	}
	return -1;
}
double TriangulationQM::Get2ndMomentTrace(Vertex* v)
{
	//M2(A) = 1/GeoArea(A) * integral_GeoArea(A) (M2_A(x,y) dx dy)
	//M2_A(x, y) = [(x - xA)^2, (x - xA)(y - yA);(x - xA)(y - yA) ,(y - yA)^2 ]
	//Here, instead of itegration the geometric center of triangles A,M_AB,M_ABC and A, M_AC, M_ABC are used.
	if (v)
	{
		Vector3D P = v->GetPoint();
		double SecondMomentTrace = 0;
		FaceIterator itf(v);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				VertexIterator itvv(f);
				Vertex* a = itvv.Next();
				Vertex* b = itvv.Next();
				Vertex* c = itvv.Next();
				if (b == v)
					std::swap(a, b);
				if (c == v)
					std::swap(a, c);
				Vector3D A = a->GetPoint();
				Vector3D B = b->GetPoint();
				Vector3D C = c->GetPoint();
				Vector3D MidAB = (A + B) / 2.0;
				Vector3D MidAC = (A + C) / 2.0;
				Vector3D MidABC = (A + B + C) / 3.0;
				Vector3D P = (MidABC + MidAB + A) / 3.0;
				Vector3D Q = (MidABC + MidAC + A) / 3.0;
				Vector3D AP = P - A;
				Vector3D AQ = Q - A;
				SecondMomentTrace += AP.abs() * AP.abs();
				SecondMomentTrace += AQ.abs() * AQ.abs();
			}
			f = itf.Next();
		}
		return SecondMomentTrace;
	}
	return -1;
}
double TriangulationQM::Get2ndMomentDeviatoricDet(Vertex* v)
{
	//M2_A(x, y) = [(x - xA)^2, (x - xA)(y - yA);(x - xA)(y - yA) ,(y - yA)^2 ]
	//D_A(x,y) = M2_A(x,y) - m * I2x2, m = trace(M2_A(x,y))/2
	//define a = x - xA, b = y = yA, then abs(Det(D_A(x,y))) = (a^2 + b^2) / 4 = dist(x,y,A) ^2 / 4;
	//D(A) = 1/GeoArea(A) * integral_GeoArea(A) (abs(Det(D_A(x,y))) dx dy)
	//Here, instead of itegration the geometric center of triangles A,M_AB,M_ABC and A, M_AC, M_ABC are used.
	return this->Get2ndMomentTrace(v) / 4.0;
}
double TriangulationQM::Get2ndMomentTraceMaxDev()
{
	if (this->Mesh2D)
	{
		int Nv = this->Mesh2D->NumVertices();
		double M2Tavg = 0;//Second moment trace average
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			double M2T = this->Get2ndMomentTrace(v);
			M2Tavg += M2T;
			v = itv.Next();
		}
		M2Tavg /= Nv;
		double MaxDev = 0;
		VertexIterator itvv(this->Mesh2D);
		Vertex* vv = itvv.Next();
		while (vv)
		{
			double M2T = this->Get2ndMomentTrace(vv);
			if (fabs(M2T - M2Tavg) > MaxDev)
				MaxDev = fabs(M2T - M2Tavg);
			vv = itvv.Next();
		}
		return MaxDev;
	}
	return -1;
}
double TriangulationQM::MaxGeoAreaDistance(Vertex* v)
{
	if (v)
	{
		Vector3D P = v->GetPoint();
		double MaxDist = 0;
		FaceIterator itf(v);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				VertexIterator itvv(f);
				Vertex* a = itvv.Next();
				Vertex* b = itvv.Next();
				Vertex* c = itvv.Next();
				if (b == v)
					std::swap(a, b);
				if (c == v)
					std::swap(a, c);
				Vector3D A = a->GetPoint();
				Vector3D B = b->GetPoint();
				Vector3D C = c->GetPoint();
				Vector3D MidAB = (A + B) / 2.0;
				Vector3D MidAC = (A + C) / 2.0;
				Vector3D MidABC = (A + B + C) / 3.0;
				double d1 = A.distance(MidAB);
				double d2 = A.distance(MidAC);
				double d3 = A.distance(MidABC);
				double maxDist = std::fmax(d1, std::fmax(d2, d3));
				if (maxDist > MaxDist)
					MaxDist = maxDist;
			}
			f = itf.Next();
		}
		return MaxDist;
	}
	return -1;
}
double TriangulationQM::MaxGeoAreaDistance(double* min_of_max)
{
	if (this->Mesh2D)
	{
		int Nv = this->Mesh2D->NumVertices();
		double* max_h_Arr = new double[Nv];
		int counter = 0;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			double hk = this->MaxGeoAreaDistance(v);
			max_h_Arr[counter] = hk;
			++counter;
			v = itv.Next();
		}
		double h = *std::max_element(max_h_Arr, max_h_Arr + Nv);
		if (min_of_max)
		{
			(*min_of_max) = *std::min_element(max_h_Arr, max_h_Arr + Nv);
		}
		delete[] max_h_Arr;
		return h;
	}
	return -1;
}
double TriangulationQM::GetPointDistributionMeasure_GeoAreaBased()
{
	if (this->Mesh2D)
	{
		double hmax = this->MaxGeoAreaDistance();
		int Nv = this->Mesh2D->NumVertices();
		this->UpdateArea();
		double A = this->Area;
		double hPrime = sqrt(sqrt(12.0) / 9.0 * A / Nv);
		return (hmax / hPrime - 1.0);
	}
	return -1;
}
double TriangulationQM::GetPointDistributionRatio_GeoAreaBased()
{
	if (this->Mesh2D)
	{
		double hmin;
		double hmax = this->MaxGeoAreaDistance(&hmin);
		return (hmax / hmin - 1.0);
	}
	return -1;
}
double TriangulationQM::GetGeoArea(Vertex* v)
{
	if (v)
	{
		Vector3D P = v->GetPoint();
		double GeoArea = 0;
		FaceIterator itf(v);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				VertexIterator itvv(f);
				Vertex* a = itvv.Next();
				Vertex* b = itvv.Next();
				Vertex* c = itvv.Next();
				if (b == v)
					std::swap(a, b);
				if (c == v)
					std::swap(a, c);
				Vector3D A = a->GetPoint();
				Vector3D B = b->GetPoint();
				Vector3D C = c->GetPoint();
				Vector3D MidAB = (A + B) / 2.0;
				Vector3D MidAC = (A + C) / 2.0;
				Vector3D MidABC = (A + B + C) / 3.0;
				Triangle t1(A, MidAB, MidABC);
				Triangle t2(A, MidAC, MidABC);
				GeoArea += t1.GetArea();
				GeoArea += t2.GetArea();
			}
			f = itf.Next();
		}
		return GeoArea;
	}
	return -1;
}
double TriangulationQM::GetCellVolumeDeviation_GeoAreaBased()
{
	if (this->Mesh2D)
	{
		double Amin = -1;
		double Amax = -1;
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			double A = this->GetGeoArea(v);//O(n)
			if (Amin < 0 || A < Amin)
				Amin = A;
			if (A > Amax)
				Amax = A;
			v = itv.Next();
		}
		return (Amax / Amin - 1.0);
	}
	return -1;
}
double TriangulationQM::GetMaxArea()
{
	if (this->Mesh2D)
	{
		double maxArea = -1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double A = this->GetArea(f);
				if (A > maxArea)
					maxArea = A;
			}
			f = itf.Next();
		}
		return maxArea;
	}
	return -1;
}
double TriangulationQM::GetMinArea()
{
	if (this->Mesh2D)
	{
		double minArea = 1.0e10;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double A = this->GetArea(f);
				if (A < minArea)
					minArea = A;
			}
			f = itf.Next();
		}
		return minArea;
	}
	return -1;
}
double TriangulationQM::GetMaxAreaMeasure()
{
	if (this->Mesh2D)
	{
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		double Amax = 0;
		double Atot = 0;
		int Nt = this->NumTriangles();
		while (f)
		{
			if (f != this->boundary)
			{
				double A = this->GetArea(f);
				if (A > Amax)
					Amax = A;
				Atot += A;
			}
			f = itf.Next();
		}
		double Aavg = Atot / double(Nt);
		return (Amax / Aavg - 1.0);
	}
	return -1;
}
double TriangulationQM::GetCircleRatioMeasure()
{
	if (this->Mesh2D)
	{
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		double qmax = 0;
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetCircumCircleRadius(f);
				double r = this->GetIncircleRadius(f);
				double q = R / r / 2.0;
				double A = this->GetArea(f);
				if (q > qmax)
					qmax = q;
			}
			f = itf.Next();
		}
		return (qmax - 1.0);
	}
	return -1;
}
double TriangulationQM::GetCircumRadDeviationMeasure()
{
	if (this->Mesh2D)
	{
		int Nt = this->NumTriangles();
		double Ravg = 0;
		double R2avg = 0;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetCircumCircleRadius(f);
				Ravg += R;
				R2avg += R * R;
			}
			f = itf.Next();
		}
		Ravg /= Nt;
		R2avg /= Nt;
		double sigma = sqrt(fabs(R2avg - Ravg * Ravg));
		return sigma / Ravg;
	}
	return -1;
}
double TriangulationQM::GetStretchFactor(Vertex* v, shared_ptr<QuadEdgeIndex> ind)//BFS
{
	if (v)
	{
		if (!ind)
			ind = make_shared<QuadEdgeIndex>(this->Mesh2D);
		int iv = ind->GetIndex(v);
		int Nv = this->NumVertices();
		vector<bool> discovered, processed;
		vector<double> distance;
		processed.resize(Nv, false);
		discovered.resize(Nv, false);
		distance.resize(Nv, -1.0);
		queue<int> Q;
		Q.push(iv);
		discovered[iv] = true;
		distance[iv] = 0;
		while (Q.size())
		{
			int iu = Q.front();
			Vertex* u = ind->GetVertex(iu);
			EdgeIterator ite(u);
			Edge* e = ite.Next();
			while (e)
			{
				Vertex* w = (e->GetOrig() == u) ? e->GetDest() : e->GetOrig();
				int iw = ind->GetIndex(w);
				if (!processed[iw])
				{
					double L = e->GetLength();
					if (!discovered[iw])
					{
						Q.push(iw);
						discovered[iw] = true;
						distance[iw] = distance[iu] + L;
					}
					else
					{
						double dist = distance[iu] + L;
						distance[iw] = std::fmin(distance[iw], dist);
					}
				}
				e = ite.Next();
			}
			processed[iu] = true;
			Q.pop();
		}
		double maxSF = 0;
		Vector3D V = v->GetPoint();
		for (unsigned iw = 0; iw < distance.size(); ++iw)
		{
			Vertex* w = ind->GetVertex(iw);
			Vector3D W = w->GetPoint();
			double vw = V.distance(W);
			double SF = distance[iw] / vw;
			if (SF > maxSF)
				maxSF = SF;
		}
		return maxSF;
	}
	return -1;
}
double TriangulationQM::GetStretchFactor()
{
	/*Ref: D2R61 Cui, S., Kanj, I. A., & Xia, G. 2011 On the stretch factor of Delaunay triangulations
	of points in convex position, Computational Geometry: Theory and applications, 44, pp 104-109*/
	if (this->Mesh2D)
	{
		double SF_max = 0;
		shared_ptr<QuadEdgeIndex> ind = make_shared<QuadEdgeIndex>(this->Mesh2D);
		VertexIterator itv(this->Mesh2D);
		Vertex* v = itv.Next();
		while (v)
		{
			double sf = this->GetStretchFactor(v, ind);
			if (sf > SF_max)
				SF_max = sf;
			v = itv.Next();
		}
		return SF_max;
	}
	return -1;
}
double TriangulationQM::GetMinCircumRadius()
{
	if (this->Mesh2D)
	{
		double minRad = -1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetCircumCircleRadius(f);
				if (minRad < 0 || R < minRad)
					minRad = R;
			}
			f = itf.Next();
		}
		return minRad;
	}
	return -1;
}
double TriangulationQM::GetMaxCircumRadius()
{
	if (this->Mesh2D)
	{
		double maxRad = -1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetCircumCircleRadius(f);
				if (R > maxRad)
					maxRad = R;
			}
			f = itf.Next();
		}
		return maxRad;
	}
	return -1;
}
double TriangulationQM::GetMinInRadius()
{
	if (this->Mesh2D)
	{
		double minRad = -1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetIncircleRadius(f);
				if (minRad < 0 || R < minRad)
					minRad = R;
			}
			f = itf.Next();
		}
		return minRad;
	}
	return -1;
}
double TriangulationQM::GetMaxInRadius()
{
	if (this->Mesh2D)
	{
		double maxRad = -1;
		FaceIterator itf(this->Mesh2D);
		Face* f = itf.Next();
		while (f)
		{
			if (f != this->boundary)
			{
				double R = this->GetIncircleRadius(f);
				if (R > maxRad)
					maxRad = R;
			}
			f = itf.Next();
		}
		return maxRad;
	}
	return -1;
}
void TriangulationQM::GetSizeRelations(Face* f, double& Lmin, double& Lmed, double& Lmax)
{
	if (this->Mesh2D && f != this->boundary)
	{
		EdgeIterator ite(f);
		Edge* e1 = ite.Next();
		Edge* e2 = ite.Next();
		Edge* e3 = ite.Next();
		double L1 = e1->GetLength();
		double L2 = e2->GetLength();
		double L3 = e3->GetLength();
		vector<double> L = { L1, L2, L3 };
		std::sort(L.begin(), L.end());
		Lmin = L[0];
		Lmed = L[1];
		Lmax = L[2];
	}
}
double TriangulationQM::GetMinContainmentRadius(Face* f)
{
	if (this->Mesh2D && f != this->boundary)
	{
		Triangle t = this->GetTriangle(f);
		double thetaMax = t.GetMaxAngle();
		if (thetaMax > pi / 2.0 + 1e-10)
		{
			double L1, L2, L3 = -2;
			this->GetSizeRelations(f, L1, L2, L3);
			return (L3 / 2.0);
		}
		else
		{
			return this->GetCircumCircleRadius(f);
		}
	}
	return -1;
}
double TriangulationQM::EdgeOrientationAverage()
{
	if (this->Mesh2D)
	{
		double avg = 0;
		int N = 0;
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetVector();
			double theta = P.GetOrientationAboutXAxis();
			avg += theta;
			++N;
			e = ite.Next();
		}
		avg /= double(N);
		return avg;
	}
	return -1;
}
double TriangulationQM::EdgeOrientationStdDeviation()
{
	if (this->Mesh2D)
	{
		double var = 0;
		double avg = this->EdgeOrientationAverage();
		int N = 0;
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetVector();
			double theta = P.GetOrientationAboutXAxis();
			var += (theta - avg) * (theta - avg);
			++N;
			e = ite.Next();
		}
		double N1 = N - 1;
		var /= N1;
		return sqrt(var);
	}
	return -1;
}
double TriangulationQM::EdgeOrientationSkewness()
{
	if (this->Mesh2D)
	{
		double avg = this->EdgeOrientationAverage();
		double sum = 0;
		int N = 0;
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetVector();
			double theta = P.GetOrientationAboutXAxis();
			sum += (theta - avg) * (theta - avg) * (theta - avg);
			++N;
			e = ite.Next();
		}
		double skewness = sum / double(N);
		return skewness;
	}
	return -1;
}
double TriangulationQM::EdgeOrientationKurtosis()
{
	if (this->Mesh2D)
	{
		double avg = this->EdgeOrientationAverage();
		double sum = 0;
		int N = 0;
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D P = e->GetVector();
			double theta = P.GetOrientationAboutXAxis();
			sum += (theta - avg) * (theta - avg) * (theta - avg) * (theta - avg);
			++N;
			e = ite.Next();
		}
		double Kurt = sum / double(N);
		return Kurt;
	}
	return -1;
}
double TriangulationQM::EdgeOrientationCoverage()
{
	if (this->Mesh2D)
	{
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		int N = this->Mesh2D->NumEdges();
		double dtheta = pi / double(N);
		vector<int> boxes;
		boxes.resize(N, 0);
		while (e)
		{
			Vector3D P = e->GetVector();
			double theta = P.GetOrientationAboutXAxis();
			int boxNum = int(fabs(theta) / dtheta);
			if (boxNum >= N)
				boxNum = N - 1;
			++boxes[boxNum];
			e = ite.Next();
		}
		double chiSquared = 0;
		for (int i = 0; i < boxes.size(); ++i)
		{
			chiSquared += (boxes[i] - 1) * (boxes[i] - 1);
		}
		return chiSquared;
	}
	return -1;
}
double TriangulationQM::EdgeOrientationCoverage_10boxes()
{
	//difference: box sizes are fixed: pi/10
	if (this->Mesh2D)
	{
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		int N = 10;//difference
		double dtheta = pi / double(N);
		vector<double> boxes;
		boxes.resize(N, 0);
		while (e)
		{
			Vector3D P = e->GetVector();
			double theta = P.GetOrientationAboutXAxis();
			int boxNum = int(fabs(theta) / dtheta);
			if (boxNum >= N)
				boxNum = N - 1;
			++boxes[boxNum];
			e = ite.Next();
		}
		double chiSquared = 0;
		for (int i = 0; i < boxes.size(); ++i)
		{
			chiSquared += (boxes[i] - 1) * (boxes[i] - 1);
		}
		return chiSquared;
	}
	return -1;
}
double TriangulationQM::GetAvgVelDotEdge(function<double(const Vector3D& P)> V[2], double& minVal, double& maxVal, double& piVal)
{
	if (this->Mesh2D)
	{
		maxVal = -1;
		minVal = 1;
		piVal = 1;
		double sum = 0;
		int counter = 0;
		EdgeIterator ite(this->Mesh2D);
		Edge* e = ite.Next();
		while (e)
		{
			Vector3D E = e->GetVector();
			Vector3D point = e->GetPoint();
			double u = V[0](point);
			double v = V[1](point);
			Vector3D vel(u, v);
			double absCosTheta = fabs(vel.abs() * E.abs() > 1e-10 ? (vel || E) / vel.abs() / E.abs() : 0);
			sum += absCosTheta;
			if (maxVal < absCosTheta)
				maxVal = absCosTheta;
			if (absCosTheta < minVal)
				minVal = absCosTheta;
			piVal *= absCosTheta;
			++counter;
			e = ite.Next();
		}
		return (sum / double(counter));
	}
	return -1;
}
bool TriangulationQM::ShiftPointsRandomly(Triangulation* T)
{
	double xMin = 9999999, yMin = 999999, xMax = -9999999, yMax = -999999;
	VertexIterator itv(T->GetMesh2D());
	Vertex* v = itv.Next();
	while (v)
	{
		if (v->GetX() < xMin)
			xMin = v->GetX();
		if (v->GetX() > xMax)
			xMax = v->GetX();
		if (v->GetY() < yMin)
			yMin = v->GetY();
		if (v->GetY() > yMax)
			yMax = v->GetY();
		v = itv.Next();
	}
	double hx = (xMax - xMin) * 0.02;
	double hy = (yMax - yMin) * 0.02;
	xMax -= hx;
	xMin += hx;
	yMax -= hy;
	yMin += hy;
	double h = sqrt(hx * hx + hy * hy);
	double resolution = min(hx, hy) / 100.0;
	EdgeIterator ite(T->GetMesh2D());
	Edge* e = ite.Next();
	double alpha = 0;//a number between -0.5 to 0.5 
	int m = 100;
	while (e)
	{
		if (e->GetLeft() != T->GetBoundary() && e->GetRight() != T->GetBoundary())
		{
			Vertex* vD = e->GetDest();
			bool use_vD = !T->IsOnBoundary(vD);
			if (use_vD)
			{
				alpha = double(abs(rand() % m));
				alpha = alpha / double(m) - 0.5;
				Vector3D O = e->GetOrig()->GetPoint();
				Vector3D D = vD->GetPoint();
				Vector3D P = D * (1.0 - alpha) + O * alpha;
				if (P(0) > xMin && P(0) < xMax && P(1) > yMin && P(1) < yMax)
				{
					vD->SetX(P(0));
					vD->SetY(P(1));
				}
			}
		}
		e = ite.Next();
	}
	bool go_on = false;
	int counter = 0;
	int maxCounter = 100;
	do
	{
		++counter;
		VertexIterator itvv(T->GetMesh2D());
		Vertex* vv = itvv.Next();
		while (vv)
		{
			if (!T->IsOnBoundary(vv))
			{
				Vector3D P = vv->GetPoint();
				VertexIterator itv2(T->GetMesh2D());
				Vertex* v2 = itv2.Next();
				while (v2)
				{
					if (vv != v2)
					{
						Vector3D P2 = v2->GetPoint();
						if (P.distance(P2) < resolution)
						{
							go_on = true;
							double backOff = h;// 2 * P.distance(P2) / resolution;
							Vector3D Q = P * (1 + backOff) + P2 * (-backOff);
							if (Q(0) > xMin && Q(0) < xMax && Q(1) > yMin && Q(1) < yMax)
							{
								vv->SetX(Q(0));
								vv->SetY(Q(1));
							}
							else {
								Q = P2 * (1 + backOff) + P * (-backOff);
								if (Q(0) > xMin && Q(0) < xMax && Q(1) > yMin && Q(1) < yMax && !T->IsOnBoundary(v2))
								{
									v2->SetX(Q(0));
									v2->SetY(Q(1));
								}
								else {
									return false;
								}
							}
							break;
						}
					}
					v2 = itv2.Next();
				}
			}
			vv = go_on ? 0 : itvv.Next();
		}
	} while (go_on && counter < maxCounter);
	if (counter == maxCounter)
		return false;
	vector<Vector3D> points;
	T->GetPoints(points);
	(*T) = DelaunayLifting::Triangulate(points);
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
	T0.Write("..\\Data\\test8.bqe");
	Triangulation T;
	T.Read("..\\Data\\test8.bqe");
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

