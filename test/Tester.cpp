#include <iostream>
#include "Tester.h"
#include "test_MathUtils.h"
#include "test_bMatrix.h"
#include "../src/Vector.h"
#include "../src/Vector3D.h"
#include "../src/SquareMatrix.h"
#include "../src/CoordTransfer.h"
#include "../src/Line.h" 
#include "../src/Plane.h"
#include "../src/Circle.h"
#include "../src/LineSegment2D.h"
#include "../src/Triangle.h"
#include "../src/Vertex.h"
#include "../src/Face.h"
#include "../src/Edge.h"
#include "../src/EdgeList.h"
#include "../src/EdgeIterator.h"
#include "../src/FaceIterator.h"
#include "../src/VertexIterator.h"
#include "../src/QuadEdge.h"
#include "../src/QuadEdgeIndex.h"
#include "../src/Triangulation.h"
#include "../src/TriangulationQM.h"
#include "../src/ConvexHull3D.h"
#include "../src/DelaunayLifting.h"
#include "../src/bPolygon.h"
#include "../src/WeatherillHassan.h"
#include "../src/RuppertShewchuk.h"
#include "Test_AdvancingFront.h"
#include "../src/ShapeFunction.h"
#include "../src/bGrid.h"
#include "../src/Node.h"
#include "../src/NodeComposite.h"
#include "../src/FVM_Grid.h"
#include "../src/ConductionConvectionProblem.h"
#include "../src/FVM.h"
#include "../src/BitMap.h"
void Tester()
{
	int NumTests = 0;
	try
	{
		std::cout << std::endl << "MathUtils test ...";
		bool correct = tester_MathUtils(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "bMatrix test ...";
		correct = tester_bMatrix(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Vector test ...";
		correct = tester_Vector(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Vector3D test ...";
		correct = tester_Vector3D(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "SquareMatrix test ...";
		correct = tester_SquareMatrix(NumTests);//+6
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "CoordTransfer test ...";
		correct = tester_CoordTransfer(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Line test ...";
		correct = tester_Line(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Plane test ...";
		correct = tester_Plane(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Circle test ...";
		correct = tester_Circle(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "LineSegment2D test ...";
		correct = tester_LineSegment2D(NumTests);//+4
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Triangle test ...";
		correct = tester_Triangle(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Vertex test ...";
		correct = tester_Vertex(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Face test ...";
		correct = tester_Face(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Edge test ...";
		correct = tester_Edge(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "EdgeList test ...";
		correct = tester_EdgeList(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "EdgeIterator test ...";
		correct = tester_EdgeIterator(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "FaceIterator test ...";
		correct = tester_FaceIterator(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "VertexIterator test ...";
		correct = tester_VertexIterator(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "QuadEdge test ...";
		correct = tester_QuadEdge(NumTests);//+7
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "QuadEdgeIndex test ...";
		correct = tester_QuadEdgeIndex(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Triangulation test ...";
		correct = tester_Triangulation(NumTests);//+16
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "TriangulationQM test ...";
		correct = tester_TriangulationQM(NumTests);//+3
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "ConvexHull3D test ...";
		correct = tester_ConvexHull3D(NumTests);//+7
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "DelaunayLifting test ...";
		correct = tester_DelaunayLifting(NumTests);//+11
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Polygon test ...";
		correct = tester_bPolygon(NumTests);//+9
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "RuppertShewchuk test ...";
		correct = tester_RuppertShewchuk(NumTests);//+11
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Tessellation_WeatherillHassan test ...";
		correct = test_WeatherillHassan(NumTests);//+5
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "AdvancingFront test ...";
		correct = tester_AdvancingFront(NumTests);//+10
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "ShapeFunction test ...";
		correct = tester_ShapeFunction(NumTests);//+3
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Grid test ...";
		correct = tester_bGrid(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Node test ...";
		correct = tester_Node(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "NodeComposite test ...";
		correct = tester_NodeComposite(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "ConductionConvectionProblem test ...";
		correct = tester_ConductionConvectionProblem(NumTests);//+2
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "FVM test ...";
		correct = tester_FVM(NumTests);//+18
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "tester_FVM_Grid test ...";
		correct = tester_FVM_Grid(NumTests);//+10
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return;
		//---------------------------------------------------------------------------------
		int Expected_Num_Tests = 141;
		if (NumTests < Expected_Num_Tests)
		{
			std::cout << std::endl  << Expected_Num_Tests - NumTests << " TESTS ARE BYPASSED!!!!";
		}
		else
		{
			std::cout << std::endl << "All tests are okay!";
		}
	}
	catch (const char* e)
	{
		std::cout << std::endl << "Code Error: " << e;
	}
	catch (std::exception& e)
	{
		std::cout << std::endl << "Standard Exception: " << e.what();
	}
}