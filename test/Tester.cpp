#include <iostream>
#include "Tester.h"
#include "Test_MathUtils.h"
#include "Test_bMatrix.h"
#include "Test_Vector.h"
#include "Test_Vector3D.h"
#include "Test_SquareMatrix.h"
#include "Test_CoordTransfer.h"
#include "Test_Line.h" 
#include "Test_Plane.h"
#include "Test_Circle.h"
#include "Test_LineSegment2D.h"
#include "Test_Triangle.h"
#include "Test_Vertex.h"
#include "Test_Face.h"
#include "Test_Edge.h"
#include "Test_EdgeList.h"
#include "Test_EdgeIterator.h"
#include "Test_FaceIterator.h"
#include "Test_VertexIterator.h"
#include "Test_QuadEdge.h"
#include "Test_QuadEdgeIndex.h"
#include "Test_Triangulation.h"
#include "Test_TriangulationQM.h"
#include "Test_ConvexHull3D.h"
#include "Test_DelaunayLifting.h"
#include "Test_bPolygon.h"
#include "Test_WeatherillHassan.h"
#include "Test_RuppertShewchuk.h"
#include "Test_AdvancingFront.h"
#include "Test_ShapeFunction.h"
#include "Test_bGrid.h"
#include "Test_Node.h"
#include "Test_NodeComposite.h"
#include "Test_FVM_Grid.h"
#include "Test_ConductionConvectionProblem.h"
#include "Test_FVM.h"
#include "Test_BitMap.h"
int Tester()
{
	int NumTests = 0;
	try
	{
		std::cout << std::endl << "MathUtils test ...";
		bool correct = tester_MathUtils(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "bMatrix test ...";
		correct = tester_bMatrix(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Vector test ...";
		correct = tester_Vector(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Vector3D test ...";
		correct = tester_Vector3D(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "SquareMatrix test ...";
		correct = tester_SquareMatrix(NumTests);//+6
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "CoordTransfer test ...";
		correct = tester_CoordTransfer(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Line test ...";
		correct = tester_Line(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Plane test ...";
		correct = tester_Plane(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Circle test ...";
		correct = tester_Circle(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "LineSegment2D test ...";
		correct = tester_LineSegment2D(NumTests);//+4
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Triangle test ...";
		correct = tester_Triangle(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct) 
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Vertex test ...";
		correct = tester_Vertex(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Face test ...";
		correct = tester_Face(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Edge test ...";
		correct = tester_Edge(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "EdgeList test ...";
		correct = tester_EdgeList(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "EdgeIterator test ...";
		correct = tester_EdgeIterator(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "FaceIterator test ...";
		correct = tester_FaceIterator(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "VertexIterator test ...";
		correct = tester_VertexIterator(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "QuadEdge test ...";
		correct = tester_QuadEdge(NumTests);//+7
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "QuadEdgeIndex test ...";
		correct = tester_QuadEdgeIndex(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Triangulation test ...";
		correct = tester_Triangulation(NumTests);//+16
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "TriangulationQM test ...";
		correct = tester_TriangulationQM(NumTests);//+3
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "ConvexHull3D test ...";
		correct = tester_ConvexHull3D(NumTests);//+7
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "DelaunayLifting test ...";
		correct = tester_DelaunayLifting(NumTests);//+11
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Polygon test ...";
		correct = tester_bPolygon(NumTests);//+9
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "RuppertShewchuk test ...";
		correct = tester_RuppertShewchuk(NumTests);//+11
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Tessellation_WeatherillHassan test ...";
		correct = test_WeatherillHassan(NumTests);//+5
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "AdvancingFront test ...";
		correct = tester_AdvancingFront(NumTests);//+10
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "ShapeFunction test ...";
		correct = tester_ShapeFunction(NumTests);//+3
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Grid test ...";
		correct = tester_bGrid(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "Node test ...";
		correct = tester_Node(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "NodeComposite test ...";
		correct = tester_NodeComposite(NumTests);//+1
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "ConductionConvectionProblem test ...";
		correct = tester_ConductionConvectionProblem(NumTests);//+2
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "FVM test ...";
		correct = tester_FVM(NumTests);//+18
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
		//---------------------------------------------------------------------------------
		std::cout << std::endl << "tester_FVM_Grid test ...";
		correct = tester_FVM_Grid(NumTests);//+10
		std::cout << (correct ? "okay!" : "ERROR!");
		if (!correct)
			return 1;
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
	return 0;
}