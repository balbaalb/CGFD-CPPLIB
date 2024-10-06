#ifndef AdvancingFrontH
#define AdvancingFrontH
#include <vector>
#include <memory>
#include <fstream>
#include "Vector3D.h"
#include "Triangulation.h"
using namespace std;
class bPolygon;
class FrontNode
{
	friend class AdvancingFront;
	Vector3D* point;
	FrontNode* next;
	FrontNode* prev;
	FrontNode(const Vector3D& P);
	FrontNode(FrontNode* rhs);
	~FrontNode();
};
class AdvancingFront
{
	int size;
	FrontNode* head;
	double LengthScale;
	vector<Vector3D> BoundaryPoints;
	vector<Vector3D> InsertedNodes;
	std::unique_ptr<AdvancingFront> child1, child2;
	bool Completed;
	void AddToTail(const Vector3D& P);
	void AddAfterHead(const Vector3D& P);
	Vector3D SuggestPoint() const;
	bool IsNewPointAcceptable(const Vector3D& C) const;
	void AdmitNewPointOverHead(const Vector3D& C);
	bool IsInsideFront(const Vector3D& C) const;//O(n)
	bool IsOnOrInsideFront(const Vector3D& C) const;//O(n)
	bool IsIntersectingFront(const Vector3D& C, bool existingNode = false) const;//O(n)
	bool IsTooCloseToFront(const Vector3D& C, double LengthScale_) const;//O(n)
	FrontNode* NodeClosestTo(const Vector3D& C);//O(n^2) //Needs reduction O(n) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	void RemoveNode(FrontNode* n);
	void CombineInsertedNodes(const AdvancingFront& Front);//O(n)
	void DivideFront(FrontNode* m);//creates child1 and child2
	void AdvanceChildren1Step();
	void MakeComplete();
	bool FrontIntegrity() const;
	void UpdateSize();
	double GetMaxFrontLength() const;
	void report();
	void GetPoints(vector<Vector3D>& points);
	static int counter;
	static ofstream file;
public:
	AdvancingFront(const vector<Vector3D>& input);//should be CCW and contain no holes
	AdvancingFront(FrontNode* h, double LengthScale_);//should be CCW and contain no holes
	~AdvancingFront();
	void Advance();//O(n^3)
	shared_ptr<bPolygon> CreatePolygonFromFront() const;//O(n)
	shared_ptr<bPolygon> CreateLastChildPolygonFromFront() const;
	AdvancingFront* LastChild();
	bool IsCompleted() const;
	Vector3D GetHeadPoint() const;
	Vector3D GetHeadNextPoint() const;
	static bool printSteps;
	static bool AdvanceOneStep;
	static Triangulation Tessellate(const vector<Vector3D>& BoundaryPoints);
};
bool tester_AdvancingFront(int& NumTests);
bool tester_AdvancingFront_1(int& NumTests);
bool tester_AdvancingFront_2(int& NumTests);
bool tester_AdvancingFront_3(int& NumTests);
bool tester_AdvancingFront_4(int& NumTests);
bool tester_AdvancingFront_5(int& NumTests);
bool tester_AdvancingFront_6(int& NumTests);
bool tester_AdvancingFront_7(int& NumTests);
bool tester_AdvancingFront_8(int& NumTests);
bool tester_AdvancingFront_9(int& NumTests);
bool tester_AdvancingFront_10(int& NumTests);//Triangle
#endif