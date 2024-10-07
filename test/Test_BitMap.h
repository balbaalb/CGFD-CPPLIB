#ifndef Test_BitMapH
#define Test_BitMapH
#include "../src/BitMap.h"
bool tester_BitMap(int& numTests)
{
	{
		bBitMap::BitMap bmp(800, 600);
		for (int i = 0; i < 800; ++i)
		{
			int x = i;
			int y = x / 2;
			bmp.Draw1Pixel(x, y, 255, 0, 0);
		}
		bmp.GenerateBMP("..\\Data\\testBmp1.bmp");
	}
	++numTests;
	return true;
}
#endif