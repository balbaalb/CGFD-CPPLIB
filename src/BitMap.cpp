#include <iostream>
#include <fstream>
//#include "BitMapFileHeader.h"
//#include "BitMapInfoHeader.h"
#include "BitMap.h"
bBitMap::BitMap::BitMap(int width, int height) : Width(width), Height(height), PixelsBites(width* height * 3), Pixels(new uint8_t[PixelsBites]{})
{
	//Create a White Canvas
	for (int i = 0; i < Width; ++i)
	{
		for (int j = 0; j < Height; ++j)
		{
			this->Draw1Pixel(i, j, 255, 255, 255);
		}
	}
}
bBitMap::BitMap::~BitMap()
{
}
void bBitMap::BitMap::Draw1Pixel(int x, int y, uint8_t red, uint8_t green, uint8_t blue)
{
	uint8_t* pPixel = Pixels.get() + (y * Width + x) * 3;
	pPixel[0] = blue;
	pPixel[1] = green;
	pPixel[2] = red;
}
bool bBitMap::BitMap::GenerateBMP(std::string bmpFileName)
{
	BitMapFileHeader fileHeader;
	fileHeader.fOffBits = sizeof(BitMapFileHeader) + sizeof(BitMapInfoHeader);
	fileHeader.fSize = fileHeader.fOffBits + PixelsBites;
	BitMapInfoHeader infoHeader;
	infoHeader.Width = Width;
	infoHeader.Height = Height;
	std::ofstream bmpFile;
	bmpFile.open(bmpFileName, std::ios::binary);
	if (!bmpFile)
		return false;
	bmpFile.write((char*)&fileHeader, sizeof(fileHeader));
	bmpFile.write((char*)&infoHeader, sizeof(infoHeader));
	bmpFile.write((char*)Pixels.get(), PixelsBites);
	bmpFile.close();
	if (!bmpFile)
		return false;
	//alternative approach:
	//https://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
	return true;
}
