#ifndef BITMAP_H
#define BITMAP_H
#include <string>
#include <cstdint>
#include <memory>
#pragma pack(2)
namespace bBitMap
{
	struct BitMapFileHeader
	{
		char HeaderBM[2]{ 'B', 'M' };
		std::int32_t fSize;
		std::int32_t Spacer{ 0 };
		std::int32_t fOffBits;
	}; 
	struct BitMapInfoHeader
	{
		std::int32_t Size{ 40 };
		std::int32_t Width;
		std::int32_t Height;
		std::int16_t Planes{ 1 };
		std::int16_t BitsPerPixel{ 24 };
		std::int32_t Compression{ 0 };
		std::int32_t SizeImage{ 0 };
		std::int32_t XRes{ 2400 };
		std::int32_t YRes{ 2400 };
		std::int32_t ClrUsed{ 0 };
		std::int32_t ClrImportant{ 0 };
	};
	class BitMap
	{
		int Width{ 0 };
		int Height{ 0 };
		int PixelsBites{ 0 };
		std::unique_ptr<uint8_t[]> Pixels{ nullptr };
	public:
		BitMap(int width, int height);
		virtual ~BitMap();
		void Draw1Pixel(int x, int y, uint8_t red, uint8_t green, uint8_t blue);
		bool GenerateBMP(std::string bmpFileName);
	};
}
#endif

