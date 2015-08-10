#pragma once

//please change this options if you want.
#define USE_SSE			//providing SSE computation
#define USE_OPENCV2		//providing OpenCV2 interface

const double PI=3.1415926535897932384626433832795;

//#define USE_BORDER_REPLICATE 1	 //BORDER_REPLICATE   = 1, //!< `aaaaaa|abcdefgh|hhhhhhh`
//#define USE_BORDER_REFLECT 1	 //BORDER_REFLECT     = 2, //!< `fedcba|abcdefgh|hgfedcb`
#define USE_BORDER_REFLECT_101 1//BORDER_REFLECT_101 = 4, //!< `gfedcb|abcdefgh|gfedcba`

//extrapolation functions
//(atE() and atS() require variables w and h in their scope respectively)
#ifdef USE_BORDER_REPLICATE 
#define atW(x) (std::max(x,0))
#define atN(y) (std::max(y,0))
#define atE(x) (std::min(x,w-1))
#define atS(y) (std::min(y,h-1))
#elif USE_BORDER_REFLECT 
#define atW(x) ( x < 0 ? std::abs(x+1) : std::abs(x))
#define atN(y) ( y < 0 ? std::abs(y+1) : std::abs(y))
#define atE(x) ( x < w ? w-1-std::abs(w-1-(x)) : w-std::abs(w-1-(x)))
#define atS(y) ( y < h ? h-1-std::abs(h-1-(y)) : h-std::abs(h-1-(y)))
#elif USE_BORDER_REFLECT_101
#define atW(x) (std::abs(x))
#define atN(y) (std::abs(y))
#define atE(x) (w-1-std::abs(w-1-(x)))
#define atS(y) (h-1-std::abs(h-1-(y)))
#endif