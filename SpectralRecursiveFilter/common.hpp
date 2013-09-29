#pragma once

//please change this options if you want.
//#define USE_SSE			//providing SSE computation
#define USE_OPENCV2		//providing OpenCV2 interface

const double PI=3.1415926535897932384626433832795;

//extrapolation functions
//(atE() and atS() require variables w and h in their scope respectively)
//AA|ABCDE|EE (?)
//#define atN(y) (0)
//#define atW(x) (0)
//#define atS(y) (h-1)
//#define atE(x) (w-1)
//CB|ABCDE|DC (Reflection-101)
#define atN(y) (abs(y))
#define atW(x) (abs(x))
#define atS(y) (h-1-abs(h-1-(y)))
#define atE(x) (w-1-abs(w-1-(x)))
