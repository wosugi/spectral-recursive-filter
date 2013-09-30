#pragma once

//please change this options if you want.
#define USE_SSE			//providing SSE computation
#define USE_OPENCV2		//providing OpenCV2 interface

const double PI=3.1415926535897932384626433832795;

//extrapolation functions
//(atE() and atS() require variables w and h in their scope respectively)
/*
/// AA|ABCDE|EE (cv::BORDER_REPLICATE)
#define atW(x) (std::max(x,0))
#define atN(y) (std::max(y,0))
#define atE(x) (std::min(x,w-1))
#define atS(y) (std::min(y,h-1))
/*/
/// CB|ABCDE|DC (cv::BORDER_REFLECT_101)
#define atW(x) (std::abs(x))
#define atN(y) (std::abs(y))
#define atE(x) (w-1-std::abs(w-1-(x)))
#define atS(y) (h-1-std::abs(h-1-(y)))
//*/
