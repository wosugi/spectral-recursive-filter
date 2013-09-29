#pragma once
#include <stdexcept>
#include <vector>
#include <cassert>
#include <cmath>
#include "common.hpp"

#ifdef USE_OPENCV2
#include <opencv2/opencv.hpp>
#endif

namespace convolution
{
//*************************************************************************************************

class gauss
{
private:
	double sx,sy; //scale of Gaussian
	std::vector<double> kx,ky; //kernel
	
public:
	gauss(double sx,double sy,double trunc=5.0):sx(sx),sy(sy)
	{
		if(sx<0.0 || sy<0.0)
			throw std::invalid_argument("\'sx\' and \'sy\' should be nonnegative!");
	
		kx=gen_kernel(sx,int(trunc*sx+0.9999999));
		ky=gen_kernel(sy,int(trunc*sy+0.9999999));
	}
	~gauss(){}

private:
	static inline std::vector<double> gen_kernel(double s,int r)
	{
		if(s==0.0)
			return std::vector<double>();

		const double gamma=1.0/(2.0*s*s);
		std::vector<double> kernel(1+r);
		for(int t=0;t<=r;++t)
			kernel[t]=exp(-gamma*t*t);

		double sum=0.0;
		for(int t=0;t<=r;++t)
			sum+=kernel[t];
		sum=2.0*sum-kernel[0];

		for(int t=0;t<=r;++t)
			kernel[t]/=sum;
		return kernel;
	}
	
	template <typename T>
	inline void filter_h(int w,int h,T* src,T* dst)
	{
		throw std::invalid_argument("Unsupported element type!");
	}
	template <typename T>
	inline void filter_v(int w,int h,T* src,T* dst)
	{
		throw std::invalid_argument("Unsupported element type!");
	}

public:
	template <typename T>
	void filter(int w,int h,T* src,T* dst)
	{
		if(w<=4.0*sx || h<=4.0*sy)
			throw std::invalid_argument("\'sx\' and \'sy\' should be less than about w/4 or h/4!");

		if(kx.empty() && ky.empty())
			return;
		else if(kx.empty())
			filter_v(w,h,src,dst);
		else if(ky.empty())
			filter_h(w,h,src,dst);
		else
		{
			filter_v(w,h,src,dst);
			filter_h(w,h,dst,dst); //filter_h() allows src==dst.
		}
	}
	
	#ifdef USE_OPENCV2 //OpenCV2 interface for easy function call.
	void filter(const cv::Mat& src,cv::Mat& dst)
	{
		//checking the format of input/output images
		if(src.size()!=dst.size())
			throw std::invalid_argument("\'src\' and \'dst\' should have the same size!");
		if(src.type()!=dst.type())
			throw std::invalid_argument("\'src\' and \'dst\' should have the same element type!");
		if(src.channels()!=1 || dst.channels()!=1)
			throw std::invalid_argument("Multi-channel images are unsupported!");
		if(src.isSubmatrix() || dst.isSubmatrix())
			throw std::invalid_argument("Subimages are unsupported!");

		switch(src.type())
		{
		case CV_32FC1:
			filter<float>(src.cols,src.rows,reinterpret_cast<float*>(src.data),reinterpret_cast<float*>(dst.data));
			break;
		default:
			throw std::invalid_argument("The element type is unsupported yet!");
			break;
		}
	}
	#endif
};

//*************************************************************************************************

template <>
inline void gauss::filter_h<float>(int w,int h,float* src,float* dst)
{
	const int r=int(kx.size())-1;
	std::vector<float> kernel(r+1);
	for(int t=0;t<=r;++t)
		kernel[t]=float(kx[t]);

	std::vector<float> buf(w); //for src==dst
	for(int y=0;y<h;++y)
	{
		float* p=&src[y*w];
		float* q=&dst[y*w];
		std::copy(p,p+w,buf.begin());
		for(int x=0;x<w;++x)
		{
			float sum=kernel[0]*p[x];
			for(int t=1;t<=r;++t)
				sum+=kernel[t]*(buf[atW(x-t)]+buf[atE(x+t)]);
			q[x]=sum;
		}
	}
}
template <>
inline void gauss::filter_v<float>(int w,int h,float* src,float* dst)
{
	const int r=int(ky.size())-1;
	std::vector<float> kernel(r+1);
	for(int t=0;t<=r;++t)
		kernel[t]=float(ky[t]);

	for(int y=0;y<h;++y)
	{
		float* q=&dst[y*w];
		for(int x=0;x<w;++x)
		{
			float sum=kernel[0]*src[x+y*w];
			for(int t=1;t<=r;++t)
				sum+=kernel[t]*(src[x+atN(y-t)*w]+src[x+atS(y+t)*w]);
			q[x]=sum;
		}
	}
}

//*************************************************************************************************
}
