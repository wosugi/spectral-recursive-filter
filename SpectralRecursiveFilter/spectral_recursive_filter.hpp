#pragma once
#include <stdexcept>
#include <vector>
#include <cassert>
#include <cmath>
#include "common.hpp"

#ifdef USE_OPENCV2
#include <opencv2/opencv.hpp>
#endif

namespace spectral_recursive_filter
{
//*************************************************************************************************

class gauss
{
private:
	static const int K=2; //order of approximation (do not change!)

private:
	double sx,sy; //scale of Gaussian
	int rx,ry; //kernel radius
	std::vector<double> spectX,spectY;
	std::vector<double> tableX,tableY; //look-up tables

public:
	gauss(double sx,double sy):sx(sx),sy(sy)
	{
		if(sx<0.0 || sy<0.0)
			throw std::invalid_argument("\'sx\' and \'sy\' should be nonnegative!");
	
		rx=estimate_radius(sx);
		ry=estimate_radius(sy);
		spectX=gen_spectrum(sx,rx);
		spectY=gen_spectrum(sy,ry);
		tableX=build_lookup_table(rx,spectX);
		tableY=build_lookup_table(ry,spectY);
	}
	~gauss(){}

private:
	static inline double phase(int r)
	{
		return 2.0*PI/(r+1+r); //DCT/DST-5
	}
	static inline int estimate_radius(double s)
	{
		//return (s<4.0) ? int(3.3333*s-0.3333+0.5) : int(3.4113*s-0.6452+0.5); //K==3
		return (s<4.0) ? int(3.0000*s-0.2000+0.5) : int(3.0000*s+0.5); //K==2
	}
	static inline std::vector<double> gen_spectrum(double s,int r)
	{
		const double phi=phase(r);
		std::vector<double> spect(K);
		for(int k=1;k<=K;k++)
			spect[k-1]=2.0*exp(-0.5*s*s*phi*phi*k*k);
		return spect;
	}
	static inline std::vector<double> build_lookup_table(int r,std::vector<double>& spect)
	{
		assert(spect.size()==K);
		const double phi=phase(r);
		std::vector<double> table(K*(1+r));
		for(int u=0;u<=r;++u)
			for(int k=1;k<=K;++k)
				table[K*u+k-1]=cos(k*phi*u)*spect[k-1];
		return table;
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
	template <typename T>
	inline void filter_sse_h(int w,int h,T* src,T* dst)
	{
		throw std::invalid_argument("Unsupported element type!");
	}
	template <typename T>
	inline void filter_sse_v(int w,int h,T* src,T* dst)
	{
		throw std::invalid_argument("Unsupported element type!");
	}

public:

#ifdef USE_SSE
	template <typename T>
	void filter(int w,int h,T* src,T* dst)
	{
		if(w<=4.0*sx || h<=4.0*sy)
			throw std::invalid_argument("\'sx\' and \'sy\' should be less than about w/4 or h/4!");
		
		//filtering is skipped if s==0.0
		if(sx==0.0 && sy==0.0)
			return;
		else if(sx==0.0)
			filter_sse_v<T>(w,h,src,dst);
		else if(sy==0.0)
			filter_sse_h<T>(w,h,src,dst);
		else
		{
			filter_sse_v<T>(w,h,src,dst);
			filter_sse_h<T>(w,h,dst,dst); //only filter_h() allows src==dst.
		}
	}
#else
	template <typename T>
	void filter(int w,int h,T* src,T* dst)
	{
		if(w<=4.0*sx || h<=4.0*sy)
			throw std::invalid_argument("\'sx\' and \'sy\' should be less than about w/4 or h/4!");
		
		//filtering is skipped if s==0.0
		if(sx==0.0 && sy==0.0)
			return;
		else if(sx==0.0)
			filter_v<T>(w,h,src,dst);
		else if(sy==0.0)
			filter_h<T>(w,h,src,dst);
		else
		{
			filter_v<T>(w,h,src,dst);
			filter_h<T>(w,h,dst,dst); //only filter_h() allows src==dst.
		}
	}
#endif
	
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
			throw std::invalid_argument("Unsupported element type!");
			break;
		}
	}
#endif
};

//*************************************************************************************************

template<>
inline void gauss::filter_h<float>(int w,int h,float* src,float* dst)
{
	const int r=rx;
	const float norm=float(1.0/(r+1+r));
	std::vector<float> table(tableX.size());
	for(int t=0;t<int(table.size());++t)
		table[t]=float(tableX[t]);
	
	const float cf11=float(table[K*1+0]*2.0/spectX[0]), cfR1=table[K*r+0];
	const float cf12=float(table[K*1+1]*2.0/spectX[1]), cfR2=table[K*r+1];

	float sum,a1,a2,b1,b2;
	float d,e,delta;
	std::vector<float> buf(w); //to allow for src==dst
	for(int y=0;y<h;++y)
	{
		std::copy(&src[w*y],&src[w*y+w],buf.begin());

		sum=buf[0];
		a1=buf[0]*table[0]; b1=buf[1]*table[0];
		a2=buf[0]*table[1]; b2=buf[1]*table[1];
		for(int u=1;u<=r;++u)
		{
			const float sumA=buf[atW(0-u)]+buf[0+u];
			const float sumB=buf[atW(1-u)]+buf[1+u];
			sum+=sumA;
			a1+=sumA*table[K*u+0]; b1+=sumB*table[K*u+0];
			a2+=sumA*table[K*u+1]; b2+=sumB*table[K*u+1];
		}

		//the first pixel (x=0)
		float* q=&dst[w*y];
		q[0]=norm*(sum+a1+a2);
		d=buf[atE(0+r+1)]-buf[atW(0-r)];
		sum+=d;

		int x=1;
		while(true) //ring buffers with the length of four
		{
			q[x]=norm*(sum+b1+b2);
			e=d; d=buf[atE(x+r+1)]-buf[atW(x-r)]; delta=e-d;
			sum+=d;
			a1+=-cf11*b1+cfR1*delta;
			a2+=-cf12*b2+cfR2*delta;
			x++; if(w<=x) break;
			
			q[x]=norm*(sum-a1-a2);
			e=d; d=buf[atE(x+r+1)]-buf[atW(x-r)]; delta=e-d;
			sum+=d;
			b1+=+cf11*a1+cfR1*delta;
			b2+=+cf12*a2+cfR2*delta;
			x++; if(w<=x) break;
			
			q[x]=norm*(sum-b1-b2);
			e=d; d=buf[atE(x+r+1)]-buf[atW(x-r)]; delta=e-d;
			sum+=d;
			a1+=-cf11*b1-cfR1*delta;
			a2+=-cf12*b2-cfR2*delta;
			x++; if(w<=x) break;
			
			q[x]=norm*(sum+a1+a2);
			e=d; d=buf[atE(x+r+1)]-buf[atW(x-r)]; delta=e-d;
			sum+=d;
			b1+=+cf11*a1-cfR1*delta;
			b2+=+cf12*a2-cfR2*delta;
			x++; if(w<=x) break;
		}
	}
}
template<>
inline void gauss::filter_v<float>(int w,int h,float* src,float* dst)
{
	const int r=ry;
	const float norm=float(1.0/(r+1+r));
	std::vector<float> table(tableY.size());
	for(int t=0;t<int(table.size());++t)
		table[t]=float(tableY[t]);

	//work space to keep raster scanning
	std::vector<float> workspace((2*K+1)*w);

	//calculating the first and second terms
	for(int x=0;x<w;++x)
	{
		float* ws=&workspace[(2*K+1)*x];
		ws[0]=src[x];
		ws[1]=src[x+w]*table[0]; ws[2]=src[x]*table[0];
		ws[3]=src[x+w]*table[1]; ws[4]=src[x]*table[1];
	}
	for(int v=1;v<=r;++v)
	{
		for(int x=0;x<w;++x)
		{
			const float sum0=src[x+w*atN(0-v)]+src[x+w*(0+v)];
			const float sum1=src[x+w*atN(1-v)]+src[x+w*(1+v)];
			float* ws=&workspace[(2*K+1)*x];
			ws[0]+=sum0;
			ws[1]+=sum1*table[K*v+0]; ws[2]+=sum0*table[K*v+0];
			ws[3]+=sum1*table[K*v+1]; ws[4]+=sum0*table[K*v+1];
		}
	}

	const float cf11=float(table[K*1+0]*2.0/spectY[0]), cfR1=table[K*r+0];
	const float cf12=float(table[K*1+1]*2.0/spectY[1]), cfR2=table[K*r+1];
	
	float *q,*p0N,*p0S,*p1N,*p1S;
	for(int y=0;y<1;++y) //the first line (y=0)
	{
		q=&dst[w*y];
		p1N=&src[w*atN(0-r  )]; p1S=&src[w*atS(0+r+1)];
		for(int x=0;x<w;++x)
		{
			float* ws=&workspace[(2*K+1)*x];
			q[x]=norm*(ws[0]+ws[2]+ws[4]);
			ws[0]+=p1S[x]-p1N[x];
		}
	}
	for(int y=1;y<h;++y) //remaining lines (with two-length ring buffers)
	{
		q=&dst[w*y];
		p0N=&src[w*atN(y-r-1)]; p0S=&src[w*atS(y+r  )];
		p1N=&src[w*atN(y-r  )]; p1S=&src[w*atS(y+r+1)];
		for(int x=0;x<w;++x)
		{
			float* ws=&workspace[(2*K+1)*x];
			q[x]=norm*(ws[0]+ws[1]+ws[3]);

			const float d0=p0S[x]-p0N[x];
			const float d1=p1S[x]-p1N[x];
			const float delta=d1-d0;

			ws[0]+=d1;
			ws[2]=cfR1*delta+cf11*ws[1]-ws[2];
			ws[4]=cfR2*delta+cf12*ws[3]-ws[4];
		}
		y++; if(h<=y) break; //to the next line
		
		q=&dst[w*y];
		p0N=&src[w*atN(y-r-1)]; p0S=&src[w*atS(y+r  )];
		p1N=&src[w*atN(y-r  )]; p1S=&src[w*atS(y+r+1)];
		for(int x=0;x<w;++x)
		{
			float* ws=&workspace[(2*K+1)*x];
			q[x]=norm*(ws[0]+ws[2]+ws[4]);

			const float d0=p0S[x]-p0N[x];
			const float d1=p1S[x]-p1N[x];
			const float delta=d1-d0;

			ws[0]+=d1;
			ws[1]=cfR1*delta+cf11*ws[2]-ws[1];
			ws[3]=cfR2*delta+cf12*ws[4]-ws[3];
		}
	}
}

//*************************************************************************************************
}
