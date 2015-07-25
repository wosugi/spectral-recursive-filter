#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <opencv2/opencv.hpp>
#include "convolution.hpp"
#include "spectral_recursive_filter.hpp"

#define CV_VERSION_NUMBER CVAUX_STR(CV_MAJOR_VERSION) CVAUX_STR(CV_MINOR_VERSION) CVAUX_STR(CV_SUBMINOR_VERSION)

#ifdef _DEBUG
#pragma comment(lib, "opencv_imgproc"CV_VERSION_NUMBER"d.lib")
#pragma comment(lib, "opencv_highgui"CV_VERSION_NUMBER"d.lib")
#pragma comment(lib, "opencv_contrib"CV_VERSION_NUMBER"d.lib")
#else
#pragma comment(lib, "opencv_imgproc"CV_VERSION_NUMBER".lib")
#pragma comment(lib, "opencv_highgui"CV_VERSION_NUMBER".lib")
#pragma comment(lib, "opencv_core"CV_VERSION_NUMBER".lib")
#endif


//Please set all the switches false when evaluating the calc time.
const bool sw_show_output=false;
const bool sw_save_output=false;

template <typename T>
double calc_snr(const cv::Mat& image1,const cv::Mat& image2,double maxval=255.0)
{
	assert(image1.size()==image2.size());
	const int w=image1.cols;
	const int h=image1.rows;
	
	double sse=0.0;
	for(int j=0;j<h;++j)
	for(int i=0;i<w;++i)
	{
		double p=image1.at<T>(j,i);
		double q=image2.at<T>(j,i);
		p=(p<0.0)?0.0:(maxval<p)?maxval:p;
		q=(q<0.0)?0.0:(maxval<q)?maxval:q;
		sse+=(q-p)*(q-p);
	}
	const double EPS=0.000001;
	if(sse<=EPS)
		return 0.0; //means the infinity

	double mse=sse/(w*h);
	double snr=-10.0*log10(mse/(maxval*maxval));
	return snr;
}

int main(int argc,char** argv)
{
	if(argc!=2)
	{
		std::cerr<<"Usage: srf [imagepath]"<<std::endl;
		return 1;
	}

	const std::string path(argv[1]);
	cv::Mat image0=cv::imread(path,0); //loaded as a grayscale image
	if(image0.empty())
	{
		std::cerr<<"input image loading failed! : "<<path<<std::endl;
		return 1;
	}

	//enabling OpenCV acceleration mode
	cv::setUseOptimized(true);
	
	//translating dynamic range [0,255] to [0,1]
	cv::Mat_<float> src=image0/255.0f;
	cv::Mat_<float> dst0(image0.size());
	cv::Mat_<float> dstA(image0.size());
	cv::Mat_<float> dstB(image0.size());

	//test various Gaussian filters over a range of [2^-0.25,2^+4.0]
	std::cerr<<cv::format("Image: (%4d,%4d) ",src.cols,src.rows)<<path<<std::endl;
	std::cerr<<"  scale :   5s-Conv    OpenCV               SRF     "<<std::endl;
	for(double octave=-0.25;octave<=4.0;octave+=0.25)
	{
		const double s=pow(2.0,octave); //scale (std. dev.) of Gaussian
		
		int64 tick0=cv::getTickCount();
		convolution::gauss cnv_gauss(s,s);
		cnv_gauss.filter(src,dst0);
		int64 tick1=cv::getTickCount();
		double time=(tick1-tick0)*1000.0/cv::getTickFrequency();
		
		int64 tickA0=cv::getTickCount();
		int sz=int(3.0*s+0.99999999);
		cv::GaussianBlur(src,dstA,cv::Size(0,0),s,0.0,cv::BORDER_REFLECT_101);
		int64 tickA1=cv::getTickCount();
		double timeA=(tickA1-tickA0)*1000.0/cv::getTickFrequency();
		double snrA=calc_snr<float>(dst0,dstA,1.0);
		
		int64 tickB0=cv::getTickCount();
		spectral_recursive_filter::gauss srf_gauss(s,s);
		srf_gauss.filter(src,dstB);
		int64 tickB1=cv::getTickCount();
		double timeB=(tickB1-tickB0)*1000.0/cv::getTickFrequency();
		double snrB=calc_snr<float>(dst0,dstB,1.0);
		
		std::string fmt="%7.4f : %7.2fms %7.2fms (%6.2fdB) %7.2fms (%6.2fdB)";
		std::cerr<<cv::format(fmt.c_str(),s,time,timeA,snrA,timeB,snrB)<<std::endl;
		
		if(sw_show_output)
		{
			cv::imshow("Input",src);
			cv::imshow("Output (conv)",dst0);
			cv::imshow("Output (ocv)",dstA);
			cv::imshow("Output (srf)",dstB);
			const float amp=50.0f;
			cv::imshow("50x amplified error (ocv)",0.5f+(dstA-dst0)*amp);
			cv::imshow("50x amplified error (srf)",0.5f+(dstB-dst0)*amp);
			cv::waitKey();
		}
		if(sw_save_output)
		{
			cv::imwrite("../ocv.png",dstA*255.0f);
			cv::imwrite("../srf.png",dstB*255.0f);
		}
	}
	return 0;
}
