#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <opencv2/opencv.hpp>

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

	cv::imshow("Input",image0);
	cv::waitKey();
	return 0;
}
