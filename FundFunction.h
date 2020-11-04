#ifndef FUND_H_
#define FUND_H_
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include <windows.h>
#include <math.h>

#define PI 3.1415926
using namespace cv;
using namespace std;


void ButterWorthLow(complex<double>*pFData, complex<double>* pF1Data, complex<double>* pS1Data, int width_T, int height_T, int D01);
void ButterWorthHigh(complex<double>*pFData, complex<double>* pF2Data, complex<double>* pS2Data, int width_T, int height_T, int D02);
void BandpassFiltering(complex<double>*pFData, complex<double>* pF3Data, complex<double>* pS3Data, int width_T, int height_T, int D03, int W);
void FastFourie(complex<double>* pSData, complex<double>* pFData, int level);//快速傅里叶变换
void Fourie(complex<double>* pSData, int width, int height, complex<double>* pFData);//傅里叶变换
void IFourie(complex<double>* pFData, int width, int height, complex<double>* pSData);//傅里叶反变换
Mat  PsColorEnhanceInFrequencyDomain(Mat img);//频域
void GetHistogram(const Mat &image, int *histogram);//计算直方图
void EqualizeHistogram(const Mat &srcImage, Mat &dstImage);//直方图均衡化


#endif