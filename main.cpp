#include"FundFunction.h"





int main() {
	Mat M = imread("k.bmp");
	imshow("basic", M);
	imshow("Deal", PsColorEnhanceInFrequencyDomain(M));
	waitKey();
	return 0;
}