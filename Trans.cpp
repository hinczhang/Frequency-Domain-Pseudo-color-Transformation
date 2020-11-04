#include"FundFunction.h"

void FastFourie(complex<double>* pSData, complex<double>* pFData, int level)
{
	int count = (int)pow(2, level);//����Ҷ�任����

	complex<double>* pW = new complex<double>[count / 2];//�����ڴ棬�洢����Ҷ�任��Ҫ��ϵ����

	for (int i = 0; i < count / 2; i++)//���㸵��Ҷ�任��ϵ��
	{
		double angle = -2 * PI*i / count;
		pW[i] = complex<double>(cos(angle), sin(angle));
	}

	/*���乤����Ҫ���ڴ�ռ�*/
	complex<double>* pWork1 = new complex<double>[count];
	complex<double>* pWork2 = new complex<double>[count];

	memcpy(pWork1, pSData, sizeof(complex<double>)*count);//��ʼ����д������

	/*�����㷨���п��ٸ���Ҷ�任*/
	for (int k = 0; k < level; k++)
	{
		for (int i = 0; i < (int)pow(2, k); i++)
		{
			int BtFlyLen = (int)pow(2, level - k);//���㳤��
			int Inter = i * BtFlyLen;
			/*�������ţ���Ȩ����*/
			for (int j = 0; j < BtFlyLen / 2; j++)
			{
				pWork2[j + Inter] = pWork1[j + Inter] + pWork1[j + Inter + BtFlyLen / 2];
				pWork2[j + Inter + BtFlyLen / 2] = (pWork1[j + Inter] - pWork1[j + Inter + BtFlyLen / 2])*pW[j*(int)pow(2, k)];
			}
		}
		/*����pWork1��pWork2������*/
		complex<double>* pTmp = pWork1;
		pWork1 = pWork2;
		pWork2 = pTmp;
	}

	/*��������*/
	for (int i = 0; i < count; i++)
	{
		int Inter = 0;
		for (int j = 0; j < level; j++)
		{
			if (i&(1 << j)) Inter += 1 << (level - j - 1);
		}
		pFData[i] = pWork1[Inter];
	}

	/*�ͷ��ڴ�*/
	delete[] pW; pW = NULL;
	delete[] pWork1; pWork1 = NULL;
	delete[]pWork2; pWork2 = NULL;
}
void Fourie(complex<double>* pSData, int width_T, int height_T, complex<double>* pFData)
{

	/*����x,y�����У������ϵĵ�������*/
	int Xlevel = (int)(log(width_T) / log(2) + 0.5);
	int Ylevel = (int)(log(height_T) / log(2) + 0.5);

	/*�з�����п��ٸ���Ҷ�任*/
	for (int i = 0; i < height_T; i++)
	{
		FastFourie(&pSData[width_T*i], &pFData[width_T*i], Xlevel);
	}

	/*ת��*/
	for (int i = 0; i < height_T; i++)
	{
		for (int j = 0; j < width_T; j++)
		{
			pSData[height_T*j + i] = pFData[width_T*i + j];
		}
	}

	/*�з�����п��ٸ���Ҷ�任*/
	for (int i = 0; i < height_T; i++)
	{
		FastFourie(&pSData[height_T*i], &pFData[height_T*i], Ylevel);
	}

	/*ת�û���*/
	for (int i = 0; i < height_T; i++)
	{
		for (int j = 0; j < width_T; j++)
		{
			pSData[width_T*i + j] = pFData[height_T*j + i];
		}
	}
	memcpy(pFData, pSData, sizeof(complex<double>)*height_T*width_T);
}
void IFourie(complex<double>* pFData, int width_T, int height_T, complex<double>* pSData)
{

	/*���乤����Ҫ���ڴ�ռ�*/
	complex<double>* pWork = new complex<double>[width_T*height_T];

	/*��Ƶ�������ȡ����*/
	for (int i = 0; i < height_T; i++)
	{
		for (int j = 0; j < width_T; j++)
		{
			complex<double>* pTmp = &pFData[width_T*i + j];
			pWork[width_T*i + j] = complex<double>(pTmp->real(), -pTmp->imag());
		}
	}

	/*���ø���Ҷ�任*/
	Fourie(pWork, width_T, height_T, pSData);

	/*�����Ĺ�����ý��*/
	for (int i = 0; i < height_T; i++)
	{
		for (int j = 0; j < width_T; j++)
		{
			complex<double>* pTmp = &pSData[width_T*i + j];
			pSData[width_T*i + j] = complex<double>(pTmp->real() / (width_T*height_T), -pTmp->imag() / (width_T*height_T));
		}
	}
	delete[] pWork; pWork = NULL;
}

void ButterWorthLow(complex<double>*pFData,complex<double>* pF1Data , complex<double>* pS1Data,int width_T,int height_T,int D01) {
	for (int i = 0; i < height_T; i++)
	{
		for (int j = 0; j < width_T; j++)
		{
			double H = i * i + j * j;//�˲�ϵ��
			H = H / D01 * D01;
			H = 1 / (1 + H);
			pF1Data[i*width_T + j] = complex<double>(H*pFData[i*width_T + j].real(), H*pFData[i*width_T + j].imag());
		}
	}
	IFourie(pF1Data, width_T, height_T, pS1Data);
}

void ButterWorthHigh(complex<double>*pFData, complex<double>* pF2Data, complex<double>* pS2Data, int width_T, int height_T, int D02) {
	for (int i = 0; i < height_T; i++)
	{
		for (int j = 0; j < width_T; j++)
		{
			double H = i * i + j * j;//�˲�ϵ��
			H = D02 * D02 / H;
			H = 1 / (1 + H);
			pF2Data[i*width_T + j] = complex<double>(H*pFData[i*width_T + j].real(), H*pFData[i*width_T + j].imag());
		}
	}
	IFourie(pF2Data, width_T, height_T, pS2Data);
}

void BandpassFiltering(complex<double>*pFData, complex<double>* pF3Data, complex<double>* pS3Data, int width_T, int height_T, int D03, int W) {
	for (int i = 0; i < height_T; i++)
	{
		for (int j = 0; j < width_T; j++)
		{
			double D_2 = sqrt(i*i + j * j);
			double H = D_2 * W*W / ((D_2 - D03 * D03)*(D_2 - D03 * D03));
			H = 1 / (1 + H);
			pF3Data[i*width_T + j] = complex<double>(H*pFData[i*width_T + j].real(), H*pFData[i*width_T + j].imag());
		}
	}
	IFourie(pF3Data, width_T, height_T, pS3Data);
}

Mat PsColorEnhanceInFrequencyDomain(Mat img)
{
	int width = img.cols;
	int height = img.rows;
	int cns = img.channels();
	unsigned char* pData = img.data;

	/*���㸵��Ҷ�任�Ŀ�ȣ�2���������ݣ�*/
	int width_T = (int)pow(2, ceil(log((double)width) / log(2)));

	/*���㸵��Ҷ�任�ĸ߶ȣ�2���������ݣ�*/
	int height_T = (int)pow(2, ceil(log((double)height) / log(2)));

	complex<double>* pSData = new complex<double>[width_T*height_T];
	complex<double>* pFData = new complex<double>[width_T*height_T];
	for (int i = 0; i < height_T*width_T; i++)
	{
		pSData[i] = complex<double>(0, 0);
	}
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			pSData[i*width_T + j] = complex<double>(pData[(i*width + j)*cns], 0);
		}
	}
	Fourie(pSData, width_T, height_T, pFData);


	/*��ͨ�˲�*/
	double D01 = 150;//��ֹƵ��
	complex<double>* pF1Data = new complex<double>[width_T*height_T];
	complex<double>* pS1Data = new complex<double>[width_T*height_T];
	ButterWorthLow(pFData, pF1Data, pS1Data, width_T, height_T, D01);
	/*��ͨ�˲�*/
	double D02 = 250;//��ֹƵ��
	complex<double>* pF2Data = new complex<double>[width_T*height_T];
	complex<double>* pS2Data = new complex<double>[width_T*height_T];
	ButterWorthHigh(pFData, pF2Data, pS2Data, width_T, height_T, D02);
	/*��ͨ�˲�*/
	double D03 = 200;//��������Ƶ��
	double W = 100;//����
	complex<double>* pF3Data = new complex<double>[width_T*height_T];
	complex<double>* pS3Data = new complex<double>[width_T*height_T];
	BandpassFiltering(pFData, pF3Data, pS3Data, width_T, height_T, D03, W);

	Mat imgR, imgG, imgB;
	imgR.create(height, width, CV_8UC1);
	imgG.create(height, width, CV_8UC1);
	imgB.create(height, width, CV_8UC1);
	unsigned char* pR = imgR.data;
	unsigned char* pG = imgG.data;
	unsigned char* pB = imgB.data;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			double real = pS1Data[i*width_T + j].real();
			double imag = pS1Data[i*width_T + j].imag();
			pR[i*width + j] = (unsigned char)max(0, min(255, sqrt(real*real + imag * imag) + 100));//r

			real = pS2Data[i*width_T + j].real();
			imag = pS2Data[i*width_T + j].imag();
			pG[i*width + j] = (unsigned char)max(0, min(255, sqrt(real*real + imag * imag) + 100));//g

			real = pS3Data[i*width_T + j].real();
			imag = pS3Data[i*width_T + j].imag();
			pB[i*width + j] = (unsigned char)max(0, min(255, sqrt(real*real + imag * imag) + 100));//b
		}
	}
	EqualizeHistogram(imgR, imgR);
	EqualizeHistogram(imgG, imgG);
	EqualizeHistogram(imgB, imgB);
	Mat imgN;
	imgN.create(height, width, CV_8UC3);
	unsigned char* pNew = imgN.data;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			pNew[(i*width + j) * 3 + 0] = pB[i*width + j];
			pNew[(i*width + j) * 3 + 1] = pG[i*width + j];
			pNew[(i*width + j) * 3 + 2] = pR[i*width + j];
		}
	}
	delete[] pSData; pSData = NULL;
	delete[] pFData; pFData = NULL;
	delete[] pF1Data; pF1Data = NULL;
	delete[] pF2Data; pF2Data = NULL;
	delete[] pF3Data; pF3Data = NULL;
	return imgN;
}
void GetHistogram(const Mat &image, int *histogram)
{
	memset(histogram, 0, 256 * sizeof(int));

	//����ֱ��ͼ
	int pixelCount = image.cols*image.rows;
	uchar *imageData = image.data;
	for (int i = 0; i <= pixelCount - 1; ++i)
	{
		int gray = imageData[i];
		histogram[gray]++;
	}
}

void EqualizeHistogram(const Mat &srcImage, Mat &dstImage)
{
	// ����ֱ��ͼ
	int histogram[256];
	GetHistogram(srcImage, histogram);

	// ����ֲ�����(Ҳ���Ǳ任����f(x))
	int numberOfPixel = srcImage.rows*srcImage.cols;
	int LUT[256];
	LUT[0] = 1.0*histogram[0] / numberOfPixel * 255;
	int sum = histogram[0];
	for (int i = 1; i <= 255; ++i)
	{
		sum += histogram[i];

		LUT[i] = 1.0*sum / numberOfPixel * 255;
	}

	// �Ҷȱ任
	unsigned char *dataOfSrc = new unsigned char[srcImage.cols*srcImage.rows];
	memcpy(dataOfSrc, srcImage.data, sizeof(unsigned char)*srcImage.cols*srcImage.rows);
	unsigned char *dataOfDst = dstImage.data;
	for (int i = 0; i <= numberOfPixel - 1; ++i)
	{
		dataOfDst[i] = LUT[dataOfSrc[i]];
	}
	delete[] dataOfSrc; dataOfSrc = NULL;
}