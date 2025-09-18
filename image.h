#pragma once
#include <opencv.hpp>
#include <iostream>
#include "tools.h"

using namespace cv;

struct Corner {
	float r;
	int i;
	int j;
	bool skipped;
};

Mat readImage(std::string fileAddress);

float* getGrayScaleFromImage(Mat* image, int& nx, int& ny, bool padding);

int z(int i, int j, int ny);

Mat grayScaleFloatToImage(Mat* image, float* in, int nx, int ny);

Mat writeCorners(Mat* image, float* out, int nx, int ny, float threshold, int typePrint, int typeFeature);

float getGrayScalePixelFromBGR(float* p);

void GaussianSmoothing(float* G, float* X, int nx, int ny);

Mat applyTransformations(const cv::Mat& inputImage, double rotationAngle, double noiseSigma, double scaleFactor);

