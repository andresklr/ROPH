#include "image.h"

Mat readImage(std::string fileAddress) {
	Mat image = imread(fileAddress, IMREAD_COLOR);
	Mat image_float;
	image.convertTo(image_float, CV_32F);
	return image_float;
}

float* getGrayScaleFromImage(Mat* image, int& nx, int& ny, bool padding) {
	nx = image->rows;
	ny = image->cols;

	int channels = image->channels();

	int nRows = nx;
	int nCols = ny;

	if (image->isContinuous())
	{
		nCols *= nRows;
		nRows = 1;
	}

	int sizeImage = nx * ny * sizeof(float);
	float* result = (float*)malloc(sizeImage);;

	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols; j++) {
			float* row = image->ptr<float>(i);
			float* p = (float*)malloc(channels * sizeof(float));
			for (int k = 0; k < channels; k++) {
				float* v = &row[((j)*channels) + k];
				p[k] = v[0];
			}
			float pGS = getGrayScalePixelFromBGR(p);
			result[z(i, j, ny)] = pGS;
			free(p);
		}
	}

	if (padding) {
		nx = nx + 2;
		ny = ny + 2;

		sizeImage = nx * ny * sizeof(float);

		float* padResult = (float*)malloc(sizeImage);

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
					padResult[z(i, j, ny)] = 0;
				}
				else {
					padResult[z(i, j, ny)] = result[z(i - 1, j - 1, ny - 2)];
				}
			}
		}

		return padResult;
	}

	return result;
}

int z(int i, int j, int ny) {
	return (i * ny + j);
}

float getGrayScalePixelFromBGR(float* p) {
	float result = (float)(0.2989 * p[2] + 0.5870 * p[1] + 0.1140 * p[0]);
	return result;
}

Mat grayScaleFloatToImage(Mat* image, float* in, int nx, int ny) {
	Mat result = Mat(nx, ny, image->type());

	int channels = result.channels();

	int nRows = result.rows;
	int nCols = result.cols;

	if (result.isContinuous())
	{
		nCols *= nRows;
		nRows = 1;
	}

	for (int i = 0; i < nRows; i++) {
		float* row = result.ptr<float>(i);
		for (int j = 0; j < nCols; j++) {
			for (int k = 0; k < channels; k++) {
				float* v = &row[(j * channels) + k];
				v[0] = in[z(i, j, ny)];
			}
		}
	}

	return result;
}

Mat writeCorners(Mat* image, float* out, int nx, int ny, float threshold, int typePrint, int typeFeature) {
	Mat result = image->clone();
	int cornersCount = 0;

	int channels = result.channels();

	int nRows = result.rows;
	int nCols = result.cols;

	if (result.isContinuous())
	{
		nCols *= nRows;
		nRows = 1;
	}

	for (int i = 0; i < nRows; i++) {
		float* row = result.ptr<float>(i);
		for (int j = 0; j < nCols; j++) {
			if (out[z(i, j, ny)] > threshold) {
				float* v = &row[(j * channels) + 0];
				switch (typePrint) {
				case (0):
					v[0] = 255;
					v[1] = 255;
					v[2] = 255;
					break;
				case (1):
					v[0] = 0;
					v[1] = 0;
					v[2] = 0;
					break;
				case (2):
					v[0] = 0;
					v[1] = 0;
					v[2] = 255;
					break;
				}

				cornersCount++;
			}
			else {
				for (int k = 0; k < channels; k++) {
					float* v = &row[(j * channels) + k];
					switch (typePrint) {
					case (0):
						v[0] = 0;
						v[1] = 0;
						v[2] = 0;
						break;
					case (1):
						v[0] = 255;
						v[1] = 255;
						v[2] = 255;
						break;
					}
				}
			}
		}
	}
	printMessage("Número de Corners:" + std::to_string(cornersCount));

	return result;
}

void GaussianSmoothing(float* G, float* X, int nx, int ny) {
	unsigned int i, j;
	//#pragma omp parallel
	for (i = 1; i < nx - 1; i++) {
		for (j = 1; j < ny - 1; j++) {
			G[z(i, j, ny)] = (
				X[z(i - 1, j - 1, ny)] + 2 * X[z(i - 1, j, ny)] + X[z(i - 1, j + 1, ny)] +
				2 * X[z(i, j - 1, ny)] + 4 * X[z(i, j, ny)] + 2 * X[z(i, j + 1, ny)] +
				X[z(i + 1, j - 1, ny)] + 2 * X[z(i + 1, j, ny)] + X[z(i + 1, j + 1, ny)]
				) / 16.0f;
		}
	}
}

Mat applyTransformations(const cv::Mat& inputImage, double rotationAngle, double noiseSigma, double scaleFactor) {
	
	cv::Mat transformedImage = inputImage.clone();
	cv::resize(transformedImage, transformedImage, cv::Size(), scaleFactor, scaleFactor, cv::INTER_LINEAR);

	cv::Point2f center(transformedImage.cols / 2.0f, transformedImage.rows / 2.0f);
	
	cv::Mat rotationMatrix = cv::getRotationMatrix2D(center, rotationAngle, 1.0);
	
	cv::warpAffine(transformedImage, transformedImage, rotationMatrix, transformedImage.size());

	cv::Mat noise = cv::Mat(transformedImage.size(), transformedImage.type());

	cv::randn(noise, 0.0, noiseSigma);
	
	cv::add(transformedImage, noise, transformedImage);

	return transformedImage;
}