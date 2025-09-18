// HarrisFast.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#include <chrono>
#include <iostream>
#include <fstream>
#include "harris.h"
#include "tools.h"
#include "image.h"
#include "fast.h"
#include "nms.h"
#include "harrisFast.h"
#include <filesystem>

using namespace std;

void main() {

	float fastThreshold = 15.0f;
	int nmsWindowSize = 3;
	//int maxCornersNMS = 6;
	float k = 0.04;
	float threshold = 1000.0f;
	float groundTruthThreshold = threshold;
	int typeExecution = 2; //0 - Harris, 1 - Fast+Harris, 2 - Both, 3 - Inverted (NMS vs Harris)
	int typeFeature = 1; // 0 - Corners, 1 - Edges, 2 - Both
	bool nms = false;
	int typePrint = 0; // 0 - Black background white features, 1 - Reversed, 2 - Black and white with red Features
	int typeCorner = 0; // 0 - Point, 1 - Cross, 2 - Circle
	int fastThreads = 5;
	int harrisThreads = 5;
	bool padding = true;
	bool preProcessing = true;
	string textNMS = nms ? "-nms-" : "";
	int numberOfExecutions = 1;
	int numberOfExecutionsNMS = 1;
	string file = "r";
	string fileFormat = ".jpg";
	string basePath = "C:\\GitHub\\Personal\\UERJME\\Dissertação\\Code\\HarrisFast\\data\\results4\\";
	//string imageName = basePath + "images\\" + file;
	string imageName = basePath + "\\" + file;
	string feature = (typeFeature == 0 ? "corners" : (typeFeature == 1 ? "edges" : "both"));
	string methodDescription = feature + "-ft-" + to_string((int)fastThreshold) + "-t-" + to_string((int)threshold);

	std::ostringstream oss;

	oss << "--- Run Parameters ---\n"
		<< "FAST Threshold: " << fastThreshold << "\n"
		<< "NMS Enabled: " << std::boolalpha << nms << "\n"
		<< "NMS Window: " << nmsWindowSize << "\n"
		<< "Harris k: " << k << "\n"
		<< "Harris Threshold: " << threshold << "\n"
		<< "File: " << file << "\n"
		// --- MODIFIED LINES START HERE ---
		<< "Execution Type: " << (typeExecution == 0 ? "Harris" :
			(typeExecution == 1 ? "Fast+Harris" :
				(typeExecution == 3 ? "Inverted (NMS vs Harris)" : "Both"))) << "\n"
		<< "Feature Type: " << (typeFeature == 0 ? "Corners" :
			(typeFeature == 1 ? "Edges" : "Both")) << "\n"
		<< "Print Type: " << (typePrint == 0 ? "Black background, white features" :
			(typePrint == 1 ? "Reversed" : "Black and white with red features")) << "\n"
		<< "Corner Type: " << (typeCorner == 0 ? "Point" :
			(typeCorner == 1 ? "Cross" : "Circle")) << "\n"
		// --- MODIFIED LINES END HERE ---
		<< "FAST Threads: " << fastThreads << "\n"
		<< "Harris Threads: " << harrisThreads << "\n"
		<< "Padding: " << std::boolalpha << padding << "\n"
		<< "Preprocessing: " << std::boolalpha << preProcessing << "\n"
		<< "----------------------";

	std::string run_parameters_string = oss.str();

	/*ofstream log_file(imageName + "-" + methodDescription + ".txt", std::ios::out | std::ios::trunc);
	streambuf* original_cout_buffer = std::cout.rdbuf();
	cout.rdbuf(log_file.rdbuf());*/

	printMessage(run_parameters_string);


	Mat image = readImage(imageName + fileFormat);

	//Mat imageTransformed = applyTransformations(image, 15, 10, 0.9);
	//
	//imwrite(imageName + "-transformed-" + fileFormat, imageTransformed);

	int nx, ny;

	float* in = getGrayScaleFromImage(&image, nx, ny, padding);

	Mat image2 = grayScaleFloatToImage(&image, in, nx, ny);

	//imwrite(imageName + "-grayscale" + fileFormat, image2);

	int sizeImage = nx * ny * sizeof(float);
	float* I = (float*)malloc(sizeImage);
	float* out = (float*)calloc(nx * ny, sizeof(float));
	float* out2 = (float*)calloc(nx * ny, sizeof(float));
	float* out3 = (float*)calloc(nx * ny, sizeof(float));
	FastCandidates candidates;
	std::vector<Corner> corners;
	std::vector<Corner> corners2;

	auto time_Start = std::chrono::high_resolution_clock::now();

	if (preProcessing) {
		GaussianSmoothing(I, in, nx, ny);
	}
	else {
		I = in;
	}


	auto time_End = std::chrono::high_resolution_clock::now();

	auto time_Pre = time_End - time_Start;

	time_Start = std::chrono::high_resolution_clock::now();

	auto time_FAST = std::chrono::high_resolution_clock::duration::zero();
	auto time_HarrisF = time_FAST;
	auto time_First = time_Start;

	for (int i = 0; i < numberOfExecutions; i++) {
		if (typeExecution == 3) {
			Harris(nx, ny, I, out2, k);
		}
		else {
			if (typeExecution > 0) {
				if (typeFeature == 0) {
					getFast4CandidatesCorners(I, out3, out2, candidates, nx, ny, fastThreshold, fastThreads);
				}
				else {
					time_Start = std::chrono::high_resolution_clock::now();
					getFast4CandidatesEdges(I, out3, out2, candidates, nx, ny, fastThreshold, fastThreads);
					time_End = std::chrono::high_resolution_clock::now();
					time_FAST += (time_End - time_Start);
				}
				time_Start = std::chrono::high_resolution_clock::now();
				HarrisFast(nx, ny, I, out2, k, candidates, out3, padding, harrisThreads);
				time_End = std::chrono::high_resolution_clock::now();
				time_HarrisF += (time_End - time_Start);
				//countScores(nx, ny, out2, &threshold);
			}
		}
	}

	time_End = std::chrono::high_resolution_clock::now();

	auto time_FastHarris = time_End - time_First;

	//time_Start = std::chrono::high_resolution_clock::now();
	//
	//if (nms) {
	//	corners2 = NonMaximumSuppression(out2, nx, ny, threshold, nmsWindowSize);
	//}
	//
	//time_End = std::chrono::high_resolution_clock::now();
	//
	//auto time_NMS2 = time_End - time_Start;

	time_Start = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < numberOfExecutions; i++) {
		if (typeExecution != 1) {
			Harris(nx, ny, I, out, k);
		}
	}

	time_End = std::chrono::high_resolution_clock::now();

	auto time_Harris = time_End - time_Start;

	time_Start = std::chrono::high_resolution_clock::now();

	if (nms) {
		for (int i = 0; i < numberOfExecutionsNMS; i++) {
			corners = NonMaximumSuppression(out, nx, ny, threshold, nmsWindowSize);
		}
	}

	time_End = std::chrono::high_resolution_clock::now();

	auto time_NMS = time_End - time_Start;


	float speedUp = ((float)(time_Harris / (std::chrono::milliseconds(1)))) / ((float)(time_FastHarris / std::chrono::milliseconds(1)));

	if (nms) {
		printMessage("tiempo de duración Harris (ms): " + to_string((time_Harris + time_NMS) / std::chrono::milliseconds(1)));
	}
	else {
		printMessage("tiempo de duración Harris (ms): " + to_string(time_Harris / std::chrono::milliseconds(1)));
	}

	printMessage("tiempo de duración FAST+Harris (FAST) (ms): " + to_string(time_FAST / std::chrono::milliseconds(1)));

	printMessage("tiempo de duración FAST+Harris (Harris) (ms): " + to_string(time_HarrisF / std::chrono::milliseconds(1)));

	printMessage("tiempo de duración FAST+Harris (Total) (ms): " + to_string(time_FastHarris / std::chrono::milliseconds(1)));

	printMessage("tiempo de duración NMS (ms): " + to_string(time_NMS / std::chrono::milliseconds(1)));

	printMessage("Speedup: " + to_string(speedUp));

	printMessage(to_string(time_Harris / std::chrono::milliseconds(1)) + ";" +
		to_string(time_FastHarris / std::chrono::milliseconds(1)) + ";" + to_string(speedUp) + ";" + file);

	convertRValues(out, out2, nx, ny, typeFeature);

	float* groundTruth;

	groundTruth = out;

	float FPH1 = 0, TPH1 = 0, TNH1 = 0, FNH1 = 0, FPH2 = 0, TPH2 = 0, TNH2 = 0, FNH2 = 0;
	float NOP = 0, NOF1 = 0, NOF2 = 0, POF1 = 0, POF2 = 0;

	Mat result = Mat(nx, ny, image2.type());

	int nRows = result.rows;
	int nCols = result.cols;

	if (result.isContinuous())
	{
		nCols *= nRows;
		nRows = 1;
	}

	NOP = nx * ny;

	for (int i = 0; i < nRows; i++)
	{
		for (int j = 0; j < nCols; j++)
		{
			if (groundTruth[z(i, j, ny)] > groundTruthThreshold) {
				if (out[z(i, j, ny)] > threshold) {
					TPH1++;
					NOF1++;
				}
				else {
					FNH1++;
				}

				if (out2[z(i, j, ny)] > threshold) {
					TPH2++;
					NOF2++;
				}
				else {
					FNH2++;
				}
			}
			else {
				if (out[z(i, j, ny)] > threshold) {
					FPH1++;
					NOF1++;
				}
				else {
					TNH1++;
				}

				if (out2[z(i, j, ny)] > threshold) {
					FPH2++;
					NOF2++;
				}
				else {
					TNH2++;
				}
			}
		}
	}

	printMessage("TPH1: " + to_string(TPH1) + " FPH1: " + to_string(FPH1) + " TNH1: " + to_string(TNH1) + " FNH1: " + to_string(FNH1));

	string confusionMatrix = "TPH2: " + to_string(TPH2) + " FPH2: " + to_string(FPH2) + " TNH2: " + to_string(TNH2) + " FNH2: " + to_string(FNH2);

	printMessage(confusionMatrix);

	float accuracyH1 = 0, precisionH1 = 0, recallH1 = 0, f1H1 = 0, accuracyH2 = 0, precisionH2 = 0, recallH2 = 0, f1H2 = 0;

	accuracyH1 = (TPH1 + TNH1) / (TPH1 + TNH1 + FPH1 + FNH1);

	precisionH1 = (TPH1) / (TPH1 + FPH1);

	recallH1 = (TPH1) / (TPH1 + FNH1);

	f1H1 = (2 * TPH1) / ((2 * TPH1) + FPH1 + FNH1);

	printMessage("accuracyH1: " + to_string(accuracyH1) +
		" precisionH1: " + to_string(precisionH1) + " recallH1: " + to_string(recallH1) +
		" f1H1: " + to_string(f1H1));

	accuracyH2 = (TPH2 + TNH2) / (TPH2 + TNH2 + FPH2 + FNH2);

	precisionH2 = (TPH2) / (TPH2 + FPH2);

	recallH2 = (TPH2) / (TPH2 + FNH2);

	f1H2 = (2 * TPH2) / ((2 * TPH2) + FPH2 + FNH2);

	string performanceMetrics = "accuracyH2: " + to_string(accuracyH2) +
		" precisionH2: " + to_string(precisionH2) + " recallH2: " + to_string(recallH2) +
		" f1H2: " + to_string(f1H2);

	printMessage(performanceMetrics);

	POF1 = NOF1 / NOP;

	POF2 = NOF2 / NOP;

	printMessage(to_string(POF1) + ";" + to_string(POF2));

	bool multiple = false;

	if (multiple) {
		for (float i = 1.0f; i <= 16; i++) {
			float t = threshold * (i / 4);

			Mat image3 = writeCorners(&image2, out, nx, ny, t, typePrint, typeFeature);

			imwrite(imageName + "-t-" + to_string(t) + fileFormat, image3);
		}
	}
	else {
		float t2 = groundTruthThreshold;
		float t3 = 1.0f;
		Mat image3;

		switch (typeExecution) {
		case(0):
			image3 = writeCorners(&image2, out, nx, ny, threshold, typePrint, typeFeature);

			imwrite(imageName + "-Harris-" + methodDescription + textNMS + fileFormat, image3);
			break;
		case(1):
			image3 = writeCorners(&image2, out2, nx, ny, threshold, typePrint, typeFeature);

			imwrite(imageName + "-FastAndHarris-" + methodDescription + textNMS + fileFormat, image3);
			break;
		case(2):
			image3 = writeCorners(&image2, out, nx, ny, threshold, typePrint, typeFeature);

			imwrite(imageName + "-Harris-" + methodDescription + textNMS + fileFormat, image3);

			image3 = writeCorners(&image2, out2, nx, ny, threshold, typePrint, typeFeature);

			imwrite(imageName + "-FastAndHarris-" + methodDescription + textNMS + fileFormat, image3);
			break;
		case(3):
			image3 = writeCorners(&image2, out, nx, ny, threshold, typePrint, typeFeature);

			imwrite(imageName + "-Harris-" + methodDescription + textNMS + fileFormat, image3);

			image3 = writeCorners(&image2, out2, nx, ny, threshold, typePrint, typeFeature);

			imwrite(imageName + "-Harris-" + methodDescription + fileFormat, image3);
			break;
		}
	}

	free(I);
	free(out);
	free(out2);
	free(out3);

	/*cout.rdbuf(original_cout_buffer);
	log_file.close();*/

	//identifyCorners(out, nx, ny, &threshold);

}
