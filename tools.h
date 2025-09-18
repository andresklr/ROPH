#pragma once
#include <string>
#include <iostream>
#include <vector>
#include "omp.h"

struct FastCandidates {
    int** data;
    int* dataSize;
    int num_sections;
};

void convertRValues(float* out, float* out2, int nx, int ny, int typeFeature);

void printMessage(std::string msg);

double getTime(time_t* start, time_t* end);

void countScores(int nx, int ny, float* out, float threshold);

void printArray(int nx, int ny, float* out);

void identifyCorners(float* out, int nx, int ny, float threshold);

float* getGroundTruth(std::string groundTruthFile);