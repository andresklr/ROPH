#pragma once
#define _HARRIS_H
#include <stdio.h>
#include <stdlib.h>
#include "image.h"

void HarrisFast(int nx, int ny, float* I, float* R, float k, FastCandidates& candidates, float* candidates3, bool padding, int num_harris_threads);

void SobelXYAndMultiply(float* Ixx, float* Ixy, float* Iyy, float* Ix, float* Iy, float* I, int nx, int ny, float* candidates, bool padding, int num_harris_threads);

void GaussAndCoarsitY(float* R, float* A, float* B, float* C, float* Ixx, float* Ixy, float* Iyy, int nx, int ny, float k, FastCandidates& candidates, int num_harris_threads);