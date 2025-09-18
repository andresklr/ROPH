#pragma once
#define _HARRIS_H

void SobelXYAndMultiplY(float* Ixx, float* Ixy, float* Iyy, float* Ix, float* Iy, float* X, int nx, int ny);

void GaussAndCoarsitY(float* R, float* A, float* B, float* C, float* Ixx, float* Ixy, float* Iyy, int nx, int ny, float k);

void Harris(int nx, int ny, float* I, float* R, float k);