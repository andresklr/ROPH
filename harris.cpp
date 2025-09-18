#include "harris.h"
#include <stdio.h>
#include <stdlib.h>
#include "image.h"

void SobelXYAndMultiplY(float* Ixx, float* Ixy, float* Iyy, float* Ix, float* Iy, float* X, int nx, int ny) {
	unsigned int i, j, index;
	//#pragma omp parallel
	for (i = 1; i < nx - 1; i++)
		for (j = 1; j < ny - 1; j++) {
			index = z(i, j, ny);
			Ix[z(i, j, ny)] = (
				-1 * X[z(i - 1, j - 1, ny)] + X[z(i - 1, j + 1, ny)] +
				-2 * X[z(i, j - 1, ny)] + 2 * X[z(i, j + 1, ny)] +
				-1 * X[z(i + 1, j - 1, ny)] + X[z(i + 1, j + 1, ny)]
				) / 8.0f;

			Iy[z(i, j, ny)] = (
				-1 * X[z(i - 1, j - 1, ny)] + -2 * X[z(i - 1, j, ny)] + -1 * X[z(i - 1, j + 1, ny)] +
				+X[z(i + 1, j - 1, ny)] + +2 * X[z(i + 1, j, ny)] + X[z(i + 1, j + 1, ny)]
				) / 8.0f;

			Ixx[index] = Ix[index] * Ix[index];
			Ixy[index] = Ix[index] * Iy[index];
			Iyy[index] = Iy[index] * Iy[index];
		}
}

void GaussAndCoarsitY(float* R, float* A, float* B, float* C, float* Ixx, float* Ixy, float* Iyy, int nx, int ny, float k) {
	unsigned int i, j;
	//#pragma omp parallel
	for (i = 1; i < nx - 1; i++)
		for (j = 1; j < ny - 1; j++)
		{
			A[z(i, j, ny)] = (
				Ixx[z(i - 1, j - 1, ny)] + 2 * Ixx[z(i - 1, j, ny)] + Ixx[z(i - 1, j + 1, ny)] +
				2 * Ixx[z(i, j - 1, ny)] + 4 * Ixx[z(i, j, ny)] + 2 * Ixx[z(i, j + 1, ny)] +
				Ixx[z(i + 1, j - 1, ny)] + 2 * Ixx[z(i + 1, j, ny)] + Ixx[z(i + 1, j + 1, ny)]
				) / 16.0f;

			B[z(i, j, ny)] = (
				Ixy[z(i - 1, j - 1, ny)] + 2 * Ixy[z(i - 1, j, ny)] + Ixy[z(i - 1, j + 1, ny)] +
				2 * Ixy[z(i, j - 1, ny)] + 4 * Ixy[z(i, j, ny)] + 2 * Ixy[z(i, j + 1, ny)] +
				Ixy[z(i + 1, j - 1, ny)] + 2 * Ixy[z(i + 1, j, ny)] + Ixy[z(i + 1, j + 1, ny)]
				) / 16.0f;

			C[z(i, j, ny)] = (
				Iyy[z(i - 1, j - 1, ny)] + 2 * Iyy[z(i - 1, j, ny)] + Iyy[z(i - 1, j + 1, ny)] +
				2 * Iyy[z(i, j - 1, ny)] + 4 * Iyy[z(i, j, ny)] + 2 * Iyy[z(i, j + 1, ny)] +
				Iyy[z(i + 1, j - 1, ny)] + 2 * Iyy[z(i + 1, j, ny)] + Iyy[z(i + 1, j + 1, ny)]
				) / 16.0f;

			R[z(i, j, ny)] = A[z(i, j, ny)] * C[z(i, j, ny)] - B[z(i, j, ny)] * B[z(i, j, ny)]
				- (k * (A[z(i, j, ny)] * A[z(i, j, ny)] + 2 * A[z(i, j, ny)] * C[z(i, j, ny)] + C[z(i, j, ny)] * C[z(i, j, ny)]));
		}
}

void Harris(int nx, int ny, float* I, float* R, float k) {
	int sizeImage = nx * ny * sizeof(float);
	float* Ix = (float*)malloc(sizeImage);;
	float* Iy = (float*)malloc(sizeImage);;
	float* Ixx = (float*)malloc(sizeImage);;
	float* Iyy = (float*)malloc(sizeImage);;
	float* Ixy = (float*)malloc(sizeImage);;
	float* A = (float*)malloc(sizeImage);;
	float* B = (float*)malloc(sizeImage);;
	float* C = (float*)malloc(sizeImage);;

	//Sobel and MultiplY
	SobelXYAndMultiplY(Ixx, Ixy, Iyy, Ix, Iy, I, nx, ny);

	//Gauss and CoarsitY
	GaussAndCoarsitY(R, A, B, C, Ixx, Ixy, Iyy, nx, ny, k);

	free(Ix);
	free(Iy);
	free(Ixx);
	free(Iyy);
	free(Ixy);
	free(A);
	free(B);
	free(C);
}

