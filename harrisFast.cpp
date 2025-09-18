#include "harrisFast.h"

void HarrisFast(int nx, int ny, float* I, float* R, float k, FastCandidates& candidates, float* candidates3, bool padding, int num_harris_threads) {

	int sizeImage = nx * ny * sizeof(float);
	float* Ix = (float*)malloc(sizeImage);
	float* Iy = (float*)malloc(sizeImage);
	float* Ixx = (float*)malloc(sizeImage);
	float* Iyy = (float*)malloc(sizeImage);
	float* Ixy = (float*)malloc(sizeImage);
	float* A = (float*)malloc(sizeImage);
	float* B = (float*)malloc(sizeImage);
	float* C = (float*)malloc(sizeImage);

	//Sobel and Multiply
	SobelXYAndMultiply(Ixx, Ixy, Iyy, Ix, Iy, I, nx, ny, candidates3, padding, num_harris_threads);

	//Gauss and CoarsitY
	GaussAndCoarsitY(R, A, B, C, Ixx, Ixy, Iyy, nx, ny, k, candidates, num_harris_threads);

	free(candidates.data);
	free(Ix);
	free(Iy);
	free(Ixx);
	free(Iyy);
	free(Ixy);
	free(A);
	free(B);
	free(C);
}

void SobelXYAndMultiply(float* Ixx, float* Ixy, float* Iyy, float* Ix, float* Iy, float* I, int nx, int ny, float* candidates, bool padding, int num_harris_threads) {
	int upperLimit = padding ? 2 : 0;
	int lowerLimitX = padding ? nx - 2 : nx;
	int lowerLimitY = padding ? ny - 2 : ny;

#pragma omp parallel for num_threads(num_harris_threads)
	for (int i = upperLimit; i < lowerLimitX; i++) {
		for (int j = upperLimit; j < lowerLimitY; j++) {
			int index = z(i, j, ny);
			if (candidates[index] > 0) {
				Ix[index] = (
					-1 * I[index - ny - 1] + I[index - ny + 1]
					- 2 * I[index - 1] + 2 * I[index + 1]
					- 1 * I[index + ny - 1] + I[index + ny + 1]
					) / 8.0f;

				Iy[index] = (
					-1 * I[index - ny - 1] - 2 * I[index + -ny] - I[index - ny + 1] +
					I[index + ny - 1] + 2 * I[index + ny] + I[index + ny + 1]
					) / 8.0f;

				Ixx[index] = Ix[index] * Ix[index];
				Iyy[index] = Iy[index] * Iy[index];
				Ixy[index] = Ix[index] * Iy[index];
			}
		}
	}
}

void GaussAndCoarsitY(float* R, float* A, float* B, float* C, float* Ixx, float* Ixy, float* Iyy, int nx, int ny, float k, FastCandidates& candidates, int num_harris_threads) {
	int index;
#pragma omp parallel for num_threads(num_harris_threads) schedule(dynamic) private(index)
	for (int section_index = 0; section_index < candidates.num_sections; section_index++) {
		for (int i = 0; i < candidates.dataSize[section_index]; i++) {
			index = candidates.data[section_index][i];

			A[index] = (
				Ixx[index - ny - 1] + 2 * Ixx[index - ny] + Ixx[index - ny + 1] +
				2 * Ixx[index - 1] + 4 * Ixx[index] + 2 * Ixx[index + 1] +
				Ixx[index + ny - 1] + 2 * Ixx[index + ny] + Ixx[index + ny + 1]
				) / 16.0f;

			C[index] = (
				Iyy[index - ny - 1] + 2 * Iyy[index - ny] + Iyy[index - ny + 1] +
				2 * Iyy[index - 1] + 4 * Iyy[index] + 2 * Iyy[index + 1] +
				Iyy[index + ny - 1] + 2 * Iyy[index + ny] + Iyy[index + ny + 1]
				) / 16.0f;

			B[index] = (
				Ixy[index - ny - 1] + 2 * Ixy[index - ny] + Ixy[index - ny + 1] +
				2 * Ixy[index - 1] + 4 * Ixy[index] + 2 * Ixy[index + 1] +
				Ixy[index + ny - 1] + 2 * Ixy[index + ny] + Ixy[index + ny + 1]
				) / 16.0f;

			R[index] = A[index] * C[index] - B[index] * B[index]
				- (k * (A[index] * A[index] + 2 * A[index] * C[index] + C[index] * C[index]));
		}
	}
}