#include "fast.h"

const int radius = 3;

void getFast4CandidatesCorners(float* in, float* out, float* out2, FastCandidates& results, int nx, int ny, float threshold, int num_fast_threads) {

	int num_sections = num_fast_threads + 2; // +2 for top and bottom

	results.data = (int**)malloc(sizeof(int*) * num_sections);

	results.dataSize = (int*)malloc(sizeof(int) * num_sections);

	results.num_sections = num_sections;

	size_t max_candidates_per_list = (nx * ny) / num_sections;

	for (int i = 0; i < num_sections; i++) {
		results.data[i] = (int*)malloc(sizeof(int) * max_candidates_per_list);
		results.dataSize[i] = 0;
	}

	auto time_Start = std::chrono::high_resolution_clock::now();

	int i, j, zi;
	float p, p1, p5, p9, p13;

	//Only p5 and p9
	for (i = 1; i < radius; i++) {
		for (j = 1; j < radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP5andP9Corners(&p, &p5, &p9, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p5,p9 and p13
	for (i = 1; i < radius; i++) {
		for (j = radius; j < ny - radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP5P9andP13Corners(&p, &p5, &p9, &p13, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p9 and p13
	for (i = 1; i < radius; i++) {
		for (j = ny - radius; j < ny - 1; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP9andP13Corners(&p, &p9, &p13, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
			else {
				out[zi] = 0;
			}
		}
	}

	//Only p1, p5 and p9
	for (i = radius; i < nx - radius; i++) {
		for (j = 1; j < radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1P5andP9Corners(&p, &p1, &p5, &p9, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	// Middle block, has p1, p5, p9 and p13 defined
#pragma omp parallel num_threads(num_fast_threads) private(zi, p, p1, p5, p9, p13)
	{
		int thread_id = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		int index = thread_id + 1;

		int total_rows = (nx - radius) - radius;
		int chunk_size = (total_rows + num_threads - 1) / num_threads;
		int start_row = radius + thread_id * chunk_size;
		int end_row = std::min(start_row + chunk_size, nx - radius);

		for (int i = start_row; i < end_row; i++) {
			for (int j = radius; j < ny - radius; j++) {
				zi = z(i, j, ny);
				out2[zi] = 0;

				if (checkP1P5P9andP13Corners(&p, &p1, &p5, &p9, &p13, in, i, j, ny, zi, threshold)) {
					results.data[index][results.dataSize[index]++] = zi;
					set3By3Neighbours(out, nx, ny, zi);
				}
			}
		}
	}

	//Only p1, p9 and p13
	for (i = radius; i < nx - radius; i++) {
		for (j = ny - radius; j < ny - 1; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1P9andP13Corners(&p, &p1, &p9, &p13, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p1 and p5
	for (i = nx - radius; i < nx - 1; i++) {
		for (j = 1; j < radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1andP5Corners(&p, &p1, &p5, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p1,p5 and p13
	for (i = nx - radius; i < nx - 1; i++) {
		for (j = radius; j < ny - radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1P5andP13Corners(&p, &p1, &p5, &p13, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p1 and p13
	for (i = nx - radius; i < nx - 1; i++) {
		for (j = ny - radius; j < ny - 1; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1andP13Corners(&p, &p1, &p13, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//auto time_Middle = std::chrono::high_resolution_clock::now();
	// //auto time_End = std::chrono::high_resolution_clock::now();
	//
	//auto time_FAST = time_End - time_Start;
	//
	//auto time_Candidates = time_End - time_Middle;
	//
	//printMessage("tiempo de duración FAST (PRÉ) (ms): " + std::to_string(time_FAST / std::chrono::milliseconds(1)));
	//printMessage("tiempo de duración Candidatos (ms): " + std::to_string(time_Candidates / std::chrono::milliseconds(1)));
}

void getFast4CandidatesEdges(float* in, float* out, float* out2, FastCandidates& results, int nx, int ny, float threshold, int num_fast_threads) {

	int num_sections = num_fast_threads + 2; // +2 for top and bottom

	results.data = (int**)malloc(sizeof(int*) * num_sections);

	results.dataSize = (int*)malloc(sizeof(int) * num_sections);

	results.num_sections = num_sections;

	size_t max_candidates_per_list = (nx * ny) / num_sections;

	for (int i = 0; i < num_sections; i++) {
		results.data[i] = (int*)malloc(sizeof(int) * max_candidates_per_list);
		results.dataSize[i] = 0;
	}

	auto time_Start = std::chrono::high_resolution_clock::now();

	int i, j, zi;
	float p, p1, p5, p9, p13;

	//Only p5 and p9
	for (i = 1; i < radius; i++) {
		for (j = 1; j < radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP5andP9Edges(&p, &p5, &p9, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p5,p9 and p13
	for (i = 1; i < radius; i++) {
		for (j = radius; j < ny - radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP5P9andP13Edges(&p, &p5, &p9, &p13, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p9 and p13
	for (i = 1; i < radius; i++) {
		for (j = ny - radius; j < ny - 1; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP9andP13Edges(&p, &p9, &p13, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
			else {
				out[zi] = 0;
			}
		}
	}

	//Only p1, p5 and p9
	for (i = radius; i < nx - radius; i++) {
		for (j = 1; j < radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1P5andP9Edges(&p, &p1, &p5, &p9, in, i, j, ny, zi, threshold)) {
				results.data[0][results.dataSize[0]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	// Middle block, has p1, p5, p9 and p13 defined
#pragma omp parallel num_threads(num_fast_threads) private(zi, p, p1, p5, p9, p13)
	{
		int thread_id = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		int index = thread_id + 1;

		int total_rows = (nx - radius) - radius;
		int chunk_size = (total_rows + num_threads - 1) / num_threads;
		int start_row = radius + thread_id * chunk_size;
		int end_row = std::min(start_row + chunk_size, nx - radius);

		for (int i = start_row; i < end_row; i++) {
			for (int j = radius; j < ny - radius; j++) {
				zi = z(i, j, ny);
				out2[zi] = 0;

				if (checkP1P5P9andP13Edges(&p, &p1, &p5, &p9, &p13, in, i, j, ny, zi, threshold)) {
					results.data[index][results.dataSize[index]++] = zi;
					set3By3Neighbours(out, nx, ny, zi);
				}
			}
		}
	}

	//Only p1, p9 and p13
	for (i = radius; i < nx - radius; i++) {
		for (j = ny - radius; j < ny - 1; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1P9andP13Edges(&p, &p1, &p9, &p13, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p1 and p5
	for (i = nx - radius; i < nx - 1; i++) {
		for (j = 1; j < radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1andP5Edges(&p, &p1, &p5, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p1,p5 and p13
	for (i = nx - radius; i < nx - 1; i++) {
		for (j = radius; j < ny - radius; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1P5andP13Edges(&p, &p1, &p5, &p13, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//Only p1 and p13
	for (i = nx - radius; i < nx - 1; i++) {
		for (j = ny - radius; j < ny - 1; j++) {
			zi = z(i, j, ny);
			out2[zi] = 0;
			if (checkP1andP13Edges(&p, &p1, &p13, in, i, j, ny, zi, threshold)) {
				results.data[num_fast_threads + 1][results.dataSize[num_fast_threads + 1]++] = zi;
				set3By3Neighbours(out, nx, ny, zi);

			}
		}
	}

	//auto time_Middle = std::chrono::high_resolution_clock::now();
	// //auto time_End = std::chrono::high_resolution_clock::now();
	//
	//auto time_FAST = time_End - time_Start;
	//
	//auto time_Candidates = time_End - time_Middle;
	//
	//printMessage("tiempo de duración FAST (PRÉ) (ms): " + std::to_string(time_FAST / std::chrono::milliseconds(1)));
	//printMessage("tiempo de duración Candidatos (ms): " + std::to_string(time_Candidates / std::chrono::milliseconds(1)));
}

void set3By3Neighbours(float* out, int nx, int ny, int zi) {
	int row = zi - ny - 1;
	out[row] = 1.0f;
	out[row + 1] = 1.0f;
	out[row + 2] = 1.0f;

	out[zi - 1] = 1.0f;
	out[zi] = 1.0f;
	out[zi + 1] = 1.0f;

	row = zi + ny - 1;
	out[row] = 1.0f;
	out[row + 1] = 1.0f;
	out[row + 2] = 1.0f;
}

bool checkP1andP5Corners(float* p, float* p1, float* p5, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);

	*p5 = getP5(in, i, j, ny);

	return diff(p, p1, threshold) && diff(p, p5, threshold);
}

bool checkP1andP5Edges(float* p, float* p1, float* p5, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);

	if (diff(p, p1, threshold)) {
		return true;
	}

	*p5 = getP5(in, i, j, ny);
	return diff(p, p5, threshold);
}

bool checkP5andP9Corners(float* p, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p5 = getP5(in, i, j, ny);

	*p9 = getP9(in, i, j, ny);

	return diff(p, p5, threshold) && diff(p, p9, threshold);
}

bool checkP5andP9Edges(float* p, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p5 = getP5(in, i, j, ny);

	if (diff(p, p5, threshold)) {
		return true;
	}

	*p9 = getP9(in, i, j, ny);
	return diff(p, p9, threshold);
}

bool checkP9andP13Corners(float* p, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p9 = getP9(in, i, j, ny);

	*p13 = getP13(in, i, j, ny);

	return diff(p, p9, threshold) && diff(p, p13, threshold);
}

bool checkP9andP13Edges(float* p, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p9 = getP9(in, i, j, ny);

	if (diff(p, p9, threshold)) {
		return true;
	}

	*p13 = getP13(in, i, j, ny);

	return diff(p, p13, threshold);
}

bool checkP1andP13Corners(float* p, float* p1, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);

	*p13 = getP13(in, i, j, ny);

	return diff(p, p1, threshold) && diff(p, p13, threshold);;
}

bool checkP1andP13Edges(float* p, float* p1, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);

	if (diff(p, p1, threshold)) {
		return true;
	}

	*p13 = getP13(in, i, j, ny);

	return diff(p, p13, threshold);;
}

bool checkP5P9andP13Corners(float* p, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p13 = getP13(in, i, j, ny);
	if (diff(p, p13, threshold))
	{
		*p5 = getP5(in, i, j, ny);
		if (diff(p, p5, threshold))
		{
			return true;
		}

		*p9 = getP9(in, i, j, ny);
		if (diff(p, p9, threshold))
		{
			return true;
		}
	}
	else {
		return checkP5andP9Corners(p, p5, p9, in, i, j, ny, zi, threshold);
	}
	return false;
}

bool checkP5P9andP13Edges(float* p, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p13 = getP13(in, i, j, ny);
	if (diff(p, p13, threshold))
	{
		return true;
	}

	*p5 = getP5(in, i, j, ny);
	if (diff(p, p5, threshold))
	{
		return true;
	}

	*p9 = getP9(in, i, j, ny);
	return diff(p, p9, threshold);
}

bool checkP1P5andP9Corners(float* p, float* p1, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);
	if (diff(p, p1, threshold))
	{
		*p5 = getP5(in, i, j, ny);
		if (diff(p, p5, threshold))
		{
			return true;
		}

		*p9 = getP9(in, i, j, ny);
		if (diff(p, p9, threshold))
		{
			return true;
		}
	}
	else {
		return checkP5andP9Corners(p, p5, p9, in, i, j, ny, zi, threshold);
	}
	return false;
}

bool checkP1P5andP9Edges(float* p, float* p1, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);
	if (diff(p, p1, threshold))
	{
		return true;
	}

	*p5 = getP5(in, i, j, ny);
	if (diff(p, p5, threshold))
	{
		return true;
	}

	*p9 = getP9(in, i, j, ny);

	return diff(p, p9, threshold);
}

bool checkP1P9andP13Corners(float* p, float* p1, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);
	if (diff(p, p1, threshold))
	{
		*p13 = getP13(in, i, j, ny);
		if (diff(p, p13, threshold))
		{
			return true;
		}

		*p9 = getP9(in, i, j, ny);
		if (diff(p, p9, threshold))
		{
			return true;
		}
	}
	else {
		return checkP9andP13Corners(p, p9, p13, in, i, j, ny, zi, threshold);
	}
	return false;
}

bool checkP1P9andP13Edges(float* p, float* p1, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p1 = getP1(in, i, j, ny);
	if (diff(p, p1, threshold))
	{
		return true;
	}

	*p13 = getP13(in, i, j, ny);
	if (diff(p, p13, threshold))
	{
		return true;
	}

	*p9 = getP9(in, i, j, ny);
	return diff(p, p9, threshold);
}

bool checkP1P5andP13Corners(float* p, float* p1, float* p5, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p5 = getP5(in, i, j, ny);
	if (diff(p, p5, threshold))
	{
		*p13 = getP13(in, i, j, ny);
		if (diff(p, p13, threshold))
		{
			return true;
		}

		*p1 = getP1(in, i, j, ny);
		if (diff(p, p1, threshold))
		{
			return true;
		}
	}
	else {
		return checkP1andP13Corners(p, p1, p13, in, i, j, ny, zi, threshold);
	}
	return false;
}

bool checkP1P5andP13Edges(float* p, float* p1, float* p5, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p5 = getP5(in, i, j, ny);
	if (diff(p, p5, threshold))
	{
		return true;
	}

	*p13 = getP13(in, i, j, ny);
	if (diff(p, p13, threshold))
	{
		return true;
	}

	*p1 = getP1(in, i, j, ny);
	return diff(p, p1, threshold);
}

bool checkP1P5P9andP13Corners(float* p, float* p1, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p13 = getP13(in, i, j, ny);
	if (diff(p, p13, threshold))
	{
		*p1 = getP1(in, i, j, ny);
		if (diff(p, p1, threshold))
		{
			return true;
		}

		*p5 = getP5(in, i, j, ny);
		if (diff(p, p5, threshold))
		{
			return true;
		}

		*p9 = getP9(in, i, j, ny);
		if (diff(p, p9, threshold))
		{
			return true;
		}
	}
	else
	{
		return checkP1P5andP9Corners(p, p1, p5, p9, in, i, j, ny, zi, threshold);
	}

	return false;
}

bool checkP1P5P9andP13Edges(float* p, float* p1, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold) {
	*p = in[zi];

	*p13 = getP13(in, i, j, ny);
	if (diff(p, p13, threshold))
	{
		return true;
	}

	*p1 = getP1(in, i, j, ny);
	if (diff(p, p1, threshold))
	{
		return true;
	}

	*p5 = getP5(in, i, j, ny);
	if (diff(p, p5, threshold))
	{
		return true;
	}

	*p9 = getP9(in, i, j, ny);
	return diff(p, p9, threshold);
}

bool diff(float* p, float* px, int threshold) {
	return *p > (*px + threshold) || *p < (*px - threshold);
}

float getP1(float* in, int i, int j, int ny) {
	return in[z(i - radius, j, ny)];
}

float getP5(float* in, int i, int j, int ny) {
	return in[z(i, j + radius, ny)];
}

float getP9(float* in, int i, int j, int ny) {
	return in[z(i + radius, j, ny)];
}

float getP13(float* in, int i, int j, int ny) {
	return in[z(i, j - radius, ny)];
}