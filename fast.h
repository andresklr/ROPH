#pragma once
#include "image.h"
#include "tools.h"

void getFast4CandidatesCorners(float* in, float* out, float* out2, FastCandidates& results, int nx, int ny, float threshold, int num_fast_threads);

void set3By3Neighbours(float* out, int nx, int ny, int zi);

void setCandidates(float* out, int nx, int ny, FastCandidates& results, bool padding);

void getFast4CandidatesEdges(float* in, float* out, float* out2, FastCandidates& results, int nx, int ny, float threshold, int num_fast_threads);

bool checkP1andP5Corners(float* p, float* p1, float* p5, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP1andP5Edges(float* p, float* p1, float* p5, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP5andP9Corners(float* p, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP5andP9Edges(float* p, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP9andP13Corners(float* p, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP9andP13Edges(float* p, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP1andP13Corners(float* p, float* p1, float* p13, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP1andP13Edges(float* p, float* p1, float* p13, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP5P9andP13Corners(float* p, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP5P9andP13Edges(float* p, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP1P5andP9Corners(float* p, float* p1, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP1P5andP9Edges(float* p, float* p1, float* p5, float* p9, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP1P9andP13Corners(float* p, float* p1, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP1P9andP13Edges(float* p, float* p1, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP1P5andP13Corners(float* p, float* p1, float* p5, float* p13, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP1P5andP13Edges(float* p, float* p1, float* p5, float* p13, float* in, int i, int j, int ny, int zi, int threshold);

bool checkP1P5P9andP13Corners(float* p, float* p1, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);
bool checkP1P5P9andP13Edges(float* p, float* p1, float* p5, float* p9, float* p13, float* in, int i, int j, int ny, int zi, int threshold);

bool diff(float* p, float* px, int threshold);

float getP1(float* in, int i, int j, int ny);

float getP5(float* in, int i, int j, int ny);

float getP9(float* in, int i, int j, int ny);

float getP13(float* in, int i, int j, int ny);