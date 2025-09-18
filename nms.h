#pragma once
#include <iostream>
#include "image.h"

using namespace std;

std::vector<Corner> NonMaximumSuppression(float* out, int nx, int ny, float threshold, int radius);

void sortCornersSerial(std::vector<Corner>& candidates);