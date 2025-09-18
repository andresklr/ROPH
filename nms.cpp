#include "nms.h"

std::vector<Corner> NonMaximumSuppression(float* out, int nx, int ny, float threshold, int windowSize) {
    cv::Mat responseMap(nx, ny, CV_32F, out);

    cv::Mat localMax;

    cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(windowSize, windowSize));

    cv::dilate(responseMap, localMax, kernel);


    std::vector<Corner> corners;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            float response = responseMap.at<float>(i, j);
            if (response > threshold && response == localMax.at<float>(i, j)) {
                Corner c;
                c.i = i;
                c.j = j;
                c.r = response;
                c.skipped = false;
                corners.push_back(c);
            }
            else {
                out[z(i, j, ny)] = 0;
            }
        }
    }

    return corners;
}

/*std::vector<Corner> NonMaximumSuppressionParallel(float* out, int nx, int ny, float threshold) {

}*/

void sortCornersSerial(std::vector<Corner>& candidates) {
	std::sort(candidates.begin(), candidates.end(), [](const Corner& a, const Corner& b) {
		return a.r > b.r; // Maior 'r' primeiro (ordem decrescente)
		});
}

void sortCornersParallel(std::vector<Corner>& candidates) {
#ifdef HAS_PARALLEL_SORT
	std::sort(std::execution::par, candidates.begin(), candidates.end(), [](const Corner& a, const Corner& b) {
		return a.r > b.r; // Maior 'r' primeiro (ordem decrescente)
		});
#else
	std::cerr << "Aviso: C++17 parallel algorithms not supported by this compiler or environment. Falling back to serial sort." << std::endl;
	// Se o paralelismo não for suportado, ainda podemos chamar a versão serial como fallback
	sortCornersSerial(candidates);
#endif
}
