#pragma once

#include <opencv2/opencv.hpp>
#include <vector>

namespace OpenLQM {
	namespace Impl {
		extern std::vector<int> roots;
		extern std::vector<cv::ml::DTrees::Node> nodes;
		extern std::vector<cv::ml::DTrees::Split> splits;
        extern int numTrees;
        extern int numClasses;
	}
}