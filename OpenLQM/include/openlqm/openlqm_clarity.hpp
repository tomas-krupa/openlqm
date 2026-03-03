/*
* Copyright 2025 Noblis, Inc.
* https://fingerprint.nist.gov/openlqm
* 
* Licensed under the Apache License, Version 2.0 (the "License"); you may not
* use this file except in compliance with the License. You may obtain a copy of
* the License at https://www.apache.org/licenses/LICENSE-2.0.
* 
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
* License for the specific language governing permissions and limitations under
* the License.
* 
* OpenLQM was developed by Noblis, Inc. under contract to the National
* Institute of Standards and Technology in 2024-2025, as a port of "LQMetric,"
* which was developed by Noblis, Inc. under contract to the Federal Bureau of
* Investigation's Criminal Justice Information Services Division in 2012-2014.
*/

#pragma once

#include <opencv2/opencv.hpp>
#include <openlqm/openlqm.hpp>

#define LQM_ATTRIBUTES_PER_SAMPLE 21
#define LQM_NUMBER_OF_CLASSES 5

namespace OpenLQM {
	namespace Core {
		void GenerateNewClarityMapRandomForest(
			int mapWidth, int mapHeight,
			cv::Mat& clarityMap,
			const unsigned char* pGrayRangeMap,
			const unsigned char* pCurvatureMap,
			const unsigned char* pDirChangeMap,
			const double* pLowFreqMagMap,
			const double* pMaxMagMap,
			const unsigned char* pGrayMedianMap,
			const unsigned char* pGrayCountMap,
			const unsigned char* pValidNeighbors,
			const double* pNormMagMap,
			const unsigned char* pQualMap
		);

		void absDeviation(const cv::Mat& inputFeatureMap, cv::Mat& outputAggFeatureMap, const int filterSize);
		double calcAbsDeviation(const cv::Mat& sData, const int filterSize);
		void CropToLargestContiguousRegionOfRidgeFlow(cv::Mat& clarityMap);
		void CropToRoi(cv::Mat& clarityMap, const std::vector<OpenLQM::Coordinate>& roi);
		void FillContourHoles(cv::Mat& input);

		// This function is a partial port of OpenCV 2.4.12 CvRTrees::predict(). This was necessary because changes in the new RTrees
		//     implementation cause cv::ml::RTrees::predict() to give different outputs than LQMetric for some inputs.
		// It only supports continuous input samples where all variables in the input are ordered (non-categorical)
		float RandomForestPredictClass(const std::vector<int>& roots, const std::vector<cv::ml::DTrees::Node>& nodes, const std::vector<cv::ml::DTrees::Split>& splits, const int numTrees, const int numClasses, cv::AutoBuffer<int>& votes, const cv::Mat& sample);
	}
}
