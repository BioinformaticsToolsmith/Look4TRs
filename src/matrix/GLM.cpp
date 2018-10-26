/*
 * glm.cpp
 *
 * Created on: May 29, 2017
 * Author: Robert Geraghty, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 *
 */

#include "GLM.h"
#include "Matrix.h"

#include <math.h>
#include <iostream>
using namespace std;

namespace matrix {

void GLM::train(Matrix& features, Matrix& labels) {
	weights = features.transpose() * features;
	weights = weights.pseudoInverse() * features.transpose() * labels;
}

Matrix GLM::predict(Matrix& features) const {
	Matrix labels;
	labels = features * weights;
	double log;
	for (int i = 0; i < labels.getNumRow(); i++) {
		log = round(1 / (1 + exp(-(labels.get(i, 0)))));
		labels.set(i, 0, log);
	}
	return labels;
}

tuple<double, double, double> GLM::accuracy(Matrix& oLabels, Matrix& pLabels) const {
	int sum = 0;
	int negSum = 0;
	int negSame = 0;
	int posSum = 0;
	int posSame = 0;
	for (int i = 0; i < oLabels.getNumRow(); i++) {
		if (oLabels.get(i, 0) == -1) {
			negSum++;
			if (oLabels.get(i, 0) == pLabels.get(i, 0)) {
				sum++;
				negSame++;
			}
		} else {
			posSum++;
			if (oLabels.get(i, 0) == pLabels.get(i, 0)) {
				sum++;
				posSame++;
			}
		}
	}

	return make_tuple<double, double, double>(
			(((double) sum * 100) / (oLabels.getNumRow())),
			(((double) posSame * 100) / (posSum)),
			(((double) negSame * 100) / (negSum)));

}

}

