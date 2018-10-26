/*
 * ScorerSat.cpp
 *
 *  Created on: Jan 6, 2017
 *      Author: Hani Zakaria Girgis, PhD
 */
#include "ScorerSat.h"

#include <iostream>
#include <math.h>

namespace satellites {

ScorerSat::ScorerSat(ChromosomeOneDigit& chromIn, int minKIn, int maxKIn,
		int halfWIn, vector<double>& compListIn,
		vector<ScorerAdjusted *> * scorerListIn, double b) :
		chrom(chromIn), compList(compListIn) {

	minK = minKIn;
	maxK = maxKIn;
	halfW = halfWIn;
	base = b;
	canDeleteZeroed = false;
	bestKList = new vector<char>(chrom.getBase()->size(), 0);
	adjustedList = new vector<int>(chrom.getBase()->size(), 0);

	scorerList = scorerListIn;

	const vector<vector<int> *> * segmentList = chrom.getSegment();

	for (int i = 0; i < segmentList->size(); i++) {
		int start = segmentList->at(i)->at(0);
		int end = segmentList->at(i)->at(1);
		if (halfW <= end - maxK + 1 - start) {
			processSegment(start, end);
		} else {
			cout << "\tSkipped ..." << endl;
		}
	}

	make_flattened();
}

ScorerSat::~ScorerSat() {
	bestKList->clear();
	delete bestKList;

	adjustedList->clear();
	delete adjustedList;

	adjustedList_flattened->clear();
	delete adjustedList_flattened;

	if (canDeleteZeroed) {
		adjustedList_zeroed->clear();
		delete adjustedList_zeroed;
	}
}

const bool ScorerSat::getCanDeleteZeroed() {
	return canDeleteZeroed;
}

void ScorerSat::processSegment(int start, int end) {
	int score = -1;
	int bestK = -1;

	for (int i = minK; i <= maxK; i++) {
		// Get the score of the nucleotide at the start
		int kMerScore = scorerList->at(i - minK)->processSegment(
				chrom.getBase(), start, end);
		if (kMerScore >= score) {
			score = kMerScore;
			bestK = i;
		}
	}
	bestKList->at(start) = bestK;
	adjustedList->at(start) = score;

	// The last scorer has the shortest list of scores because it has the longest k
	int firstEnd = end - maxK + 1;
	for (int h = start + 1; h <= firstEnd; h++) {
		int score = -1;
		int bestK = -1;

		for (int i = minK; i <= maxK; i++) {
			int kMerScore = scorerList->at(i - minK)->moveOneNucleotide();

			if (kMerScore >= score) {
				score = kMerScore;
				bestK = i;
			}
		}
		bestKList->at(h) = bestK;
		adjustedList->at(h) = score;
	}

	// Handle the last maxK nucleotides
	for (int h = firstEnd + 1, j = scorerList->size() - 2;
			j >= 0 && (h <= end - minK + 1); j--, h++) {
		int score = -1;
		int bestK = -1;
		for (int i = 0; i <= j; i++) {
			int kMerScore = scorerList->at(i)->moveOneNucleotide();

			if (kMerScore >= score) {
				score = kMerScore;
				bestK = i + minK;
			}
		}
		bestKList->at(h) = bestK;
		adjustedList->at(h) = score;
	}

	// Handle the last minK nucleotides
	int lastScore = adjustedList->at(end - minK + 1);
	int lastBestK = bestKList->at(end - minK + 1);
	for (int h = end - minK + 2; h <= end; h++) {
		if (lastBestK < 0) {
			throw "The extension amount cannot be negative.\n";
		}
		lastBestK = lastBestK - 1;
		bestKList->at(h) = lastBestK;
		adjustedList->at(h) = lastScore;
	}

	// These scorers are ready to be used for another segment
	for (int i = minK; i <= maxK; i++) {
		scorerList->at(i - minK)->clear();
	}
}

void ScorerSat::make_flattened() {
	adjustedList_flattened = new vector<int>(*adjustedList);
	for (auto& score : (*adjustedList_flattened)) {
		if (score != 0) {
			score = round(log(score) / log(base));
		} else {
			score = 0;
		}
	}
}

void ScorerSat::make_zeroed(double mean) {
	adjustedList_zeroed = new vector<int>(*adjustedList);
	for (int& score : *(adjustedList_zeroed)) {
		if (score < mean) {
			score = 0;
		}
	}
	canDeleteZeroed = true;
}

vector<int>* ScorerSat::getScores() const {
	return adjustedList;
}

vector<int>* ScorerSat::getFlatScores() const {
	return adjustedList_flattened;
}

vector<int>* ScorerSat::getZeroScores() const {
	return adjustedList_zeroed;
}

ChromosomeOneDigit* ScorerSat::get_chrom() const {
	return &chrom;
}

vector<double>* ScorerSat::get_compList() const {
	return &compList;
}

vector<char> * ScorerSat::getBestKList() const {
	return bestKList;
}

double ScorerSat::avgScore() {
	double n = (double) accumulate(adjustedList->begin(), adjustedList->end(),
			0) / adjustedList->size();
	return n;
}

double ScorerSat::standardDeviation(double mean) {
	double total = 0;
	for (int score : (*adjustedList)) {
		double temp = score - mean;
		total += temp * temp;
	}
	return sqrt(total / adjustedList->size());
}

int ScorerSat::getMax() {
	if (adjustedList_flattened->size() > 0) {
		int s = adjustedList_flattened->at(0);
		for (int i = 1; i < adjustedList_flattened->size(); i++) {
			if (adjustedList_flattened->at(i) > s) {
				s = adjustedList_flattened->at(i);
			}
		}
		return s;
	} else {
		throw InvalidStateException("The score list is empty\n");
	}
}

void ScorerSat::printScores(string outputFile, bool canAppend) {
	ofstream outScores;
	if (canAppend) {
		outScores.open(outputFile.c_str(), ios::out | ios::app);
	} else {
		outScores.open(outputFile.c_str(), ios::out);
	}

	int step = 50;
	outScores << chrom.getHeader() << endl;
	int len = adjustedList->size();
	for (int i = 0; i < len; i = i + step) {
		int e = (i + step - 1 > len - 1) ? len - 1 : i + step - 1;
		for (int k = i; k <= e; k++) {
			outScores << adjustedList->at(k) << " ";
		}
		outScores << endl;
	}

	outScores << endl;
	outScores.close();
}

} /* namespace nonltr */
