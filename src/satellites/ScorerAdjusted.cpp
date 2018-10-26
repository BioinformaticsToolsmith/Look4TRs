/*
 * ScorerAdjusted.cpp
 *
 *  Created on: Jan 8, 2017
 *      Author: Hani Zakaria Girgis, PhD
 *      --- half window --- center --- half window ---
 */

#include "ScorerAdjusted.h"

namespace nonltr {

ScorerAdjusted::ScorerAdjusted(int kIn, int halfWIn, vector<double>& compListIn) :
		compList(compListIn) {
	k = kIn;
	halfW = halfWIn;
	table = new KmerHashTable<unsigned long, int>(k, initValue);
	hasher = table->getHashMaker();
	hashList = new vector<unsigned long>();
	expectedScoreList = new vector<double>();
	center = -1;
	lastIndexInHashList = -1;
}

ScorerAdjusted::~ScorerAdjusted() {
	delete table;

	hashList->clear();
	delete hashList;

	expectedScoreList->clear();
	delete expectedScoreList;
}

void ScorerAdjusted::fillExpectedScoreList(int start, int end) {
	double expectedInWindow = log2(2 * halfW + 1);

	// Add logs of the first k
	for (int i = start; i < start + k; i++) {
		expectedInWindow += compList.at(seq->at(i));
	}
	expectedScoreList->push_back(expectedInWindow);

	// Subtract the first one and add the last one
	for (int i = start + k; i <= end; i++) {
		expectedInWindow -= compList.at(seq->at(i - k));
		expectedInWindow += compList.at(seq->at(i));
		expectedScoreList->push_back(expectedInWindow);
	}
}

int ScorerAdjusted::processSegment(const string* seqIn, int start, int end) {
	seq = seqIn;

	// Hash the whole segment at once
	hasher->hash(seq->c_str(), start, end - k + 1, hashList);

	lastIndexInHashList = hashList->size() - 1;

	// Calculate the expected count of each word in the window
	fillExpectedScoreList(start, end);

	// Make sure that the sizes of the observed and the expected scores are the same
	if (hashList->size() != expectedScoreList->size()) {
		string msg("Error: size of hash list != size of expected score list\n");
		throw InvalidStateException(msg);
	}

	//if (halfW <= end - k + 1 - start) {
	if (halfW < hashList->size()) {
		// Fill the first nucleotide and the right half window
		for (int i = 0; i < halfW + 1; i++) {
			table->increment(hashList->at(i));
		}

		// At this point the score of the first nucleotide is available
		center = 0;

		// Return the adjusted center of the center
		return getAdjustedScoreOfCenter();

	} else {
		string msg("Error: This segment is too short.\n");
		throw InvalidInputException(msg);
	}
}

int ScorerAdjusted::moveOneNucleotide() {
	// ToDO: Write a precondition checking the center and the lastIndex
	center++;
	if (center <= lastIndexInHashList) {
		int wStart = center - halfW;
		int wEnd = center + halfW;

		// Delete the first nucleotide in the window if applicable
		if (wStart - 1 >= 0) {
			table->decrement(hashList->at(wStart - 1));
		}

		// Add the last nucleotide if you can
		if (wEnd <= lastIndexInHashList) {
			table->increment(hashList->at(wEnd));
		}

		// Return the adjusted center of the center
		return getAdjustedScoreOfCenter();
	} else {
		string msg("Error: Center is beyond the current segment.\n");
		throw InvalidStateException(msg);
	}
}

int ScorerAdjusted::getScoreOfCenter() {
	return table->valueOf(hashList->at(center));
}

int ScorerAdjusted::getAdjustedScoreOfCenter() {
	double expected = round(pow(2, expectedScoreList->at(center)));
	int adjusted = table->valueOf(hashList->at(center)) - expected - 1;
	return (adjusted < 0) ? 0 : (adjusted * k);
}

void ScorerAdjusted::clear() {
	// Clear the contents within
	// --- half window --- center --- half window ---
	int wStart = center - halfW;
	int wEnd = center + halfW;

	if (wStart < 0) {
		wStart = 0;
	}

	if (wEnd > lastIndexInHashList) {
		wEnd = lastIndexInHashList;
	}

	for (int i = wStart; i <= wEnd; i++) {
		table->insert(hashList->at(i), this->initValue);
	}

	// Clear the hash list
	hashList->clear();
	expectedScoreList->clear();

	// Reset variables
	center = -1;
	lastIndexInHashList = -1;

}

} /* namespace nonltr */
