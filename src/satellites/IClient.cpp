/**
 * Author: Vincent Wells, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Refactored by Hani Z. Girgis, PhD, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 */

#include "IClient.h"

using namespace std;
using namespace satellites;

namespace satellites {

IClient::IClient(HMM * hmmIn, vector<double>& compListIn, int minKIn,
		int maxKIn, int halfWIn, double baseIn) :
		compList(compListIn) {

	hmm = hmmIn;
	compList = compListIn;
	minK = minKIn;
	maxK = maxKIn;
	halfW = halfWIn;
	base = baseIn;
	tableList = new vector<ScorerAdjusted*>();

	// This list of table is used by ScorerSat
	for (int k = minK; k <= maxK; k++) {
		tableList->push_back(new ScorerAdjusted(k, halfW, compList));
	}
}

IClient::~IClient() {
	Util::deleteInVector(tableList);
	delete tableList;
}

ScorerSat * IClient::makeScorer(ChromosomeOneDigit * chrom) {
	auto scorer = new ScorerSat(*chrom, minK, maxK, halfW, compList, tableList,
			base);
	// This is where we print the scores
	// scorer->printScores(chrom->getHeader(), false);
	return scorer;
}

HMM * IClient::getHMM() {
	return hmm;
}


/**
 * Parameter chrom: is a single digit chromosome that we search for STR
 * Parameter scorer: is the corresponding scores for the chromosome
 * Parameter chromSats: is a list of the found satellites
 */
void IClient::decode(ChromosomeOneDigit * chrom, ScorerSat * scorer,
		vector<ILocation*>* chromSats) {

	auto segmentList = chrom->getSegment();
	for (int i = 0; i < segmentList->size(); i++) {
		auto segment = segmentList->at(i);
		int segStart = segment->at(0);
		int segEnd = segment->at(1);
		hmm->decode(segStart, segEnd, scorer->getFlatScores(), *chromSats);
	}
}


}
