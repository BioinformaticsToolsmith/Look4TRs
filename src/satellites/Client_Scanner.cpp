#include "Client_Scanner.h"

using namespace std;
using namespace satellites;

namespace satellites {

	Client_Scanner::Client_Scanner(const vector<ChromosomeOneDigit*>* chromListIn,
		HMM* hmm, vector<double>& compList, int minK, int maxK, int halfW,
		double base, double idnIn, int smoothingWindowIn, int mtfIn,
		Predictor<int> * predIn, int minRegIn, bool will_merge_in) :
	IClient(hmm, compList, minK, maxK, halfW, base) {

		chromList = chromListIn;
		hmm_sats = new vector<vector<ILocation*>*>();
		idn = idnIn;
		smoothingWindow = smoothingWindowIn;
		mtf = mtfIn;
		pred = predIn;
		minReg = minRegIn;
		will_merge = will_merge_in;
	}

	Client_Scanner::~Client_Scanner() {
		for (auto satList : *hmm_sats) {
			Util::deleteInVector(satList);
		}
		Util::deleteInVector(hmm_sats);
		hmm_sats->clear();
		delete hmm_sats;
	}

/**
 * This method uses the HMM to detect satellites
 */
	void Client_Scanner::get_hmm_sats(
		vector<tuple<ILocation*, ChromosomeOneDigit *, string, string, double> >& output) {
		for (int i = 0; i < chromList->size(); i++) {

		// Score the chromosome
			auto scorer = makeScorer(chromList->at(i));
			vector<char> * bestKList = scorer->getBestKList();

		// Find STR using the HMM
			vector<ILocation*> * chromSats = new vector<ILocation*>();
			decode(chromList->at(i), scorer, chromSats);

			hmm_sats->push_back(chromSats);

		// Collect the results from the cores
			vector<tuple<ILocation*, ChromosomeOneDigit *, string, string, double>> coreCollect(
				chromSats->size());

		// Find the repeated motif in this STR
			for (int j = 0; j < chromSats->size(); j++) {
				auto sat = chromSats->at(j);

			// Extend the end
				int bestK = bestKList->at(sat->getEnd());
				if (bestK < 0) {
					cerr << "Client_Scanner::get_hmm_sats - ";
					cerr << "the extension amount cannot be negative";
					cerr << endl;
					throw std::exception();
				} else if (bestK > 0 ) {
					sat->setEnd(sat->getEnd() + bestK - 1);
				}
				
			// Convert digits to nucleotides
				string candidate = oneDigitToNucleotide(chromList->at(i)->getBase(),
					sat->getStart(), sat->getLength());

			// Without the following condition the motif discovery module
			// will fail when the region is smaller than twice the smoothing
			// window
				if (sat->getLength() > 2 * smoothingWindow
					&& sat->getLength() > (minReg / 2.0)) {
					if (mtf) {
						FindMotif * findMotif;
						string candidateSample;

						candidateSample = candidate.substr(0, 5000);

					// We search for a micro region
						findMotif = new FindMotif(candidateSample, idn, pred);
						

						if (findMotif->getIsFound()) {
							coreCollect.at(j) = make_tuple(sat, chromList->at(i),
								candidate, findMotif->getFoundMotif(),
								findMotif->getIdentityScore());
						} else {
							coreCollect.at(j) = make_tuple(sat, chromList->at(i),
								candidate, string("-"), 0.0);
						}
						delete findMotif;
					} else {
						coreCollect.at(j) = make_tuple(sat, chromList->at(i),
							candidate, string("-"), 0.0);
					}

				} else {
					coreCollect.at(j) = make_tuple(sat, chromList->at(i), candidate,
						string("-"), 0.0);
				}



			}

		// Collect and filter the results. A region is considered
		// only when the identity score between the exact repeat (synthetic)
		// and the candidate region is above the threshold
			for (int j = 0; j < coreCollect.size(); j++) {
				if (!mtf || std::get<4>(coreCollect.at(j)) >= idn) {
					output.push_back(coreCollect.at(j));
				}
			}

		// Free up memory
			delete scorer;
		}
	}

/**
 * Parameter chrom: is a single digit chromosome that we search for STR
 * Parameter scorer: is the corresponding scores for the chromosome
 * Parameter chromSats: is a list of the found satellites
 */
	void Client_Scanner::decode(ChromosomeOneDigit * chrom, ScorerSat * scorer,
		vector<ILocation*>* chromSats) {

		auto segmentList = chrom->getSegment();
		for (int i = 0; i < segmentList->size(); i++) {
			std:;vector<ILocation *> mergingChromSats;

			auto segment = segmentList->at(i);
			int segStart = segment->at(0);
			int segEnd = segment->at(1);

			hmm->decode(segStart, segEnd, scorer->getFlatScores(), mergingChromSats);

			if(will_merge){
				Util::merge(chromSats, minReg);
			}

			chromSats->insert(std::end(*chromSats), std::begin(mergingChromSats), std::end(mergingChromSats));
		}
	}


	string Client_Scanner::oneDigitToNucleotide(const string * seq, int index,
		int len) {
		string seed("");
		for (int i = index; i < index + len; i++) {
			char n = seq->at(i);

			if (n == 0) {
				seed.append(1, 'A');
			} else if (n == 1) {
				seed.append(1, 'C');
			} else if (n == 2) {
				seed.append(1, 'G');
			} else if (n == 3) {
				seed.append(1, 'T');
			} else if (n == 78) {
				seed.append(1, 'C');
			}
			else{
				cerr << "Client_Scanner::oneDigitToNucleotide - ";
				cerr << "undefined base code: " << (int) n << endl;
				throw std::exception();
			}
		}
		return seed;
	}

}
