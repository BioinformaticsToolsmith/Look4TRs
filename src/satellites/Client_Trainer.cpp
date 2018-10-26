/**
 * Author: Hani Z. Girgis, PhD, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 * Based on code by Vincent Wells, The Bioinformatics Toolsmith Laboratory, The University of Tulsa
 *
 * This class uses the mean shift to train the HMM.
 */

// ToDo:
// Make this class work on one chromosome only
#ifndef CLIENT_TRAINER_CPP_
#define CLIENT_TRAINER_CPP_

#include "Client_Trainer.h"
#include "../nonltr/ChromosomeTR.h"
#include "../cluster/Progress.h"

 using namespace std;
 using namespace satellites;

 namespace satellites {



 	Client_Trainer::Client_Trainer(ChromosomeTR* chrom, HMM* hmmIn,
 		int minKIn, int maxKIn, int halfWIn, double baseIn,
 		int min_reg,
 		vector<double>& compListIn) :
 	IClient(hmmIn, compListIn, minKIn, maxKIn, halfWIn, baseIn) {
 		minK = minK;


 		randChrom = chrom;
 		scorer = makeScorer(randChrom);

 		trim();
 		train();
 	}

 	Client_Trainer::~Client_Trainer() {
 		Util::deleteInVector(chromSats);
 		// Util::deleteInVector(trimmedRegions);
 		// delete trimmedRegions;
 		delete chromSats;
 		delete scorer;

 	}

/**
 * This methods does the following:
 * 	(1) Puts tandem repeats in the random chromosome
 * 	(2) Trains the HMM
 */
 void Client_Trainer::train() {
	// // Train the HMM
 // 	cout << "\tAbout to train the HMM ..." << endl;

 	// hmm->train(scorer->getFlatScores(), randChrom->getSegment(), trimmedRegions);
 	hmm->train(scorer->getFlatScores(), randChrom->getSegment(), randChrom->getRegionList());

 	hmm->normalize();


 	chromSats = new vector<ILocation*>;
 	decode(randChrom, scorer, chromSats);

 	vector<char> * bestKList = scorer->getBestKList();
 	for (int j = 0; j < chromSats->size(); j++) {
 		auto sat = chromSats->at(j);

			// Extend the end
 		int bestK = bestKList->at(sat->getEnd());
 		if (bestK - 1 < 0) {
 			cerr << "Client_Scanner::get_hmm_sats - ";
 			cerr << "the extension amount cannot be negative";
 			cerr << endl;
 			throw std::exception();
 		}
 		sat->setEnd(sat->getEnd() + bestK - 1);
 	}
 	Util::merge(chromSats);

 	
 	// for(int i = 0; i < msRegionList->size(); i++){
 	// 	std::cout << msRegionList->at(i)->toString() << std::endl;
 	// 	std::cout << randChrom->getSubStr(msRegionList->at(i)->getStart(), msRegionList->at(i)->getEnd()) << std::endl;
 	// }
 }

/**
* This method will trim the trimmedRegions based on minK.
*/
void Client_Trainer::trim(){
	// trimmedRegions = new vector<ILocation *>(randChrom->getRegionList()->size());
	// for(int i = 0; i < randChrom->getRegionList()->size(); i++){
	// 	trimmedRegions->at(i) = new Location(*(randChrom->getRegionList()->at(i)));
	// }


	// for(auto it = trimmedRegions->begin(); it != trimmedRegions->end(); it++){
	// 	if((*it)->getEnd() - minK + 1 < (*it)->getStart()){
 // 			std::cerr << "Cannot trim! End - minK is smaller than start!" << std::endl;
 // 			throw std::exception();
	// 		continue;
	// 	}
	// 	(*it)->setEnd((*it)->getEnd() - minK + 1);
	// }
}

void Client_Trainer::trainPredictor(Predictor<int> * pred, int threshold){

	std::vector<std::string> * repeatList = randChrom->getRepeatList();
	for(int i = 0; i < repeatList->size(); i++){
		motif::FindMotif * f = new motif::FindMotif(repeatList->at(i), threshold, pred);
		delete f;
		if(pred->get_is_trained()){
			return;
		}

	}
	if(!pred->get_is_trained()){
		std::cerr << "Predictor is not trained! The ChromosomeTR has too few repeats!" << std::endl;
		throw std::exception();
	}
}


vector<ILocation*>* Client_Trainer::getChromSats(){
	return chromSats;
}



/*
	The returns the precision obtain from training. Precision is TP / GT.
*/
	double Client_Trainer::getSensitivity(){

		std::vector<Location *> * intersection = Util::locationIntersect(chromSats, randChrom->getRegionList());

		double TP = Util::sumTotalLength(intersection);
		double GT = Util::sumTotalLength(randChrom->getRegionList());

		double sens = 0;
		if(GT == 0){
			sens = 0;
		} else {
			sens = TP / GT;
		}

		Util::deleteInVector(intersection);
		delete intersection;
		return sens;
	}

/*
	The returns the sensitivity obtain from training. Sensitivity is TP / (TP + FP).
*/
	double Client_Trainer::getPrecision(){
		

		std::vector<Location *> * subtraction = Util::locationSubtract(chromSats, randChrom->getRegionList());

		double FP = Util::sumTotalLength(subtraction);


		std::vector<Location *> * intersection = Util::locationIntersect(chromSats, randChrom->getRegionList());

		double TP = Util::sumTotalLength(intersection);

		double prec = 0;

		if(TP + FP == 0){
			prec = 0;
		} else {
			prec = TP / (TP + FP);
		}

		Util::deleteInVector(intersection);
		Util::deleteInVector(subtraction);
		delete intersection;
		delete subtraction;
		return prec;
	}

	double Client_Trainer::getFMeasure(){
		double sens = getSensitivity();
		double prec = getPrecision();
		if(sens + prec == 0){
			return 0;
		}
		return 2 * sens * prec / (sens + prec);

	}


}

#endif
