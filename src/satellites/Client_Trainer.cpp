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



  Client_Trainer::Client_Trainer(ChromosomeTR* trainChrom_in, ChromosomeTR* testChrom_in, HMM* hmmIn,
   int minKIn, int maxKIn, int halfWIn, double baseIn,
   int min_reg,
   vector<double>& compListIn) :
  IClient(hmmIn, compListIn, minKIn, maxKIn, halfWIn, baseIn) {
    minK = minK;
    
    
    trainChrom = trainChrom_in;
    testChrom = testChrom_in;
    train_scorer = makeScorer(trainChrom);
    test_scorer = makeScorer(testChrom);
    
    train();
  }
  
  Client_Trainer::~Client_Trainer() {
    Util::deleteInVector(chromSats);
    delete chromSats;
    delete train_scorer;
    delete test_scorer;
    
  }
  
  /**
   * This methods does the following:
   * 	(1) Puts tandem repeats in the random chromosome
   * 	(2) Trains the HMM
   */
  void Client_Trainer::train() {
    // Train the HMM
    hmm->train(train_scorer->getFlatScores(), trainChrom->getSegment(), trainChrom->getRegionList());
    
    hmm->normalize();
    
    
    chromSats = new vector<ILocation*>;
    decode(testChrom, test_scorer, chromSats);
    
    
    
    vector<char> * bestKList = test_scorer->getBestKList();
    for (int j = 0; j < chromSats->size(); j++) {
      auto sat = chromSats->at(j);
      
      // Extend the end
      int bestK = bestKList->at(sat->getEnd());
      if (bestK < 0) {
       std::cout << sat->toString() << std::endl;
       for(int k = sat->getStart(); k < sat->getEnd(); k++){
         std::cout << bestKList->at(k) << " ";
       }
       std::cout << std::endl;

       std::string sstring = testChrom->getSubStr(sat->getStart(), sat->getEnd());
       for(auto it = sstring.begin(); it != sstring.end(); it++){
         std::cout << (int) *it << " ";
       }
       std::cout << std::endl;
       cerr << "Client_Trainer::get_hmm_sats - ";
       cerr << "the extension amount cannot be negative: ";
       cerr << bestK;
       cerr << endl;
       throw std::exception();
     } else if (bestK > 0){
       sat->setEnd(sat->getEnd() + bestK - 1);
     }
   }
   Util::merge(chromSats);

 }


 void Client_Trainer::trainPredictor(Predictor<int> * pred, int threshold){

  std::vector<std::string> * repeatList = testChrom->getRepeatList();
  for(int i = 0; i < repeatList->size(); i++){
    motif::FindMotif * f = new motif::FindMotif(repeatList->at(i), threshold, pred);
    delete f;
    if(pred->get_is_trained()){
     return;
   }

 }
 while(!pred->get_is_trained()){
  std::string artificalRepeat = testChrom->getRandTR();
  std:: cout << artificalRepeat << std::endl;
  motif::FindMotif * f = new motif::FindMotif(artificalRepeat, threshold, pred);
  delete f;
}
}



vector<ILocation*>* Client_Trainer::getChromSats(){
  return chromSats;
}



  /*
    The returns the precision obtain from training. Precision is TP / GT.
  */
double Client_Trainer::getSensitivity(){

  std::vector<Location *> * intersection = Util::locationIntersect(chromSats, testChrom->getRegionList());

  double TP = Util::sumTotalLength(intersection);
  double GT = Util::sumTotalLength(testChrom->getRegionList());

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


  std::vector<Location *> * subtraction = Util::locationSubtract(chromSats, testChrom->getRegionList());

  double FP = Util::sumTotalLength(subtraction);


  std::vector<Location *> * intersection = Util::locationIntersect(chromSats, testChrom->getRegionList());

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
