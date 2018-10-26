
#include "ChromosomeTR.h"
#include <random>
#include "../motif/FindMotif.h"
#include "../mutate/HandleSeq.h"
#include <random>
#include <math.h>


ChromosomeTR::ChromosomeTR(int nIn, ChromosomeOneDigit* oChromIn,
	char unreadIn, int length, int minMotifIn, int maxMotifIn, int init_regIn): ChromosomeRandom(nIn, oChromIn, unreadIn, length){
	regionList = new std::vector<ILocation *>();
	repeatList = new std::vector<std::string>();
	init_reg = init_regIn;
	minMotif = minMotifIn;
	maxMotif = maxMotifIn;
	shuffle();

}


ChromosomeTR::~ChromosomeTR(){
	Util::deleteInVector(regionList);
	delete repeatList;
	delete regionList;
}

void ChromosomeTR::removeLocations(const std::vector<Location *> * deleteList)
{
	if(deleteList->back()->getStart() > rBase->size() || deleteList->back()->getEnd() > rBase->size()){
		std::cerr << "The delete list has an out of bound region" << deleteList->back()->toString() << "!" << std::endl;
		throw std::exception();
	}

	auto region = deleteList->begin();
	std::default_random_engine generator(0);
	std::uniform_int_distribution<int> exactRepeatDistr(0, 3);


	for(int index = 0; index < rBase->size(); index++){

		if(region != deleteList->end()){

			if(index < (*region)->getStart()){
				index = (*region)->getStart() - 1;
			} else if (index > (*region)->getEnd()){
				index--;
				region++;
			} else if(index >= (*region)->getStart()){
				rBase->at(index) = (char) exactRepeatDistr(generator);
			} 

		} else {
			break;
		}
	}
}


void ChromosomeTR::replaceSeq(std::string& subSeq, Location& l){
	int originalSize = rBase->size();
	rBase->replace(l.getStart(), l.getEnd() - l.getStart(), subSeq);
	if(rBase->size() < originalSize){
		std::cerr << "The base got smaller!" << std::endl;
		throw std::exception();
	} else if (rBase->size() > originalSize){
		std::cerr << "The base got bigger!" << std::endl;
		throw std::exception();
	}
}

void ChromosomeTR::replaceSeq(std::string& subSeq, int start, int end){
	Location l(start, end);
	replaceSeq(subSeq, l);
}

std::string ChromosomeTR::getSubStr(int start, int end){
	return rBase->substr(start, end - start);
}

void ChromosomeTR::shuffle(){

	std::default_random_engine generator (0);


 	std::uniform_int_distribution<int> exactRepeatDistr(init_reg, 20 * init_reg); // init_reg, 20 * init_reg
 	std::uniform_int_distribution<int> wordSizeDistr(minMotif, maxMotif);

 	int mutationRate = 0;
 	int globalMutationSum = 0;
 	int numRepeats  = 0;

 	for(int i = 0; i < getSegment()->size(); i++){

 		int start = getSegment()->at(i)->at(0);
 		int end = getSegment()->at(i)->at(1) + 1; // The end is inclusive
 		std::string a = getBase()->substr(start, end - start);
 		int mutationSum = 0;
 		// Choose between minK and maxK every time
 		int wordSize;
 		std::uniform_int_distribution<int> findPoint(0, 20);



 		std::reverse(a.begin(), a.end());
 		std::string b ("");
 		for(int seqLen = 0; a.size() > 0;){

 			if(a.size() > 20 * init_reg && double(mutationSum) / seqLen < 0.05 && findPoint(generator) == 0){
 				numRepeats++;

 				int exactLength =  exactRepeatDistr(generator);
 				while(exactLength > a.size()){
 					exactLength =  exactRepeatDistr(generator);
 				}

 				wordSize = wordSizeDistr(generator);
 				if (exactLength / wordSize < 2){
 					std::uniform_int_distribution<int> newExactLengthDist(3, 20);
 					exactLength = wordSize * newExactLengthDist(generator);
 				}


 				int TRSize = (exactLength / wordSize) * wordSize;
 				std::string motif = a.substr(a.size() - wordSize, wordSize);
 				std::reverse(motif.begin(), motif.end());

 				std::string newSubSeq = motif::FindMotif::makeExact(motif, TRSize);

 				HandleSeq handleSeq(1, true);

 				auto results = handleSeq.mutate(newSubSeq, mutationRate);
 				newSubSeq = results.second;
 				newSubSeq = newSubSeq.substr(0, min(newSubSeq.size(), a.size()));

 				mutationRate++;
 				mutationRate %= 26;
 			

 				TRSize = newSubSeq.size();
 				b += newSubSeq;
 				repeatList->push_back(newSubSeq);

 				for(int j = 0; j < TRSize; j++){
 					a.pop_back();
 				}
 				regionList->push_back(new Location(b.size() - TRSize + start, b.size() + start));
 				seqLen += TRSize;
 				mutationSum += TRSize;

 			} else {
 				char newChar = a.back();
 				a.pop_back();
 				b += newChar;
 				seqLen++;
 			}

 		}
 		// std::cout << start << " " << end << std::endl;
 		replaceSeq(b, start, start + b.size());
 		globalMutationSum += mutationSum;
 	}

 	setK( ceil( log( globalMutationSum / (double) numRepeats) / log(4) ) - 1 ); 
 }


 void ChromosomeTR::setK(int a){
 	if (a <= 0){
 		std::cerr << "Invalid choice of K (" << a << ")!" << std::endl;
 		throw std::exception();
 	}
 	K = a;
 }

 int ChromosomeTR::getK(){
 	return K;
 }


 void ChromosomeTR::printRBase(){
 	std::cout << *(rBase) << std::endl;
 }

 int ChromosomeTR::size(){
 	return rBase->size();
 }

 std::vector<ILocation *> * ChromosomeTR::getRegionList(){
 	return regionList;
 }

 std::vector<std::string> * ChromosomeTR::getRepeatList(){
 	return repeatList;
 }