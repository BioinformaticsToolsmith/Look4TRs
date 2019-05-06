
#include "ChromosomeSelfTR.h"
#include <random>
#include "../motif/FindMotif.h"
#include "../mutate/HandleSeq.h"
#include <random>
#include <math.h>

ChromosomeSelfTR::ChromosomeSelfTR(int nIn, ChromosomeOneDigit *oChromIn,
	char unreadIn, int length, int minMotifIn, int maxMotifIn, int init_regIn, int seed_in) : ChromosomeTR(nIn, oChromIn, unreadIn, length, seed_in)
{
	regionList = new std::vector<ILocation *>();
	repeatList = new std::vector<std::string>();
	motifList = new std::vector<std::string>();
	init_reg = init_regIn;
	minMotif = minMotifIn;
	maxMotif = maxMotifIn;
	baseCountOChrom = oChromIn->getBaseCount();
	int sum = 0;
	for(int i = 0; i < baseCountOChrom->size(); i++){
		sum += baseCountOChrom->at(i);
	}
	for(int i = 0; i < baseCountOChrom->size(); i++){
		baseCountOChrom->at(i) = (baseCountOChrom->at(i) * 100) / sum;
	}
	shuffle();
}

ChromosomeSelfTR::~ChromosomeSelfTR()
{
	delete motifList;
}


std::string ChromosomeSelfTR::getRandTR()
{

	std::default_random_engine generator(seed);
	seed++;
	std::uniform_int_distribution<int> segmentGen(0, getSegment()->size() - 1);
	int segmentChoice = segmentGen(generator);

	int start = getSegment()->at(segmentChoice)->at(0);
	int end = getSegment()->at(segmentChoice)->at(1) + 1;

	std::uniform_int_distribution<int> motifIndex(0, end - start - 1 - maxMotif);
	std::uniform_int_distribution<int> exactRepeatDistr(init_reg, 20 * init_reg);
	std::uniform_int_distribution<int> wordSizeDistr(minMotif, maxMotif);

	int wordSize = wordSizeDistr(generator);
	int motifStart = motifIndex(generator);
	std::string motif = getBase()->substr(start, end - start).substr(motifStart, wordSize);
	int exactLength = exactRepeatDistr(generator);

	int TRSize = (exactLength / wordSize) * wordSize;

	std::string newSubSeq = motif::FindMotif::makeExact(motif, TRSize);

	HandleSeq handleSeq(1, true);

	auto results = handleSeq.mutate(newSubSeq, mutationRate, baseCountOChrom);
	newSubSeq = results.second;

	mutationRate++;
	mutationRate = max(mutationRate % 26, 0);
	return newSubSeq;
}

void ChromosomeSelfTR::shuffle()
{

	std::default_random_engine generator(seed);
	seed++;

	std::uniform_int_distribution<int> exactRepeatDistr(init_reg, 20 * init_reg); // init_reg, 20 * init_reg
	std::uniform_int_distribution<int> wordSizeDistr(minMotif, maxMotif);

	int globalMutationSum = 0;
	int numRepeats = 0;

	for (int i = 0; i < getSegment()->size(); i++)
	{

		int start = getSegment()->at(i)->at(0);
		int end = getSegment()->at(i)->at(1) + 1; // The end is inclusive
		std::string a = getBase()->substr(start, end - start);
		int mutationSum = 0;
		// Choose between minK and maxK every time
		int wordSize;
		std::uniform_int_distribution<int> findPoint(0, 20);

		std::reverse(a.begin(), a.end());
		std::string b("");
		for (int seqLen = 0; a.size() > 0;)
		{

			if (a.size() > 20 * init_reg && double(mutationSum) / seqLen < 0.05 && findPoint(generator) == 0)
			{
				numRepeats++;

				int exactLength = exactRepeatDistr(generator);
				while (exactLength > a.size())
				{
					exactLength = exactRepeatDistr(generator);
				}

				wordSize = wordSizeDistr(generator);
				if (exactLength / wordSize < 2)
				{
					std::uniform_int_distribution<int> newExactLengthDist(3, 20);
					exactLength = wordSize * newExactLengthDist(generator);
				}

				int TRSize = (exactLength / wordSize) * wordSize;
				std::string motif = a.substr(a.size() - wordSize, wordSize);
				std::reverse(motif.begin(), motif.end());

				motifList->push_back(motif);

				std::string newSubSeq = motif::FindMotif::makeExact(motif, TRSize);

				HandleSeq handleSeq(1, true);

				auto results = handleSeq.mutate(newSubSeq, mutationRate, baseCountOChrom);
				newSubSeq = results.second;
				newSubSeq = newSubSeq.substr(0, min(newSubSeq.size(), a.size()));

				mutationRate++;
				mutationRate = max(mutationRate % 26, 0);

				TRSize = newSubSeq.size();
				b += newSubSeq;
				repeatList->push_back(newSubSeq);

				for (int j = 0; j < TRSize; j++)
				{
					a.pop_back();
				}
				regionList->push_back(new Location(b.size() - TRSize + start, b.size() + start));
				seqLen += TRSize;
				mutationSum += TRSize;
			}
			else
			{
				char newChar = a.back();
				a.pop_back();
				b += newChar;
				seqLen++;
			}
		}
		replaceSeq(b, start, start + b.size());
		globalMutationSum += mutationSum;
	}
	setK(ceil(log(globalMutationSum / (double)numRepeats) / log(4)) - 1);
}


void ChromosomeSelfTR::printBedData(std::string file, std::string chr)
{
	ofstream outSequence;
	outSequence.open(file.c_str(), ios::out);

	if (motifList->size() != repeatList->size() || repeatList->size() != regionList->size())
	{
		std::cerr << "The lengths aren't equal!" << std::endl;
		throw std::exception();
	}
	outSequence << "Start\tEnd\tMotif\tMutationRate\tTR" << std::endl;
	for (int i = 0; i < motifList->size(); i++)
	{
		outSequence << chr << "\t" << regionList->at(i)->getStart() << "\t" <<
		regionList->at(i)->getEnd() << "\t" <<
		Util::oneDigitToNuc(motifList->at(i)) << "\t" << i % 25 << "\t" <<
		Util::oneDigitToNuc(repeatList->at(i)) << std::endl;
	}

	outSequence.close();
}
