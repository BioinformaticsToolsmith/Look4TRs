
#include "ChromosomeReadTR.h"
#include <random>
#include "../motif/FindMotif.h"
#include "../mutate/HandleSeq.h"
#include <random>
#include <math.h>
#include <algorithm>    // std::min


ChromosomeReadTR::ChromosomeReadTR(int nIn, ChromosomeOneDigit *oChromIn,
	char unreadIn, int length, std::string fa_file, std::string bed_file, int seed_in) : ChromosomeTR(nIn, oChromIn, unreadIn, length, seed_in)
{
	regionList = new std::vector<ILocation *>();
	repeatList = new std::vector<std::string>();

	readTRList = oChromIn->getSequenceFromLocations(bed_file);

	shuffle();
}

ChromosomeReadTR::ChromosomeReadTR(int nIn, ChromosomeOneDigit *oChromIn,
	char unreadIn, int length, std::string fa_file, std::string bed_file,
	ChromosomeOneDigit * read_from_chrom, int seed_in) : ChromosomeTR(nIn, oChromIn, unreadIn, length, seed_in)
{
	regionList = new std::vector<ILocation *>();
	repeatList = new std::vector<std::string>();

	readTRList = read_from_chrom->getSequenceFromLocations(bed_file);

	
	shuffle();
}

ChromosomeReadTR::~ChromosomeReadTR()
{
	delete readTRList;
}




std::string ChromosomeReadTR::getRandTR()
{
	std::string result = readTRList->at(TR_file_index);
	TR_file_index++;
	TR_file_index = TR_file_index % readTRList->size();
	return result;
}

void ChromosomeReadTR::shuffle()
{
	int originalSize = rBase->size();

	std::default_random_engine generator(seed);
	seed++;

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

			if (double(mutationSum) / seqLen < 0.05 && findPoint(generator) == 0)
			{
				numRepeats++;

				std::string newSubSeq = getRandTR();


				int TRSize = newSubSeq.size();
				int trimmedTRSize =  TRSize > a.size() ? a.size() : TRSize;
				newSubSeq = newSubSeq.substr(0, trimmedTRSize);

				b += newSubSeq;
				repeatList->push_back(newSubSeq);

				for (int j = 0; j < trimmedTRSize; j++)
				{
					a.pop_back();
				}
				regionList->push_back(new Location(b.size() - trimmedTRSize + start - 1, b.size() + start - 1));
				seqLen += trimmedTRSize;
				mutationSum += trimmedTRSize;
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

	if(rBase->size() != originalSize){
		std::cerr << "Not original size!" << std::endl;
		throw std::exception();
	}
}


void ChromosomeReadTR::printBedData(std::string file, std::string chrom_name)
{
	ofstream outSequence;
	outSequence.open(file.c_str(), ios::out);

	if (repeatList->size() != regionList->size())
	{
		std::cerr << "The lengths aren't equal!" << std::endl;
		throw std::exception();
	}
	for (int i = 0; i < repeatList->size(); i++)
	{
		outSequence << chrom_name << "\t" << regionList->at(i)->getStart() << "\t" <<
		regionList->at(i)->getEnd() << "\t" <<
		Util::oneDigitToNuc(repeatList->at(i)) << std::endl;
	}

	outSequence.close();
}
