/*
 * ChromosomeRandom.cpp
 *
 *  Created on: Feb 4, 2013
 *  Author: Hani Zakaria Girgis, PhD
 *  Modified by Robert
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <random>

#include "KmerHashTable.h"
#include "../exception/InvalidInputException.h"
#include "../exception/InvalidStateException.h"
#include "../utility/Util.h"
#include "ChromosomeRandom.h"

 using namespace std;
 using namespace exception;
 using namespace utility;

 namespace nonltr {

 	ChromosomeRandom::ChromosomeRandom(int nIn, ChromosomeOneDigit* oChromIn,
 		char unreadIn) {
 		canDelete = false;

 		initializer(nIn, oChromIn, unreadIn, oChrom->getBase()->size());
 	}

 	ChromosomeRandom::ChromosomeRandom(int nIn, ChromosomeOneDigit* oChromIn,
 		char unreadIn, int length) {
 		canDelete = false;
 		initializer(nIn, oChromIn, unreadIn, length);
 	}

 	void ChromosomeRandom::initializer(int nIn, ChromosomeOneDigit* oChromIn,
 		char unreadIn, int length) {
	// Check the order
 		if (nIn < 0) {
 			string msg("The Markov order must be non-negative. ");
 			msg.append("The order received is: ");
 			msg.append(Util::int2string(nIn));
 			msg.append(".");
 			throw InvalidInputException(msg);
 		}

	// n here is the length of the word, i.e. the order + 1
 		n = nIn + 1;
 		unread = unreadIn;

	// Initialize the random sequence
 		oChrom = oChromIn;
 		oBase = oChrom->getBase();
 		randLength = length;
 		rBase = new string(randLength, unread);

	// Initialize the table
 		table = new KmerHashTable<int, double>(n, 1);

	//	Print codes for translating the one digit form of the sequence to the
	//	lettered A,C,G,T form.
 		printCodes = new map<char, char>();
 		printCodes->insert(map<char, char>::value_type((char) 0, 'A'));
 		printCodes->insert(map<char, char>::value_type((char) 1, 'C'));
 		printCodes->insert(map<char, char>::value_type((char) 2, 'G'));
 		printCodes->insert(map<char, char>::value_type((char) 3, 'T'));
 		printCodes->insert(map<char, char>::value_type('n', 'n'));
 		printCodes->insert(map<char, char>::value_type('N', 'N'));

 		randSegmentList = new vector<vector<int> *>();

 		countWords();
 		convertToProbabilities();
 		makeSegmentList();
 		generateRandomSequence();
 	}

 	ChromosomeRandom::~ChromosomeRandom() {
 		printCodes->clear();
 		delete printCodes;
 		delete table;
 		delete rBase;

 		Util::deleteInVector(randSegmentList);
 		delete randSegmentList;
 	}

 	void ChromosomeRandom::countWords() {
 		const vector<vector<int> *> * segmentList = oChrom->getSegment();
 		int segmentCount = segmentList->size();

 		const char * oCharBase = oBase->c_str();

 		for (int i = 0; i < segmentCount; i++) {
 			table->wholesaleIncrement(oCharBase, segmentList->at(i)->at(0),
 				segmentList->at(i)->at(1) - n + 1);
 		}
 	}

 	void ChromosomeRandom::convertToProbabilities() {
 		int alphaCount = 4;
 		int tableSize = table->getMaxTableSize();
 		for (int i = 0; i < tableSize; i += alphaCount) {
 			double sum = 0.0;
 			for (int j = 0; j < alphaCount; j++) {
 				sum += table->valueOf(i + j);
 			}
 			for (int j = 0; j < alphaCount; j++) {
 				table->insert(i + j, table->valueOf(i + j) / sum);
 			}
 		}
 	}

/**
 * Make a list of segment based on the desired length and the original chromosome
 */
 void ChromosomeRandom::makeSegmentList() {
	// Get the original segments
 	const vector<vector<int> *> * segmentList = oChrom->getSegment();
 	int segmentCount = segmentList->size();

 	int total = 0;
 	for (int i = 0; i < segmentCount; i++) {
 		int s = segmentList->at(i)->at(0);
 		int e = segmentList->at(i)->at(1);

 		if (e < randLength) {
 			randSegmentList->push_back(new vector<int>(*segmentList->at(i)));
 		} else {
 			if (s < randLength) {
 				vector<int> * segment = new vector<int>();
 				segment->push_back(s);
 				segment->push_back(randLength - 1);

 				randSegmentList->push_back(segment);
 			}
 			break;
 		}
 	}

 	// Handle the case when the real chromsome is shorter than
 	// the desired length
 	auto lastSegment = randSegmentList->back();
 	if(lastSegment->at(1) < randLength-1){
 		lastSegment->at(1) = randLength-1;
 	}


	//Post condition
 	for (int i = 0; i < randSegmentList->size(); i++) {
 		auto seg = randSegmentList->at(i);
 		Location loc(seg->at(0), seg->at(1));
 	}

	// Post condition
 	for (int i = 1; i < randSegmentList->size(); i++) {
 		auto seg1 = randSegmentList->at(i - 1);
 		auto seg2 = randSegmentList->at(i);
 		if (seg2->at(0) <= seg1->at(1)) {
 			cerr << "ChromosomeRandom::makeSegmentList(): "
 			<< "Error while constructing segments" << endl;
 			throw std::exception();
 		}
 	}

 }

 void ChromosomeRandom::generateRandomSequence() {
	// Get the original sequence

	// Alphabet count
 	int alphaCount = 4;

	//CHANGE THIS AFTER TESTING
	// srand(1);
 	std::default_random_engine generator;

	// Generate random segments
 	int segmentCount = randSegmentList->size();
 	for (int i = 0; i < segmentCount; i++) {
 		int s = randSegmentList->at(i)->at(0);
 		int e = randSegmentList->at(i)->at(1);

 		if (e - s + 1 > n) {
 			string order("");

			// The first order is based on the original sequence.
 			for (int w = s; w < s + n - 1; w++) {
 				(*rBase)[w] = oBase->at(w);
 				order.append(1, oBase->at(w));
 			}

 			for (int h = s + n - 1; h <= e; h++) {
				// Subsequent orders are based on the random sequence.
 				vector<vector<int> > lottery;
 				int chanceSoFar = 0;
 				for (int k = 0; k < alphaCount; k++) {
 					string temp = order;
					temp.append(1, k); // Change to k
					if (table->valueOf(temp.c_str()) > 0) {
						int periodStart = chanceSoFar;
						int periodEnd = periodStart
						+ (100 * table->valueOf(temp.c_str()));
						chanceSoFar = periodEnd + 1;
						vector<int> entry;
						entry.push_back(k);
						entry.push_back(periodStart);
						entry.push_back(periodEnd);
						lottery.push_back(entry);
					} else {
						string msg("This word must exist in the table: ");
						msg.append(temp);
						msg.append(".");
						throw InvalidStateException(msg);
					}
				}

				if (lottery.size() > 0) {
					// int randInt = rand() % chanceSoFar;
					std::uniform_int_distribution<int> distribution(0,chanceSoFar - 1);
					int randInt = distribution(generator);



					for (int tt = 0; tt < alphaCount; tt++) {
						vector<int> entry = lottery.at(tt);
						if (randInt >= entry.at(1) && randInt <= entry.at(2)) {
							(*rBase)[h] = entry.at(0);
							order.append(1, entry.at(0));
							break;
						}
					}
					lottery.clear();
				} else {
					string msg("The lottery vector cannot be empty.");
					throw InvalidStateException(msg);
				}
				order = order.substr(1);
			}
		}
	}

	// Make sure that the generated sequence has the desired length
	if (randLength != rBase->size()) {
		cerr << "The original sequence and the random sequence ";
		cerr << "do not have the same size." << endl;
		cerr << "The desired size is: " << randLength << endl;
		cerr << "Generated sequence size is: " << rBase->size() << endl;
	}
}


/**
 * Returns the segments of the original chromosome
 */
 const vector<vector<int> *> * ChromosomeRandom::getSegment() {
 	return randSegmentList;
 }


/**
 * Returns the random sequence
 */
 const string* ChromosomeRandom::getBase() {
 	return rBase;
 }

/**
 * Returns the header indicating the order of the Markov chain
 */
 string ChromosomeRandom::getHeader() {
 	string header = oChrom->getHeader();
 	return header;
 }

 void ChromosomeRandom::printEffectiveSequence(string outputFile) {
 	int totalSize = rBase->size();
 	string * effectiveRBase = new string("");
 	for (int i = 0; i < totalSize; i++) {
 		char b = rBase->at(i);
 		if (b != unread) {
 			effectiveRBase->append(1, b);
 		}
 	}

	// Make sure that the effective sequence is shorter than the original
	// length
 	if (effectiveRBase->size() > totalSize) {
 		cerr << "The effective length must be <= the original length." << endl;
 		cerr << "Generated sequence size is: " << totalSize << endl;
 		cerr << "The effective size is: " << effectiveRBase->size() << endl;

 	}

 	printSequence(outputFile, effectiveRBase);

 	delete effectiveRBase;
 }

 void ChromosomeRandom::printSequence(string outputFile) {
 	printSequence(outputFile, rBase);
 }

 void ChromosomeRandom::printSequence(string outputFile, string * baseToPrint) {
 	cout << "Printing chromosome to file ..." << endl;
 	ofstream outSequence;
 	outSequence.open(outputFile.c_str(), ios::out);

 	int step = 50;

 	outSequence << getHeader() << endl;
 	int len = baseToPrint->size();

 	for (int i = 0; i < len; i = i + step) {
 		int e = (i + step - 1 > len - 1) ? len - 1 : i + step - 1;
 		for (int k = i; k <= e; k++) {
 			outSequence << printCodes->at(baseToPrint->at(k));
 		}
 		outSequence << endl;
 	}
 	outSequence << endl;

 	outSequence.close();
 }

} /* namespace nonltr */
