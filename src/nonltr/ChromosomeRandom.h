/*
 * ChromosomeRandom.h
 *
 *  Created on: Feb 4, 2013
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef CHROMOSOMERANDOM_H_
#define CHROMOSOMERANDOM_H_

#include <map>

#include "IChromosome.h"
#include "ChromosomeOneDigit.h"
#include "KmerHashTable.h"

namespace nonltr {

class ChromosomeRandom: public nonltr::ChromosomeOneDigit {

private:
	int n;
	int randLength;
	char unread;
	bool canDelete;
	IChromosome * oChrom;
	const string * oBase;
	bool isChromOneDigit;
	map<char, char> * printCodes;
	KmerHashTable<int, double> * table;
	ChromosomeOneDigit* oChromOneDigit;
	vector<vector<int> *> * randSegmentList;


	void countWords();
	void convertToProbabilities();
	void printTable();
	void makeSegmentList();
	void generateRandomSequence();
protected:
	string * rBase;

public:
	ChromosomeRandom(int, ChromosomeOneDigit*, char);
	ChromosomeRandom(int, ChromosomeOneDigit*, char, int);
	void initializer(int, ChromosomeOneDigit*, char, int);
	virtual ~ChromosomeRandom();

	virtual const string* getBase();
	virtual const vector<vector<int> *> * getSegment();
	virtual string getHeader();
	virtual void printSequence(string);
	void printSequence(string, string *);
	void printEffectiveSequence(string);
};

} /* namespace nonltr */
#endif /* CHROMOSOMERANDOM_H_ */
