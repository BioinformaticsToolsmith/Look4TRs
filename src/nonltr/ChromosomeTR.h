#ifndef _CHROMSOMETR_H_
#define _CHROMSOMETR_H_

#include <string>
#include <random>
#include "../utility/Location.h"
#include "ChromosomeRandom.h"

namespace nonltr {

	class ChromosomeTR : public ChromosomeRandom {

	public:
		ChromosomeTR(int, ChromosomeOneDigit *, char, int, int);
		virtual ~ChromosomeTR();
		void removeLocations(const std::vector<Location *> *);
		virtual void replaceSeq(std::string&, Location&);
		virtual void replaceSeq(std::string&, int, int);
		virtual std::string getSubStr(int, int);
		virtual std::string getRandTR() = 0;
		virtual void shuffle() = 0;
		int size();
		void printRBase();
		void setK(int);
		int getK();
		virtual void printBedData(std::string, std::string chr = "") = 0;
		std::vector<ILocation *> * getRegionList();
		std::vector<std::string> * getRepeatList();
	protected:
		std::vector<ILocation *> * regionList;
		std::vector<std::string> * repeatList;
		std::vector<std::string> * readTRList;
		std::vector<int> * baseCountOChrom;
		int K;
		int seed = 0;

	};


}
#endif
