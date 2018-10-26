#ifndef _CHROMSOMETR_H_
#define _CHROMSOMETR_H_

#include <string>
#include "../utility/Location.h"
#include "ChromosomeRandom.h"

namespace nonltr {

	class ChromosomeTR : public ChromosomeRandom {

	public:
		ChromosomeTR(int, ChromosomeOneDigit*, char, int, int, int, int);
		~ChromosomeTR();
		void removeLocations(const std::vector<Location *> *);
		void replaceSeq(std::string&, Location&);
		void replaceSeq(std::string&, int, int);
		std::string getSubStr(int, int);
		void shuffle();
		int size();
		void printRBase();
		int getK();
		std::vector<ILocation *> * getRegionList();
		std::vector<std::string> * getRepeatList();
	private:
		std::vector<ILocation *> * regionList;
		std::vector<std::string> * repeatList;
		void setK(int);
		int init_reg;
		int minMotif, maxMotif;
		int K;
	};

}
#endif