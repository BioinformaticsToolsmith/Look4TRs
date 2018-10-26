#ifndef _RegionList_H_
#define _RegionList_H_ value

#include <deque>
#include <string>
#include <vector>
#include <unordered_map>
#include "Location.h"

namespace utility {
class RegionList {
private:
	std::vector<std::string> headerList;
	unordered_map<std::string, std::deque<utility::Location *> *> * regionList;
public:
	RegionList();
	RegionList(std::string);
	~RegionList();
	void readFile(std::string);
	void pushLocation(std::string, utility::Location *);
	std::deque<utility::Location *> * getLocationList(std::string);
	bool hasLocationList(std::string);
	std::vector<std::string> getRegionList();

};
}
#endif
