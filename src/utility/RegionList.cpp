#include "RegionList.h"
#include "Util.h"

namespace utility {
RegionList::RegionList() {
	regionList = new unordered_map<std::string,
			std::deque<utility::Location *> *>();

}

RegionList::RegionList(std::string file) {
	Util::checkFile(file);
	regionList = new unordered_map<std::string,
			std::deque<utility::Location *> *>();
	Util::readCoordinates(file, regionList);

}

RegionList::~RegionList() {
	for (auto it = regionList->begin(); it != regionList->end(); ++it) {
		while (it->second->size() > 0) {
			delete it->second->at(0);
			it->second->pop_front();
		}
		delete it->second;
	}
	delete regionList;
}

void RegionList::readFile(std::string file) {
	Util::checkFile(file);
	Util::readCoordinates(file, regionList);
}

void RegionList::pushLocation(std::string header, utility::Location * loc) {
	if (regionList->count(header) == 0) {
		std::deque<utility::Location *> * vec = new std::deque<
				utility::Location *>();
		regionList->emplace(header, vec);
	}
	regionList->at(header)->push_back(loc);
}

std::deque<utility::Location *> * RegionList::getLocationList(
		std::string header) {
	return regionList->at(header);
}

bool RegionList::hasLocationList(std::string header) {
	return regionList->count(header) != 0;
}

std::vector<std::string> RegionList::getRegionList() {
	headerList.clear();
	for (auto it = regionList->begin(); it != regionList->end(); ++it) {

		headerList.push_back(it->first);
	}
	return headerList;
}

}
