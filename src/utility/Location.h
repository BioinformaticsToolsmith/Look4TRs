/*
 * Location.h
 *
 *  Created on: Dec 19, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef LOCATION_H_
#define LOCATION_H_

#include "ILocation.h"

#include <string>

 using namespace std;

 namespace utility {

 	class Location: public ILocation {
 	private:
 		int start;
 		int end;
 		void initialize(int, int);
 		void check();

 	public:
 		Location(int, int);
 		Location(const ILocation&);
 		Location(const ILocation*);
 		virtual ~Location();

 		int getEnd() const;
 		int getStart() const;
 		void setEnd(int);
 		void setStart(int);
 		int getLength();
 		string toString();

 	};

 	// std::ostream& operator<< (std::ostream& out, Location& a){
 	// 	return out;
 	// }

 }

#endif /* LOCATION_H_ */
