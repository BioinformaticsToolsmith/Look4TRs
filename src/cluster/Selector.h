/*
 * Author: Benjamin T James
 * The Bioinformatics Toolsmith Laboratory, the University of Tulsa
 */

#ifndef SELECTOR_H
#define SELECTOR_H
#include <vector>
#include <unordered_map>
#include "Point.h"

template<class T>
struct pra {
	Point<T>* first;
	Point<T>* second;
	double val;
	pra<T>(Point<T>* f, Point<T>* s, double v) : first(f), second(s), val(v) {};
	pra<T>() {};
};

template<class T>
class Selector {
public:
	Selector(std::vector<Point<T>*> v, size_t sample_size_,size_t max_pts_from_one_) : points(v), sample_size(sample_size_), max_pts_from_one(max_pts_from_one_) {

	}
	void set_training(pair<vector<pra<T> >,
			vector<pra<T> > > &tr) { training = tr; }
	void set_testing(pair<vector<pra<T> >,
			vector<pra<T> > > &tr) { testing = tr; }
	void select(double cutoff);
	static double align(Point<T>*a, Point<T>*b);
	pair<vector<pra<T> >,
	vector<pra<T> > > get_training() const { return training; }
	pair<vector<pra<T> >,
	vector<pra<T> > > get_testing() const { return testing; }
private:
	vector<std::pair<Point<T>*, Point<T>*> > split(double cutoff);

	vector<pra<T> > get_align(vector<std::pair<Point<T>*,Point<T>*> >&) const;

	pair<vector<pra<T> > ,
	vector<pra<T> > > get_labels(vector<pra<T> > &, double cutoff) const;

	const size_t sample_size, max_pts_from_one;
	std::vector<Point<T>*> points;
	pair<vector<pra<T> >,
	vector<pra<T> > > training, testing;
};
#endif
