
#include "ChromosomeTR.h"
#include <random>
#include "../motif/FindMotif.h"
#include "../mutate/HandleSeq.h"
#include <random>
#include <math.h>

ChromosomeTR::ChromosomeTR(int nIn, ChromosomeOneDigit *oChromIn,
	char unreadIn, int length, int seed_in) : ChromosomeRandom(nIn, oChromIn, unreadIn, length)
{
	regionList = new std::vector<ILocation *>();
	repeatList = new std::vector<std::string>();
	baseCountOChrom = oChromIn->getBaseCount();
	int sum = 0;
	for(int i = 0; i < baseCountOChrom->size(); i++){
		sum += baseCountOChrom->at(i);
	}
	for(int i = 0; i < baseCountOChrom->size(); i++){
		baseCountOChrom->at(i) = (baseCountOChrom->at(i) * 100) / sum;
	}
}

ChromosomeTR::~ChromosomeTR(){
	Util::deleteInVector(regionList);
	delete repeatList;
	delete regionList;
}

void ChromosomeTR::removeLocations(const std::vector<Location *> *deleteList)
{
	if (deleteList->back()->getStart() > rBase->size() || deleteList->back()->getEnd() > rBase->size())
	{
		std::cerr << "The delete list has an out of bound region" << deleteList->back()->toString() << "!" << std::endl;
		throw std::exception();
	}

	auto region = deleteList->begin();
	std::default_random_engine generator(0);
	std::uniform_int_distribution<int> exactRepeatDistr(0, 3);

	for (int index = 0; index < rBase->size(); index++)
	{

		if (region != deleteList->end())
		{

			if (index < (*region)->getStart())
			{
				index = (*region)->getStart() - 1;
			}
			else if (index > (*region)->getEnd())
			{
				index--;
				region++;
			}
			else if (index >= (*region)->getStart())
			{
				rBase->at(index) = (char)exactRepeatDistr(generator);
			}
		}
		else
		{
			break;
		}
	}
}

void ChromosomeTR::replaceSeq(std::string &subSeq, Location &l)
{
	if (l.getEnd() > rBase->size())
	{
		std::cerr << "Can't replace beyond the random chromosome! ";
		std::cerr << l.getEnd() << " " << rBase->size() << std::endl;
		throw std::exception();
	}
	int originalSize = rBase->size();
	rBase->replace(l.getStart(), l.getEnd() - l.getStart(), subSeq);
	if (rBase->size() < originalSize)
	{
		std::cerr << "The base got smaller!" << std::endl;
		throw std::exception();
	}
	else if (rBase->size() > originalSize)
	{
		std::cerr << "The base got bigger!" << std::endl;
		throw std::exception();
	}
}

void ChromosomeTR::replaceSeq(std::string &subSeq, int start, int end)
{
	Location l(start, end);
	replaceSeq(subSeq, l);
}

std::string ChromosomeTR::getSubStr(int start, int end)
{
	return rBase->substr(start, end - start);
}



void ChromosomeTR::setK(int a)
{
	if (a <= 0)
	{
		std::cerr << "Invalid choice of K (" << a << ")!" << std::endl;
		throw std::exception();
	}
	K = a;
}

int ChromosomeTR::getK()
{
	return K;
}

void ChromosomeTR::printRBase()
{
	for (auto it = rBase->begin(); it != rBase->end(); it++)
	{
		if (*it == (char)0)
		{
			std::cout << 'A';
		}
		else if (*it == (char)1)
		{
			std::cout << 'C';
		}
		else if (*it == (char)2)
		{
			std::cout << 'G';
		}
		else if (*it == (char)3)
		{
			std::cout << 'T';
		}
		else
		{
			std::cerr << "Unexptected sequence! " << (int)*it << std::endl;
			throw std::exception();
		}
	}
	std::cout << std::endl;
}


int ChromosomeTR::size()
{
	return rBase->size();
}

std::vector<ILocation *> *ChromosomeTR::getRegionList()
{
	return regionList;
}

std::vector<std::string> *ChromosomeTR::getRepeatList()
{
	return repeatList;
}
