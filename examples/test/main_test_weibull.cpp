// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/stiffness.h"
#include "../../features/sample.h"
#include "../../utilities/parser.h"

#include <fstream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
	CommandLineParser parser("Test the randomness of a weibull-distributed stiffness") ;
	parser.addFlag("--renew-base", false, "renew the base of results") ;
	parser.addValue("--seed", 1, "random seed") ;
	parser.parseCommandLine(argc, argv) ;
	bool renew = parser.getFlag("--renew-base") ;
	std::srand( parser.getValue("--seed") ) ;


        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new WeibullDistributedElasticStiffness(20e9, 0.2, 0.2) ) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(64) ;
        f.setDeltaTime(1) ;
	
	f.step() ;

        std::vector<DelaunayTriangle *> trg = f.get2DMesh()->getConflictingElements(dynamic_cast<Rectangle *>(&rect)) ;
	Point p ;

	Vector area(trg.size()) ;
	Vector stiffness(trg.size()) ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		stiffness[i] = trg[i]->getBehaviour()->getTensor(p)[0][0] ;
		area[i] = trg[i]->area() ;
	}
	double total = area.sum() ;

	double min = stiffness.min() ;
	double max = stiffness.max() ;
	double avg = 0. ;
	for(size_t i = 0 ; i < trg.size() ; i++)
		avg += stiffness[i]*area[i] ;
	avg /= total ;

	double stdev = 0. ;
	for(size_t i = 0 ; i < trg.size() ; i++)
		stdev += area[i]*(stiffness[i]-avg)*(stiffness[i]-avg) ;
	stdev /= total ;
	stdev = sqrt(stdev) ;

	std::ofstream out ;
	if(renew)
		out.open("../examples/test/test_weibull_base", std::ios::out) ;
	else
		out.open("../examples/test/test_weibull_current", std::ios::out) ;


	out << min << "\t" << max << std::endl ;
	out << avg << "\t" << stdev << std::endl ;

	Vector histogram(10) ;
	for(size_t i = 0 ; i < trg.size() ; i++)
	{
		for(size_t j = 0 ; j < histogram.size() ; j++)
		{
			if( stiffness[i] < 16e9+j*1e9)
				histogram[j] += area[i] ;
		}
	}

	for(size_t j = 0 ; j < histogram.size() ; j++)
		out << 16e9+j*1e9 << "\t" << histogram[j]/total << std::endl ;

	return 0 ;
}

