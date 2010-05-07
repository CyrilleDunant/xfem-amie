// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.cpp"
#include "../physics/homogenization/properties_base.h"
#include "../physics/homogenization/scheme_base.h"
#include "../physics/homogenization/cracked_homogenization.h"
#include "../physics/homogenization/elastic_homogenization.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#define DEBUG 

using namespace Mu ;

int main(int argc, char *argv[])
{

	Vector epsilon(100) ;
	Vector bulkx(100) ;
	Vector shearx(100) ;
	Vector epsilonx(100) ;

	for(int i = 1 ; i < epsilon.size() ; i++)
	{
		epsilon[i] = 0.01*i ;
	}

	Material agg ;
	agg.push_back(Properties(TAG_YOUNG_MODULUS,70.)) ;
	agg.push_back(Properties(TAG_POISSON_RATIO,0.2)) ;

	agg.push_back(Properties(TAG_ELLIPSE_A,1.)) ;
	agg.push_back(Properties(TAG_ELLIPSE_B,0.75)) ;

	agg.push_back(Properties(TAG_CRACK_DENSITY,0.01)) ;

	BudianskyScheme bud ;
	agg.findMissing(bud.inputList()) ;

	agg.print() ;

	Material agg_cracked(bud.homogenize(agg)) ;

	std::cout << "---------" << std::endl ;

	agg_cracked.print() ;


/*

	bulkx[0] = kmu.getValue(0) ;
	shearx[0] = kmu.getValue(1) ;

	Properties eps(CRACK_DENSITY,0.) ;
	Properties ell(ELLIPSE_SHAPE,std::make_pair(1.,0.75)) ;

	for(int i = 1 ; i < epsilon.size() ; i++)
	{
		eps.setValue(0,epsilon[i]) ;

		Material agg(kmu) ;
		agg.push_back(eps) ;
		agg.push_back(ell) ;

		std::vector<Material> mat ;
		mat.push_back(agg) ;

		Material agg_cracked = BudianskyScheme().apply(mat).second ;

		bulkx[i] = agg_cracked[0].getValue(0) ;
		shearx[i] = agg_cracked[0].getValue(1) ;
		epsilonx[i] = agg_cracked[1].getValue(0) ;

		std::cout << i << ";" << epsilonx[i] << ";" << bulkx[i] << ";" << shearx[i] << std::endl ;

	}*/
	
	return 0 ;
}
