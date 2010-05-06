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


int main(int argc, char *argv[])
{

	Vector epsilon(57) ;
	Vector bulkx(57) ;
	Vector shearx(57) ;
	Vector epsilonx(57) ;

	for(int i = 1 ; i < epsilon.size() ; i++)
	{
		epsilon[i] = 0.01*i ;
	}

	Properties enu(HOOKE, std::make_pair(70.,0.2)) ;
	Properties kmu = enu.convert(BULK_SHEAR).second ;

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

	}
	
	return 0 ;
}
