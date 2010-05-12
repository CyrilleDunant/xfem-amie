// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.cpp"
#include "../physics/homogenization/properties_base.h"
#include "../physics/homogenization/converter.h"
#include "../physics/homogenization/scheme_base.h"
#include "../physics/homogenization/cracked_homogenization.h"
#include "../physics/homogenization/elastic_homogenization.h"
#include "../physics/homogenization/expansion_homogenization.h"
#include "../physics/homogenization/octave_manager.h"

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
	BudianskyDryCrackScheme crack ;
	MoriTanaka elastic ;
	HashinScheme alpha ;
	GeneralConverter young(TAG_YOUNG_MODULUS) ;
	GeneralConverter poisson(TAG_POISSON_RATIO) ;
	GeneralConverter volume(TAG_VOLUME) ;
	AdditionConverter total(TAG_VOLUME_TOTAL) ;
	GeneralConverter fraction(TAG_VOLUME_FRACTION) ;


	Material dummy(MAT_DUMMY) ;
	dummy.push_back(Properties(TAG_ELLIPSE_A,1.)) ;
	dummy.push_back(Properties(TAG_ELLIPSE_B,0.75)) ;
	dummy.push_back(Properties(TAG_CRACK_DENSITY,0.1)) ;
	dummy.findMissing(crack.inputList()) ;

	Material dummy_cracked(crack.homogenize(dummy)) ;
	double factor = dummy_cracked.val(TAG_CRACK_DENSITY,-1) / dummy.val(TAG_CRACK_DENSITY,-1) ;
//	double factor = 1. ;
	std::cout << dummy_cracked.val(TAG_CRACK_DENSITY,-1) << std::endl ;

	Vector reacted(300) ;
	Vector expansion(300) ;
	Vector bulk(300) ;

	for(int i = 0 ; i < 301 ; i++)
	{
		// reset all
		crack.reset() ;
		elastic.reset() ;
		alpha.reset() ;
		young.reset() ;
		poisson.reset() ;
		volume.reset() ;
		total.reset() ;
		fraction.reset() ;

		// initialize
		Material aggregate(MAT_AGGREGATE) ;
		Material cement(MAT_CEMENT) ;

		// mass
		aggregate.push_back(Properties(TAG_MASS,3.789)) ;
		cement.push_back(Properties(TAG_MASS,1.844)) ;

		// volume
		aggregate.merge(Material(volume.homogenize(aggregate))) ;
		cement.merge(Material(volume.homogenize(cement))) ;

		// volume total
		aggregate.merge(Material(total.homogenize(cement,aggregate))) ;
		cement.merge(Material(total.homogenize(cement,aggregate))) ;

		// volume fraction
		aggregate.merge(Material(fraction.homogenize(aggregate))) ;
		cement.merge(Material(fraction.homogenize(cement))) ;

		// elastic properties
		aggregate.findMissing(elastic.inputList()) ;
		cement.findMissing(elastic.inputList()) ;

		// crack aggregate
		aggregate.push_back(Properties(TAG_ELLIPSE_A,1.)) ;
		aggregate.push_back(Properties(TAG_ELLIPSE_B,0.75)) ;
		aggregate.push_back(Properties(TAG_CRACK_DENSITY,(double) i / 100. / 100. / factor)) ;
		aggregate.findMissing(crack.inputList()) ;
		aggregate.merge(Material(crack.homogenize(aggregate))) ;
		double epsilon = aggregate.val(TAG_CRACK_DENSITY,-1) ;

		// expansion coefficient
		aggregate.push_back(Properties(TAG_EXPANSION_COEFFICIENT,epsilon/3)) ;
		cement.push_back(Properties(TAG_EXPANSION_COEFFICIENT,0.)) ;

		// expansion
		Material concrete(elastic.homogenize(cement,aggregate)) ;

		// elastic properties
		concrete.merge(Material(alpha.homogenize(cement,aggregate))) ;
		concrete.merge(Material(young.homogenize(concrete))) ;
		concrete.merge(Material(poisson.homogenize(concrete))) ;

		reacted[i] = epsilon ;
		expansion[i] = concrete.val(TAG_EXPANSION_COEFFICIENT,-1) ;
		bulk[i] = concrete.val(TAG_YOUNG_MODULUS,-1) ;

		if(i == 10)
			std::cout << aggregate.val(TAG_VOLUME_FRACTION,-1) << ";" ;

	}

	std::cout << expansion[300] << std::endl ;
	
	std::string resultfile = "result.m" ;
	OctaveManager result(resultfile.c_str(), std::ios::out) ;

	result.writeArray("R",reacted) ;
	result.writeArray("X",expansion) ;
	result.writeArray("E",bulk) ;

	result.writePlot("R","X",1) ;
	result.writePlot("R","E",2) ;

	result.close() ;

	return 0 ;
}
