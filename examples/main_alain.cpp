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

	Material unhydrated("UNHYDRATED") ;
	Material water("WATER") ;

	MeanScheme density(false,true,TAG_DENSITY) ;

	Material cement("CEMENT") ;
	cement = cement + water*"0.6_MASS"*"1._DENSITY" + unhydrated*"1.2_MASS"*"3.1_DENSITY" ;
	cement = cement * "14._YOUNG_MODULUS" * "0.2_POISSON_RATIO" ;
	cement.build(&density) ;
	cement.print() ;

/*
	GeneralConverter volume(TAG_VOLUME) ;
	GeneralConverter fraction(TAG_VOLUME_FRACTION) ;
	AdditionConverter volume_total(TAG_VOLUME) ;
	MoriTanaka elastic ;
	SelfConsistentMultiLayerExpansion asr(0.5) ;


	// initialize all materials
	Material unhydrated ;
	unhydrated.add(TAG_MASS,1.2) ;
	unhydrated.add(TAG_MASS_TOTAL,1.8) ;
	unhydrated.add(TAG_DENSITY,3.1) ;
	unhydrated.findMissing(density.inputList()) ;

	Material water ;
	water.add(TAG_MASS,0.6) ;
	water.add(TAG_MASS_TOTAL,1.8) ;
	water.add(TAG_DENSITY,1.) ;
	water.findMissing(density.inputList()) ;

	Material cement(density.homogenize(water,unhydrated)) ;
	cement.add(TAG_MASS,1.8) ;
	cement.add(TAG_MASS_TOTAL,5.6) ;
	cement.add(TAG_YOUNG_MODULUS,14.) ;
	cement.add(TAG_POISSON_RATIO,0.3) ;
	cement.merge(volume.homogenize(cement)) ;
	cement.findMissing(elastic.inputList()) ;

	Material aggregate ;
	aggregate.add(TAG_YOUNG_MODULUS,59.) ;
	aggregate.add(TAG_POISSON_RATIO,0.3) ;
	aggregate.add(TAG_MASS,3.8) ;
	aggregate.add(TAG_MASS_TOTAL,5.6) ;
	aggregate.add(TAG_DENSITY,2.2) ;
	aggregate.merge(volume.homogenize(aggregate)) ;
	aggregate.findMissing(elastic.inputList()) ;

	Material homogenized(volume_total.homogenize(cement,aggregate)) ;

	cement.add(TAG_VOLUME_TOTAL,homogenized.val(TAG_VOLUME,-1)) ;
	aggregate.add(TAG_VOLUME_TOTAL,homogenized.val(TAG_VOLUME,-1)) ;
	
	homogenized.merge(elastic.homogenize(cement,aggregate)) ;

	Material asrGel ;
	asrGel.add(TAG_YOUNG_MODULUS,31.) ;
	asrGel.add(TAG_POISSON_RATIO,0.28) ;
	asrGel.add(TAG_VOLUME_TOTAL,cement.val(TAG_VOLUME_TOTAL,-1)) ;

	// case 1 : asrGel + aggregate + cement


	Material gel(asrGel) ;
	Material agg(aggregate) ;
	Material non(aggregate) ;
	Material cem(cement) ;
	Material hom(homogenized) ;

	gel.rename("GEL") ;
	agg.rename("AGG") ;
	non.rename("NON") ;
	cem.rename("CEM") ;
	hom.rename("HOM") ;

	gel.add(TAG_VOLUME,agg.val(TAG_VOLUME,-1)*0.01*0.25) ;
	agg.replace(Properties(TAG_VOLUME,agg.val(TAG_VOLUME,-1)*0.99*0.25)) ;
	non.replace(Properties(TAG_VOLUME,non.val(TAG_VOLUME,-1)*0.75)) ;
	cem.replace(Properties(TAG_VOLUME,cem.val(TAG_VOLUME,-1))) ;
	hom.add(TAG_VOLUME,cem.val(TAG_VOLUME_TOTAL,-1)*0.5) ;

	std::vector<Material> mat ;
	mat.push_back(gel) ;
	mat.push_back(agg) ;
	mat.push_back(cem) ;
	mat.push_back(non) ;
	mat.push_back(hom) ;
	Material tot(volume_total.homogenize(mat)) ;
	tot.print() ;
	for(size_t i = 0 ; i < mat.size() ; i++)
	{
		mat[i].replace(Properties(TAG_VOLUME_TOTAL,tot.val(TAG_VOLUME,-1))) ;	
		mat[i].findMissing(asr.inputList()) ;
		mat[i].print() ;
	}

	Material expanded(asr.homogenize(mat)) ;

	expanded.print() ;*/

	return 0 ;
}
