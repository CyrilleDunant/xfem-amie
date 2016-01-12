// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../features/sample.h"
#include "../../features/polygonSample.h"
#include "../../physics/stiffness.h"
#include "../../physics/logarithmic_creep_with_external_parameters.h"
#include "../../physics/material_laws/mechanical_material_laws.h"
#include "../../physics/material_laws/material_laws.h"
#include "../../utilities/parser.h"
#include "../../utilities/itoa.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../geometry/sampler/gradient_sampler.h" 
#include "../../geometry/sampler/regular_sampler.h" 
#include "../../utilities/mineral.h" 


#include <dirent.h>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;


int main( int argc, char *argv[] )
{
        Sample rect(nullptr, 0.04,0.04,0,0) ;
	rect.setBehaviour(new Stiffness( 10e9, 0.2 ) ) ;
        
	FeatureTree f(&rect) ;
	f.setSamplingNumber(16) ;

	Inclusion * inc = new Inclusion( 0.01,0,0 ) ;
	inc->setBehaviour(new Stiffness( 20e9, 0.2 ) ) ;
//	f.addFeature( &rect, inc ) ;

	f.setSampler( &rect, new GradientSampler( Point(0.02,0.01), Point(-0.02,0), 0.5, 2.) ) ;
	
	f.step() ;

	TriangleWriter trg( "toto", &f, 1. ) ;
	trg.getField( "C11" ) ;
	trg.write() ;

	
	return 0 ;
}
