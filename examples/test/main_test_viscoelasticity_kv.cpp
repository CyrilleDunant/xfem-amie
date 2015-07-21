// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../main.h"
#include "../../features/features.h"
#include "../../physics/viscoelasticity.h"
#include "../../physics/materials/paste_behaviour.h"
#include "../../utilities/writer/triangle_writer.h"
#include "../../features/sample.h"

#include <fstream>
#include <omp.h>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

int main(int argc, char *argv[])
{
        Sample rect(nullptr, 0.04,0.04,0,0) ;
	Matrix C = ElasticOnlyPasteBehaviour(10e9).param ;
	rect.setBehaviour(new Viscoelasticity(GENERALIZED_KELVIN_VOIGT, C, C, C*25 ) ) ;

	FeatureTree f(&rect) ;
	f.setSamplingNumber(0) ;
	int deltaTime = 1 ;
        f.setDeltaTime(1) ;
	
	f.step() ;

	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT_AFTER) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM_AFTER ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -1e6 ) ) ;

	std::ofstream out ;
	out.open("../examples/test/test_viscoelasticity_kv_current", std::ios::out) ;

	for(size_t i = 0 ; i < 30 ; i++)
	{
                f.setDeltaTime(deltaTime++) ;
		f.step() ;
		out << f.getCurrentTime() << "\t" << f.getAverageField( REAL_STRESS_FIELD, -1,1 )[1] << "\t" << f.getAverageField( STRAIN_FIELD, -1,1 )[1] << std::endl ;
	}

	return 0 ;
}
