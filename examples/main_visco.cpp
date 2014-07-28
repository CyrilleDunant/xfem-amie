// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/kelvinvoight.h"
#include "../physics/maxwell.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/parallel_behaviour.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/orthotropicstiffness.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/homogenization/homogenization_base.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/features.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/itoa.h"

#include <fstream>
#include <omp.h>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Mu ;

int main(int argc, char *argv[])
{
    double timestep = 1. ;
    double appliedLoadEta = -1e6 ;

    Sample box(nullptr, 0.1,0.1,0.,0.) ;

    FeatureTree F(&box) ;
    F.setSamplingNumber(4) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
    F.setDeltaTime(timestep) ;
    F.setMinDeltaTime(timestep*1e-9) ;

    ElasticOnlyPasteBehaviour dummy ;
    LogarithmicCreep creep( dummy.param,dummy.param*2, 10 ) ;
    box.setBehaviour( &creep );

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0, 2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;

    F.step() ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, appliedLoadEta));

    std::string filename = "logcreep" ;
    std::ofstream output ;
    output.open( filename.c_str(), std::ios::out ) ;

    while(F.getCurrentTime() < 10000.)
	{
        F.setDeltaTime(++timestep) ;
        F.setMinDeltaTime(timestep*1e-9) ;
        F.step() ;

        Vector strain = F.getAverageField(STRAIN_FIELD, -1, 1.) ;
        Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1.) ;
        output << F.getCurrentTime() << "\t" << strain[1] << "\t" << stress[1] << std::endl ;
	}

    return 0 ;
}

