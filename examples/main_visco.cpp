// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/orthotropicstiffness.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
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
#include "../utilities/parser.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"

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
	CommandLineParser parser("2D viscoelastic simulation") ;
	parser.addString("--paste-behaviour", "", "path to a *.ini file containing the behaviour of the cement paste") ;
	parser.addString("--aggregate-behaviour", "", "path to a *.ini file containing the behaviour of the aggregates") ;
	parser.addValue("--duration", 365, "duration of the creep simulation") ;
	parser.addValue("--set-delta-time-increment", 0.1, "increases the time step at each step") ;
	parser.parseCommandLine( argc, argv ) ;

	std::string iniPaste = parser.getString("--paste-behaviour") ;
	std::string iniAgg = parser.getString("--aggregate-behaviour") ;
	double maxTime = parser.getValue("--duration") ;
	double incr = parser.getValue("--set-delta-time-increment") ;

	Sample box(nullptr, 0.08,0.08,0.,0.) ;

	Form * paste = new ViscoElasticOnlyPasteBehaviour() ;
	if(iniPaste.size() > 0)
		paste = ConfigParser::getBehaviour( iniPaste, paste, SPACE_TWO_DIMENSIONAL, true ) ;

	Form * aggregate = new ViscoElasticOnlyAggregateBehaviour() ;
	if(iniAgg.size() > 0)
		aggregate = ConfigParser::getBehaviour( iniAgg, aggregate, SPACE_TWO_DIMENSIONAL, true ) ;
	

	FeatureTree F(&box) ;
	box.setBehaviour( paste ) ;

	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setSamplingNumber(128) ;
	F.setDeltaTime(1) ;
	F.setMinDeltaTime(1e-9) ;
	F.setSamplingRestriction( 4 ) ;

	std::vector<Feature *> finc = PSDGenerator::get2DConcrete( &F, aggregate, 6000, 0.008, 0.000001, new PSDBolomeA(), new InclusionGenerator(), 6000*10) ; 

	parser.setFeatureTree( &F ) ;

	double d = F.getDeltaTime() ;

	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_XI, BOTTOM_LEFT_AFTER ) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_ETA, BOTTOM_AFTER ) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6 ) ) ;

	while(F.getCurrentTime() < maxTime)
	{
		F.step() ;
		Vector strain = F.getAverageField( GENERALIZED_VISCOELASTIC_STRAIN_FIELD, -1, 1 ) ;
		Vector stress = F.getAverageField( REAL_STRESS_FIELD, -1, 1 ) ;
		std::cout << F.getCurrentTime() << "\t" ;
		for(size_t j = 0 ; j < strain.size() ; j++)
			std::cout << strain[j] << "\t" ;
		for(size_t j = 0 ; j < stress.size() ; j++)
			std::cout << stress[j] << "\t" ;
		std::cout << std::endl ;

		d += incr ;
		F.setDeltaTime( d ) ;

	}
	return 0 ;
}

