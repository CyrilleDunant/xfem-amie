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
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/parallel_behaviour.h"
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
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"

#include <fstream>
#include <omp.h>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

void checkTriangle(DelaunayTriangle * trg)
{
	Point p ;
//		triangles[i]->getBehaviour()->getTensor(p).print() ;
	Matrix m = trg->getBehaviour()->getViscousTensor(p) ;
	int k = 0 ;
	for(size_t r = 0 ; r < 3 ; r++)
	{
		for(size_t c = 0 ; c < 3 ; c++)
		{
			if(m[r][c] > 1)
				k++ ;
			if(m[r+3][c] > 1)
				k++ ;
			if(m[r][c+3] > 1)
				k++ ;
		}
	}
	if(k > 0)
		std::cout << k << " bad components" << std::endl ;
}

int main(int argc, char *argv[])
{
	double deltat = atof(argv[1]) ;

	Sample box(nullptr, 0.2,0.4,0.,0.) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(16) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	double time_step = deltat ;
	F.setDeltaTime(time_step) ;
	F.setMinDeltaTime(1e-9) ;
	F.setSamplingRestriction( SAMPLE_RESTRICT_16 ) ;

	LogarithmicCreepWithExternalParameters paste("young_modulus = 20e9, poisson_ratio = 0.3, creep_modulus = 50e9, creep_poisson = 0.3, creep_characteristic_time = 0.1") ;
	LogarithmicCreepWithExternalParameters aggregates("young_modulus = 60e9, poisson_ratio = 0.2") ;

	box.setBehaviour( &paste );

	std::vector<Feature *> inc = PSDGenerator::get2DConcrete(&F, &aggregates, 10, 0.075, 0.0002, new GranuloFromCumulativePSD("../examples/data/bengougam/granulo_luzzone", CUMULATIVE_PERCENT), CIRCLE, 1., M_PI, 1000000, 0.8, new Rectangle(0.5,1.,0,0)) ;
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		if(inc[i]->intersects(dynamic_cast<Rectangle*>(&box)))
			F.setSamplingFactor( inc[i], 3. ) ;
		else
			F.setSamplingFactor( inc[i], 2. ) ;
	}
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;

	F.step() ;

	int pasteCachePosition = F.get2DMesh()->generateCache( new ViscosityHigherThanElementChecker(1.,3,3) ) ;
	int aggregatesCachePosition = F.get2DMesh()->generateCache( new ViscosityLowerThanElementChecker(1.,3,3) ) ;

	double agg = F.get2DMesh()->getArea(aggregatesCachePosition) ;
	std::cout << "effective surface covered by aggregates: " << agg << std::endl ;

	BoundingBoxDefinedBoundaryCondition * stress = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, 1e6) ;

	F.addBoundaryCondition(stress) ;

	std::fstream out ;
	std::string file = "superposition_bengougam_tension_" ;
	file.append(argv[1]) ;
	out.open(file.c_str(), std::ios::out) ;
	int i = 0 ;
	bool down = false ;
	bool up = false ;
	int incr = 10 ;

	while(F.getCurrentTime() < 1e5)
	{
		F.setDeltaTime(deltat) ;
		bool goOn = F.step() ;
		if(!goOn)
		{
			TriangleWriter trg("tata", &F, 1.) ;
			trg.getField(STRAIN_FIELD) ;
			trg.getField(PRINCIPAL_REAL_STRESS_FIELD) ;
			trg.getField(TWFT_STIFFNESS) ;
			trg.getField(TWFT_VISCOSITY) ;
			trg.write() ;
			return 0 ;
		}

	Vector strain = F.get2DMesh()->getField(STRAIN_FIELD, -1, 1.) ;
	Vector stress = F.get2DMesh()->getField(REAL_STRESS_FIELD, -1, 1.) ;
	Vector stress_paste = F.get2DMesh()->getField(REAL_STRESS_FIELD, pasteCachePosition, -1, 1.) ;
	Vector stress_aggregates = F.get2DMesh()->getField(REAL_STRESS_FIELD, aggregatesCachePosition, -1, 1.) ;


	out << F.getCurrentTime() << "\t" << strain[1]*1e6 << "\t" << strain[4]*1e6 << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress_paste[0] << "\t" << stress_paste[1] <<  "\t" << stress_aggregates[0] << "\t" << stress_aggregates[1] << std::endl ;
	std::cout << F.getCurrentTime() << "\t" << strain[0]*1e6 << "\t" << strain[1]*1e6 << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress_paste[0] << "\t" << stress_paste[1] <<  "\t" << stress_aggregates[0] << "\t" << stress_aggregates[1] << std::endl ;
	i++ ;
    }

    return 0 ;

}

