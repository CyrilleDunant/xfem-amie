// Author: Alain Giorla <alain.b.giorla@gmail.com>, (C) 2005-2014
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

int main(int argc, char *argv[])
{
    LinearInterpolatedExternalMaterialLaw temperatureProfile(std::make_pair("x","temperature"), "temperature_profile.txt") ;
    LinearInterpolatedExternalMaterialLaw neutronProfile(std::make_pair("x","neutron_fluence"), "neutron_fluence_profile.txt") ;
    ThermalExpansionMaterialLaw thermalExpansion("temperature=293") ;
    RadiationInducedExpansionMaterialLaw aggregateExpansion("radiation_expansion_delay=0.15, neutron_fluence_correction = 0.57, maximum_radiation_expansion=0.0073") ;

    LogarithmicCreepWithExternalParameters paste("young_modulus = 12e9, poisson_ratio = 0.2, creep_modulus = 31.5e12, creep_poisson = 0.2, creep_characteristic_time = 4500, thermal_expansion_coefficient=8e-6") ;
    paste.addMaterialLaw(&temperatureProfile);
    paste.addMaterialLaw(&neutronProfile);
    paste.addMaterialLaw(&thermalExpansion);
    LogarithmicCreepWithExternalParameters aggregates("young_modulus = 70e9, poisson_ratio = 0.3, thermal_expansion_coefficient=12e-6") ;
    aggregates.addMaterialLaw(&temperatureProfile);
    aggregates.addMaterialLaw(&neutronProfile);
    aggregates.addMaterialLaw(&thermalExpansion);
    aggregates.addMaterialLaw(&aggregateExpansion);
    LogarithmicCreepWithExternalParameters concrete("young_modulus = 17e9, poisson_ratio = 0.2, creep_modulus = 31.5e12, creep_poisson = 0.2, creep_characteristic_time = 4500, thermal_expansion_coefficient=10e-6") ;
    concrete.addMaterialLaw(&temperatureProfile);
    concrete.addMaterialLaw(&neutronProfile);
    concrete.addMaterialLaw(&thermalExpansion);

    Sample box(nullptr, 1.,0.3,0.,0.) ;
    Sample left(nullptr, 0.3,0.3, -0.35, 0.) ;
    box.setBehaviour( &concrete ) ;
    left.setBehaviour( &paste ) ;

    FeatureTree F(&box) ;
    F.setSamplingNumber(32) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
    int time_step = 1 ;
    F.setDeltaTime(time_step) ;
    F.setMinDeltaTime(1e-9) ;
    F.setSamplingFactor(&box, 0.5);
    F.setSamplingFactor(&left, 5.);
    F.addRefinementZone(new Rectangle(0.3,0.3,-0.35,0.));


    F.addFeature(nullptr, &left) ;
    std::vector<Feature *> inclusions = PSDGenerator::get2DConcrete(&F, &aggregates, 100, 0.025, 0.00001,new PSDBolomeA(), CIRCLE, 1., M_PI, 100000, 0.8, new Rectangle(0.295,0.295,-0.35,0.) ) ;
    for(size_t i = 0 ; i < inclusions.size() ;i++)
        F.setSamplingFactor(inclusions[i], 3.);

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;

    F.step() ;

    std::cout << F.getCurrentTime() << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1.)[1] << std::endl ;
    TriangleWriter trg("tata", &F, 1.) ;
    trg.getField(STRAIN_FIELD) ;
    trg.getField(REAL_STRESS_FIELD) ;
    trg.write() ;

    return 0 ;
}

