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

int main(int argc, char *argv[])
{
    Function f = 1.-f_exp(Function("-0.2 t *")) ;
    Function g = 10. * f_sin(Function("t")) ;

    SpaceTimeDependentExternalMaterialLaw * temperature = new SpaceTimeDependentExternalMaterialLaw("temperature", f, true) ;
    SpaceTimeDependentExternalMaterialLaw * temperature2 = new SpaceTimeDependentExternalMaterialLaw("temperature", g, true) ;
//    LinearInterpolatedExternalMaterialLaw * temperature = new LinearInterpolatedExternalMaterialLaw(std::make_pair("t","temperature"),"temperature_history.txt") ;
    LinearInterpolatedExternalMaterialLaw * humidity = new LinearInterpolatedExternalMaterialLaw(std::make_pair("t","relative_humidity"),"humidity_history.txt") ;
    ThermalExpansionMaterialLaw * thermalExpansion = new ThermalExpansionMaterialLaw("temperature = 293") ;
    DryingShrinkageMaterialLaw * dryingShrinkage = new DryingShrinkageMaterialLaw() ;

    // dimensions
    Sample box(nullptr, 0.2,0.2,0.,0.) ;

    FeatureTree F(&box) ;
    // number of points on an edge
    F.setSamplingNumber(20) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
    F.setDeltaTime(1) ;

    // mechanical behaviour for cement paste
    LogarithmicCreepWithExternalParameters paste("young_modulus = 12e9, poisson_ratio = 0.2, thermal_expansion_coefficient = 10e-6, drying_shrinkage_coefficient = 0.0001") ;
    paste.addMaterialLaw(temperature);
    paste.addMaterialLaw(humidity);
    paste.addMaterialLaw(thermalExpansion);
    paste.addMaterialLaw(dryingShrinkage);
    box.setBehaviour( &paste );

    // mechanical behaviour for aggregates
    LogarithmicCreepWithExternalParameters aggregate("young_modulus = 70e9, poisson_ratio = 0.2, thermal_expansion_coefficient = 8e-6") ;
    aggregate.addMaterialLaw(temperature);
    aggregate.addMaterialLaw(thermalExpansion);
    int number_of_aggregates = 10 ;
    double aggregate_max_radius = 0.008 ;
    double aggregate_spacing = 0.0001 ;
    PSDGenerator::get2DConcrete( &F, &aggregate, number_of_aggregates, aggregate_max_radius, aggregate_spacing, new ConstantSizeDistribution()) ;


    // boundary conditions
//    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6)  );
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0.00001, 1)  );
//    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0., 3)  );
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0., 0)  );
    F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0., 2)  );

    // make calculation
    std::fstream file ;
    file.open("output3", std::ios::out) ;
    MultiTriangleWriter trg("first_calculation_header","first_calculation", &F, 1.) ;
    while(F.getCurrentTime() < 28)
    {

        F.step() ;
        Vector strain = F.getAverageField(STRAIN_FIELD, -1, 1.) ;
        Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1.) ;
        file << F.getCurrentTime() << "\t" ;
        file << strain[0] << "\t" ;
        file << strain[1] << "\t" ;
        file << strain[2] << "\t" ;
        file << stress[0] << "\t" ;
        file << stress[1] << "\t" ;
        file << stress[2] << std::endl ;

        trg.getField(STRAIN_FIELD) ;
        trg.getField(REAL_STRESS_FIELD) ;
        trg.getField(TWFT_STIFFNESS) ;
        trg.getField("temperature") ;
        trg.append() ;
        trg.writeSvg(); ;
    }

    // get stress and strain

    // write results



    return 0 ;
}

