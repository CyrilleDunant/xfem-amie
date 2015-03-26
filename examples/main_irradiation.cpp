// Author: Alain Giorla <alain.b.giorla@gmail.com>, (C) 2005-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../polynomial/vm_function_base.h"
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

    ThermalExpansionMaterialLaw thermalExpansion("temperature=293") ;

    std::string file = "mortar_visco_mid3" ;
    std::fstream out;
    out.open(file.c_str(), std::ios::out) ;


    Function f_neutron_fluence("0.01 t *") ;
    Function f_drying_shrinkage("-0.00191 x * 0.00092  -") ;
    Function f_weight_loss("0.015 x 0.0075 * +") ;
    Function f_ageing_creep("y t x / 1 + *") ;
    SpaceTimeDependentExternalMaterialLaw neutronRadiation("neutron_fluence", f_neutron_fluence) ;
    RadiationInducedExpansionMaterialLaw radiationExpansion ;
    LinearInterpolatedExternalMaterialLaw radiationDamage( std::make_pair("neutron_fluence","young_modulus"), "../examples/data/irradiation/limestone_aggregate_modulus" ) ;

    SimpleDependentExternalMaterialLaw dryingShrinkage("imposed_deformation","weight_loss", f_drying_shrinkage, ADD) ;
    SimpleDependentExternalMaterialLaw weightLoss("weight_loss", "neutron_fluence", f_weight_loss ) ;

    VariableDependentExternalMaterialLaw ageingCreep("creep_modulus",  f_ageing_creep) ;
    ageingCreep.setAsX("creep_characteristic_time") ;
    ageingCreep.setAsY("creep_modulus_28days") ;

    LogarithmicCreepWithExternalParameters limestone("young_modulus=72.4e9, poisson_ratio=0.25, imposed_deformation=0, thermal_expansion_coefficient=6.35e-6, temperature=318, maximum_radiation_expansion = 0.0187, radiation_expansion_delay=0.0008, neutron_fluence_correction=3.0411") ;
    limestone.addMaterialLaw(&thermalExpansion) ;
    limestone.addMaterialLaw(&neutronRadiation);
    limestone.addMaterialLaw(&radiationExpansion);
//    limestone.addMaterialLaw(&radiationDamage);
    LogarithmicCreepWithExternalParameters paste("young_modulus=22e9, poisson_ratio=0.2, imposed_deformation=0, thermal_expansion_coefficient=10e-6, temperature=318, creep_modulus=4.5e9, creep_characteristic_time=0.25, creep_poisson = 0.2, creep_modulus_28days=4.5e9") ;
//    paste.setLogCreepAccumulator(LOGCREEP_CONSTANT);
    paste.addMaterialLaw( &ageingCreep ) ;
    paste.addMaterialLaw(&thermalExpansion) ;
    paste.addMaterialLaw(&neutronRadiation);
    paste.addMaterialLaw(&weightLoss);
    paste.addMaterialLaw(&dryingShrinkage);

    Sample box(nullptr, 0.0635,0.0127,0.0635/2,0.) ;
    box.setBehaviour(&paste);

    FeatureTree F(&box) ;
    F.setSamplingNumber(32) ; // calculation with 96
    F.setSamplingRestriction(SAMPLE_RESTRICT_4);
    F.setOrder(LINEAR_TIME_LINEAR) ;
    double time_step = 0.1 ;
    F.setDeltaTime(time_step) ;
    F.setMinDeltaTime(1e-9) ;

    std::vector<Feature *> inclusions = PSDGenerator::get2DConcrete(&F, &limestone, 10000, 0.00625, 0.00001, new GranuloFromCumulativePSD("../examples/data/irradiation/limestone_aggregate_psd", CUMULATIVE_PERCENT), CIRCLE, 1., M_PI, 1000000, 0.6, new Rectangle(0.0635, 0.279,0.0635/2,0.) ) ;
    std::vector<Geometry *> placed ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
    {
        if(box.intersects(inclusions[i]) || box.in(inclusions[i]->getCenter()))
        {
            F.setSamplingFactor(inclusions[i], 3.) ;
            placed.push_back(new Circle(inclusions[i]->getRadius(), inclusions[i]->getCenter())) ;
        }
    }
    int count = placed.size() ;
    std::vector<Feature *> small = PSDGenerator::get2DConcrete(&F, &limestone, 1000, 0.0001, 0.00003, new ConstantSizeDistribution(), CIRCLE, 1., M_PI, 500000, 0.6, new Rectangle(0.0635, 0.0127,0.0635/2,0.), placed ) ;
    count += small.size() ;
    std::cout << count << " aggregates in slice" << std::endl ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0, 1)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0, 3)) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;

//    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT_AFTER, 1e6));

    double aggregateArea = 0. ;
    std::vector<bool> iAgg ;
    std::vector<bool> iCem ;
    int done = 1 ;

    while(VirtualMachine().eval(f_neutron_fluence,0.0635/2.,0.,0.,F.getCurrentTime()) < 10.)
//    while(F.getCurrentTime() < 400)
    {
//        F.setDeltaTime(0.01*done);
        F.step() ;
        Vector strain = F.getAverageField(STRAIN_FIELD, -1, 1.) ;
        Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1.) ;

        Vector strainAgg(3) ;
        Vector strainCem(3) ;
        Vector stressAgg(3) ;
        Vector stressCem(3) ;
        Vector strainMaxAgg(3) ;
        Vector strainMaxCem(3) ;
        Vector stressMaxAgg(3) ;
        Vector stressMaxCem(3) ;
        Vector tmp(3) ;
        done++ ;

        if(aggregateArea < POINT_TOLERANCE)
        {
            for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
            {
                if(i->getBehaviour()->param[0][0] > 30e9)
                {
                    aggregateArea += i->area() ;
                    iAgg.push_back(true) ;
                    iCem.push_back(false);
                } else {
                    iCem.push_back(true) ;
                    iAgg.push_back(false);
                }
            }
            std::cout << aggregateArea / (0.0635*0.0127) << std::endl ;
        }

        for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
        {
            if(iAgg[i.getPosition()])
            {
                double a = i->area() ;
                i->getState().getAverageField(STRAIN_FIELD, tmp, nullptr, -1, 1.) ;
                strainAgg += tmp*a ;
                for(size_t j = 0 ; j < 3 ; j++)
                {
                    if(tmp[j] > strainMaxAgg[j])
                        strainMaxAgg[j] = tmp[j] ;
                }
                i->getState().getAverageField(REAL_STRESS_FIELD, tmp, nullptr, -1, 1.) ;
                stressAgg += tmp*a ;
                for(size_t j = 0 ; j < 3 ; j++)
                {
                    if(tmp[j] > stressMaxAgg[0])
                        stressMaxAgg[0] = tmp[0] ;
                }
            }
        }

        strainAgg /= aggregateArea ;
        stressAgg /= aggregateArea ;


        for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
        {
            if(iCem[i.getPosition()])
            {
                double a = i->area() ;
                i->getState().getAverageField(STRAIN_FIELD, tmp, nullptr, -1, 1.) ;
                strainAgg += tmp*a ;
                for(size_t j = 0 ; j < 3 ; j++)
                {
                    if(tmp[j] > strainMaxAgg[j])
                        strainMaxAgg[j] = tmp[j] ;
                }
                i->getState().getAverageField(REAL_STRESS_FIELD, tmp, nullptr, -1, 1.) ;
                stressAgg += tmp*a ;
                for(size_t j = 0 ; j < 3 ; j++)
                {
                    if(tmp[j] > stressMaxAgg[0])
                        stressMaxAgg[0] = tmp[0] ;
                }
            }
        }
        strainCem /= box.area() - aggregateArea ;
        stressCem /= box.area() - aggregateArea ;

        std::cout << F.getCurrentTime() << "\t" << VirtualMachine().eval(f_neutron_fluence,0.0635/2.,0.,0.,F.getCurrentTime()) << "\t" << strain[0] << "\t" << strain[1] << "\t" << stress[0] << "\t" << stress[1] <<
                     "\t" << strainAgg[0] << "\t" << strainAgg[1] << "\t" << stressAgg[0] << "\t" << stressAgg[1] <<
                     "\t" << strainMaxAgg[0] << "\t" << strainMaxAgg[1] << "\t" << stressMaxAgg[0] << "\t" << stressMaxAgg[1] <<
                     "\t" << strainCem[0] << "\t" << strainCem[1] << "\t" << stressCem[0] << "\t" << stressCem[1] <<
                     "\t" << strainMaxCem[0] << "\t" << strainMaxCem[1] << "\t" << stressMaxCem[0] << "\t" << stressMaxCem[1] << std::endl ;

        out << F.getCurrentTime() << "\t" << VirtualMachine().eval(f_neutron_fluence,0.0635/2.,0.,0.,F.getCurrentTime()) << "\t" << strain[0] << "\t" << strain[1] << "\t" << stress[0] << "\t" << stress[1] <<
                     "\t" << strainAgg[0] << "\t" << strainAgg[1] << "\t" << stressAgg[0] << "\t" << stressAgg[1] <<
                     "\t" << strainMaxAgg[0] << "\t" << strainMaxAgg[1] << "\t" << stressMaxAgg[0] << "\t" << stressMaxAgg[1] <<
                     "\t" << strainCem[0] << "\t" << strainCem[1] << "\t" << stressCem[0] << "\t" << stressCem[1] <<
                     "\t" << strainMaxCem[0] << "\t" << strainMaxCem[1] << "\t" << stressMaxCem[0] << "\t" << stressMaxCem[1] << std::endl ;

    }
/*    TriangleWriter writer("tata", &F, 1.) ;
    writer.getField("young_modulus") ;
    writer.getField("creep_modulus") ;
    writer.write() ;*/

    return 0 ;
}

