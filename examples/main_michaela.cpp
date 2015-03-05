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
#include "../physics/material_laws/humidity_material_laws.h"
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

double tMax=10;

struct ElleuchRelativeHumidityHistoryExternalMaterialLaw : public ExternalMaterialLaw
{
    ElleuchRelativeHumidityHistoryExternalMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
    {
        double t = s.getNodalCentralTime() ;
        double k = s.get("permeability_coefficient", defaultValues) ;
        double h = 1. ;
        if(t < 90)
            h = 0.6 +  0.4*exp(-t/k) ;
        else
        {
            h = 0.6 + 0.4*exp(-90/k) ;

            if(t < 116)
                h -= 0.4*(t-90)/(116-90) ;
            else
                h -= 0.4 ;

            if(t > 148)
                h -= 0.2*(1-exp(-(t-148)/k)) ;
        }

        s.set("relative_humidity", h) ;
    }
};

struct ElleuchHydratingModulusExternalLaw : public ExternalMaterialLaw
{
    ElleuchHydratingModulusExternalLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { }

    virtual void preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
    {
        double t = s.getNodalCentralTime() ;
        double E0 = s.get("young_modulus_90d", defaultValues) ;
        s.set("young_modulus", E0*(0.2+0.8*exp(-t/28))) ;
    }

};

int main(int argc, char *argv[])
{
    std::string rig = "g1" ;//(argv[1]) ; // g1 g2 g3 g4 or g5
//material behaviour
    bool simple = true ;//(std::string(argv[2]) == std::string("simple")) ; // full
    bool visco = false ;//(std::string(argv[3]) == std::string("visco")) ; // elastic
    bool plastic = false ;//(std::string(argv[4]) == std::string("plastic")) ; // nodamage
    bool brittle = false ;//(std::string(argv[4]) == std::string("brittle")) ;

//environment history functions
   // std::string temperatureHistory = "../examples/data/elleuch/temperature_history_"+rig ;
    //std::string fluenceHistory = "../examples/data/elleuch/fluence_history_"+rig ;
    std::string outputFile = "testing/output_e" ;

//    for(int i = 1 ; i < argc ; i++)
  //  {
    //    outputFile.append("_") ;
      //  outputFile.append(argv[i]) ;
    //}
    std::string trgFile = outputFile+"_trg" ;

    std::fstream out ;
    out.open(outputFile.c_str(), std::ios::out) ;

    LinearInterpolatedExternalMaterialLaw * temperature = new LinearInterpolatedExternalMaterialLaw(std::make_pair("t","temperature"),"temperature_history_misa.txt") ;
     LinearInterpolatedExternalMaterialLaw * humidity = new LinearInterpolatedExternalMaterialLaw(std::make_pair("t","relative_humidity"),"RH_misa.txt") ;
      LinearInterpolatedExternalMaterialLaw * fluence = new LinearInterpolatedExternalMaterialLaw(std::make_pair("t","neutron_fluence"),"fluence_misa.txt") ;

//applying of the function to material law
   // LinearInterpolatedExternalMaterialLaw temperature(std::make_pair("t", "temperature"), temperatureHistory) ;
   // LinearInterpolatedExternalMaterialLaw fluence(std::make_pair("t", "neutron_fluence"), fluenceHistory) ;
   // ElleuchRelativeHumidityHistoryExternalMaterialLaw humidity("permeability_coefficient = 25") ;

    ThermalExpansionMaterialLaw thermalExpansion("temperature=293") ;
    DryingShrinkageMaterialLaw dryingShrinkage("relative_humidity=0.6") ;
    RadiationInducedExpansionMaterialLaw radiationExpansion ;
    LinearInterpolatedExternalMaterialLaw radiationDamage(std::make_pair("neutron_fluence","young_modulus"), "../examples/data/elleuch/aggregate_damage") ;
    CreepArrheniusMaterialLaw creepArrhenius("temperature=293, creep_modulus = 3.6e9, creep_characteristic_time=0.5, creep_activation_energy = 1.666e-4") ;
    CreepRelativeHumidityMaterialLaw creepHumidity("creep_humidity_coefficient = 5") ;

//damage behaviour
    FractureCriterion * pasteCriterion = nullptr ;
    DamageModel * damageModel = nullptr ;
    FractureCriterion * serpentineCriterion = nullptr ;
    if(plastic || brittle)
    {
        if(plastic)
        {
            pasteCriterion = new SpaceTimeNonLocalMaximumStress(5e6) ;
            serpentineCriterion = new SpaceTimeNonLocalMaximumStress(50e6) ;
            damageModel = new SpaceTimeFiberBasedIsotropicLinearDamage(0.01, 1e-6, 0.9) ;

        }
        else
        {
            pasteCriterion = new SpaceTimeNonLocalMaximumStress(17e6) ;
            serpentineCriterion = new SpaceTimeNonLocalMaximumStress(50e6) ;
            damageModel = new SpaceTimeFiberBasedIsotropicLinearDamage(0.01, 1e-6, 0.9) ;//cisla?
        }
        //co to je?velikost trhliny?
        pasteCriterion->setMaterialCharacteristicRadius(0.0001) ;
        serpentineCriterion->setMaterialCharacteristicRadius(0.0001) ;
    }

//cement paste properties, if not given visco-prop=>just elastic
    LogarithmicCreepWithExternalParameters paste("young_modulus = 20.1e9, young_modulus_90d = 20.1e9, poisson_ratio = 0.2, thermal_expansion_coefficient=9e-6, drying_shrinkage_coefficient=0.0002, imposed_deformation=0", pasteCriterion, damageModel ) ; //, creep_modulus = 3.6e9, creep_characteristic_time=0.5, creep_poisson = 0.2
    if(visco)
    {
        paste.addMaterialParameter("creep_modulus", 3.6e9);
        paste.addMaterialParameter("creep_characteristic_time", 0.5);
        paste.addMaterialParameter("creep_poisson", 0.2);
    }
    paste.addMaterialLaw(temperature);
    paste.addMaterialLaw(fluence);
    paste.addMaterialLaw(humidity);
    paste.addMaterialLaw(&thermalExpansion);
    paste.addMaterialLaw(&dryingShrinkage);
    if(visco)
    {
        paste.addMaterialLaw(&creepArrhenius);
        paste.addMaterialLaw(&creepHumidity);
    }

    LogarithmicCreepWithExternalParameters serpentine("young_modulus = 79e9, poisson_ratio = 0.2, thermal_expansion_coefficient=7.5e-6, radiation_expansion_delay=84.268, maximum_radiation_expansion=0.0972, neutron_fluence_correction = 0.0099, imposed_deformation=0", serpentineCriterion, damageModel) ;
    serpentine.addMaterialLaw(temperature);
    serpentine.addMaterialLaw(fluence);
    serpentine.addMaterialLaw(humidity);
    serpentine.addMaterialLaw(&thermalExpansion);
    serpentine.addMaterialLaw(&radiationExpansion);
    serpentine.addMaterialLaw(&radiationDamage);

//creation of sample and its dimensions
    Sample box(nullptr, 0.02,0.02,0.,0.) ;//0.025
    box.setBehaviour(&paste) ;

    int sampling = 1 ;//25
    if(simple)
        sampling = 1 ;

    FeatureTree F(&box) ;
    F.setSamplingNumber(sampling) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
   // F.setSamplingRestriction(SAMPLE_RESTRICT_8);//stop the cal if not enough elements on boundary
    double time_step = 0.5 ;
    F.setMaxIterationsPerStep(100) ;
    F.setDeltaTime(time_step) ;
    F.setMinDeltaTime(1e-6) ;
    F.setSamplingFactor(&box, 0.25);//number of elements multiplied by x..(&box, x)

    int nagg = 5000 ;
    if(simple)
        nagg = 5;//number of aggregate,5000

    std::vector<Feature *> allInclusions = PSDGenerator::get2DConcrete(&F, &serpentine, nagg, 0.005, 0.00002, new GranuloFromCumulativePSD("aggregate_psd_m_el", CUMULATIVE_PERCENT), CIRCLE, 1., M_PI, 100000, 0.5, new Rectangle(0.1,0.1,0.,0.)) ;//calculation, aggregate, number of agg, max grain, min spacing,...curve of particle distribution...shape, parameters for elips, number of tries to create the aggregate box, percentage of agg area in sample, dimensions of sample=>the sample for model will be cut from this one

    std::vector<Circle *> inInclusions ;
    for(size_t i = 0 ; i < allInclusions.size() ; i++)
    {
        if(box.intersects(allInclusions[i]) || box.in(allInclusions[i]->getCenter()))
            inInclusions.push_back(new Circle(allInclusions[i]->getRadius(), allInclusions[i]->getCenter())) ;
    }

//boundary conditions of box
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;
    //just for the stress load
    BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, 0) ;
    F.addBoundaryCondition(load);


    F.step() ;
    F.getAssembly()->setEpsilon(1e-20);//precision of calculation..square from epsilon
//deviding of elements to paste and agg ones

        int pasteCachePosition = F.get2DMesh()->generateCache( F.getFeature(0) ) ;
std::vector<Geometry *> allinc ;
for(size_t i = 0 ; i < allInclusions.size() ; i++)
    allinc.push_back( dynamic_cast<Geometry*>(allInclusions[i])) ;
int aggCachePosition = F.get2DMesh()->generateCache( allinc ) ;

    double areaPaste = F.get2DMesh()->getArea( pasteCachePosition ) ;
    double areaAgg = F.get2DMesh()->getArea( aggCachePosition ) ;
    std::cout << "area " << pasteCachePosition  << "   " << aggCachePosition << std::endl ;
	exit(0) ;



    std::vector<DelaunayTriangle *> elements = F.get2DMesh()->getConflictingElements(dynamic_cast<Rectangle *>(F.getFeature(0))) ;
    std::vector<size_t> iAggregates ;
    std::vector<size_t> iPaste ;


//forgot
    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * state = dynamic_cast< GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables * >(elements[0]->getStatePointer()) ;
    std::map<std::string, double> dummy ;
//if i want to apply stress
//    BoundingBoxDefinedBoundaryCondition * stress = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, 0., 1) ;
//    F.addBoundaryCondition(stress) ;


//calculation

    while(F.getCurrentTime() < 1)
    {
        double t = F.getCurrentTime();
        double stress_load;
       if(t >= tMax - time_step)

          stress_load = -15e6;

       else
       {
           stress_load = 0;
       }



//        stress->setData(0.00001*F.getCurrentTime()) ;
        F.setDeltaTime(time_step);
        load->setData(stress_load);//load in calculation
        bool goOn = true ;
        do {
            goOn = F.step() ;

           // TriangleWriter writer("tata", &F, 1.) ;
           // writer.getField(STRAIN_FIELD) ;
           // writer.getField(REAL_STRESS_FIELD) ;
           // writer.getField(TWFT_STIFFNESS) ;
           // writer.getField(TWFT_DAMAGE) ;
           // writer.getField(TWFT_VISCOSITY) ;
           // writer.write();
//end of calculation if not converged
            if(!F.solverConverged())
            {
                goOn = true ;
                std::cout << "STOP" ;
                break ;
            }

        } while(!goOn) ;

        if(!F.solverConverged())
            break ;
//saving to file
        Vector strain = F.getAverageField(STRAIN_FIELD, -1, 1.) ;
        Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1.) ;
        Vector strainPaste = F.get2DMesh()->getField(STRAIN_FIELD, pasteCachePosition, -1, 1.) ;
        Vector strainAgg = F.get2DMesh()->getField(STRAIN_FIELD, aggCachePosition, -1, 1.) ;

        double damagePaste = 0. ;
        double damageAggregates = 0. ;
        if(brittle || plastic)
        {
            for(size_t i = 0 ; i < iPaste.size() ; i++)
                damagePaste += elements[iPaste[i]]->area()*elements[iPaste[i]]->getBehaviour()->getDamageModel()->getState().max() ;
            for(size_t i = 0 ; i < iAggregates.size() ; i++)
                damageAggregates += elements[iAggregates[i]]->area()*elements[iAggregates[i]]->getBehaviour()->getDamageModel()->getState().max() ;
        }
        damagePaste /= areaPaste ;
        damageAggregates /= areaAgg ;

        std::cout << F.getCurrentTime() << "\t" << state->get("temperature", dummy) << "\t" << state->get("neutron_fluence", dummy) << "\t" << state->get("relative_humidity", dummy) <<  "\t"  << state->get("stress_load", dummy) <<  "\t" << strain[0] << "\t" << strain[1] << "\t" << stress[0] << "\t" << stress[1]  << "\t" << damageAggregates << "\t" << damagePaste << std::endl ;
        out << F.getCurrentTime() << "\t" << state->get("temperature", dummy) << "\t" << state->get("neutron_fluence", dummy) << "\t" << state->get("relative_humidity", dummy) << "\t" << state->get("stress_load", dummy) <<  "\t" << strain[0] << "\t" << strain[1] << "\t" << stress[0] << "\t" << stress[1] << "\t" << strainPaste[0] << "\t" << strainPaste[1] << "\t" << strainAgg[0] << "\t" << strainAgg[1]<< "\t" << damageAggregates << "\t" << damagePaste << std::endl ;

    }
 MultiTriangleWriter trg("elleuch_misa","elleuch_misa_2", &F, 1.) ;

    trg.getField(STRAIN_FIELD) ;
    trg.getField(REAL_STRESS_FIELD) ;
    trg.getField(TWFT_STIFFNESS) ;
    trg.getField(TWFT_DAMAGE) ;
    trg.append() ;
    trg.writeSvg(); ;
    return 0 ;
}

