// Author: Alain Giorla <alain.b.giorla@gmail.com>, (C) 2005-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/material_laws/humidity_material_laws.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/features.h"
#include "../utilities/granulo.h"
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

double tMax(std::string r)
{
    if(r == std::string("g1"))
        return 235 ;
    if(r == std::string("g2"))
        return 270 ;
    if(r == std::string("g3"))
        return 310 ;
    if(r == std::string("g4"))
        return 350 ;
    if(r == std::string("g5"))
        return 480 ;
    return 0. ;
}

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
    std::string rig(argv[1]) ; // g1 g2 g3 g4 or g5

    bool simple = (std::string(argv[2]) == std::string("simple")) ; // full
    bool visco = (std::string(argv[3]) == std::string("visco")) ; // elastic
    bool plastic = (std::string(argv[4]) == std::string("plastic")) ; // nodamage
    bool brittle = (std::string(argv[4]) == std::string("brittle")) ;

    std::string temperatureHistory = "../examples/data/elleuch/temperature_history_"+rig ;
    std::string fluenceHistory = "../examples/data/elleuch/fluence_history_"+rig ;
    std::string outputFile = "../examples/data/elleuch/output" ;
    for(int i = 1 ; i < argc ; i++)
    {
        outputFile.append("_") ;
        outputFile.append(argv[i]) ;
    }
    std::string trgFile = outputFile+"_trg" ;

    std::fstream out ;
    out.open(outputFile.c_str(), std::ios::out) ;

    LinearInterpolatedExternalMaterialLaw temperature(std::make_pair("t", "temperature"), temperatureHistory) ;
    LinearInterpolatedExternalMaterialLaw fluence(std::make_pair("t", "neutron_fluence"), fluenceHistory) ;
    ElleuchRelativeHumidityHistoryExternalMaterialLaw humidity("permeability_coefficient = 25") ;

    ThermalExpansionMaterialLaw thermalExpansion("temperature=293") ;
    DryingShrinkageMaterialLaw dryingShrinkage("relative_humidity=1.") ;
    RadiationInducedExpansionMaterialLaw radiationExpansion ;
    LinearInterpolatedExternalMaterialLaw radiationDamage(std::make_pair("neutron_fluence","young_modulus"), "data_elleuch/aggregate_damage") ;
    CreepArrheniusMaterialLaw creepArrhenius("temperature=293, creep_modulus = 3.6e9, creep_characteristic_time=0.5, creep_activation_energy = 1.666e-4") ;
    CreepRelativeHumidityMaterialLaw creepHumidity("creep_humidity_coefficient = 5") ;

    FractureCriterion * pasteCriterion = nullptr ;
    DamageModel * damageModel = nullptr ;
    FractureCriterion * serpentineCriterion = nullptr ;
    if(plastic || brittle)
    {
        if(plastic)
        {
            pasteCriterion = new SpaceTimeNonLocalMaximumStress(17e6) ;
            serpentineCriterion = new SpaceTimeNonLocalMaximumStress(50e6) ;
            damageModel = new SpaceTimeFiberBasedIsotropicLinearDamage(0.01, 1e-6, 0.9) ;

        }
        else
        {
            pasteCriterion = new SpaceTimeNonLocalMaximumStrain(0.00038) ;
            serpentineCriterion = new SpaceTimeNonLocalMaximumStrain(0.00039) ;
            damageModel = new SpaceTimeFiberBasedIsotropicLinearDamage(0.01, 1e-6, 0.9) ;
        }
        pasteCriterion->setMaterialCharacteristicRadius(0.0001) ;
        serpentineCriterion->setMaterialCharacteristicRadius(0.0001) ;
    }

    LogarithmicCreepWithExternalParameters paste("young_modulus = 43e9, young_modulus_90d = 43e9, poisson_ratio = 0.2, thermal_expansion_coefficient=9e-6, drying_shrinkage_coefficient=0.00039, imposed_deformation=0", pasteCriterion, damageModel ) ; //, creep_modulus = 3.6e9, creep_characteristic_time=0.5, creep_poisson = 0.2
    if(visco)
    {
        paste.addMaterialParameter("creep_modulus", 3.6e9);
        paste.addMaterialParameter("creep_characteristic_time", 0.5);
        paste.addMaterialParameter("creep_poisson", 0.2);
    }
    paste.addMaterialLaw(&temperature);
    paste.addMaterialLaw(&fluence);
    paste.addMaterialLaw(&humidity);
    paste.addMaterialLaw(&thermalExpansion);
    paste.addMaterialLaw(&dryingShrinkage);
    if(visco)
    {
        paste.addMaterialLaw(&creepArrhenius);
        paste.addMaterialLaw(&creepHumidity);
    }

    LogarithmicCreepWithExternalParameters serpentine("young_modulus = 130e9, poisson_ratio = 0.2, thermal_expansion_coefficient=7.5e-6, radiation_expansion_delay=84.268, maximum_radiation_expansion=0.0972, neutron_fluence_correction = 0.0099, imposed_deformation=0", serpentineCriterion, damageModel) ;
    serpentine.addMaterialLaw(&temperature);
    serpentine.addMaterialLaw(&fluence);
    serpentine.addMaterialLaw(&humidity);
    serpentine.addMaterialLaw(&thermalExpansion);
    serpentine.addMaterialLaw(&radiationExpansion);
    serpentine.addMaterialLaw(&radiationDamage);

    Sample box(nullptr, 0.025,0.02,0.,0.) ;
    box.setBehaviour(&paste) ;

    int sampling = 256 ;
    if(simple)
        sampling = 64 ;

    FeatureTree F(&box) ;
    F.setSamplingNumber(sampling) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
    F.setSamplingRestriction(SAMPLE_RESTRICT_4);
    double time_step = 0.5 ;
    F.setMaxIterationsPerStep(100) ;
    F.setDeltaTime(time_step) ;
    F.setMinDeltaTime(1e-6) ;

    int nagg = 5000 ;
    if(simple)
        nagg = 500 ;

    std::vector<Feature *> allInclusions = PSDGenerator::get2DConcrete(&F, &serpentine, nagg, 0.006, 0.00002, new GranuloFromCumulativePSD("data_elleuch/aggregate_psd", CUMULATIVE_PERCENT), nullptr, 100000, 0.67, new Rectangle(0.1,0.1,0.,0.)) ;
    std::vector<Circle *> inInclusions ;
    for(size_t i = 0 ; i < allInclusions.size() ; i++)
    {
        if(box.intersects(allInclusions[i]) || box.in(allInclusions[i]->getCenter()))
            inInclusions.push_back(new Circle(allInclusions[i]->getRadius(), allInclusions[i]->getCenter())) ;
    }

    double tmax = tMax(rig) ;

    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;

    F.step() ;
    F.getAssembly()->setEpsilon(1e-20);

    double areaAggregates = 0. ;
    double areaPaste = 0. ;
    for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
    {
        bool inPaste = true ;
        for(size_t j = 0 ; j < inInclusions.size() ; j++)
        {
            if(inInclusions[j]->in(i->getCenter()))
            {
                areaAggregates += i->area() ;
                inPaste = false ;
            }
        }
        if(inPaste)
        {
            areaPaste += i->area() ;
        }
    }


    GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & state = dynamic_cast< GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & >(F.get2DMesh()->begin()->getState()) ;
    std::map<std::string, double> dummy ;

//    BoundingBoxDefinedBoundaryCondition * stress = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, 0., 1) ;
//    F.addBoundaryCondition(stress) ;

    while(F.getCurrentTime() < tmax)
    {
//        stress->setData(0.00001*F.getCurrentTime()) ;
        F.setDeltaTime(time_step);
        bool goOn = true ;
        do {
            goOn = F.step() ;
            TriangleWriter writer("trg_elleuch", &F, 1.) ;
            writer.getField(PRINCIPAL_STRAIN_FIELD) ;
            writer.getField(PRINCIPAL_REAL_STRESS_FIELD) ;
            writer.getField(TWFT_STIFFNESS) ;
            writer.getField(TWFT_DAMAGE) ;
            writer.getField(TWFT_VISCOSITY) ;
            writer.write();
		exit(0) ;

            if(!F.solverConverged())
            {
                goOn = true ;
                std::cout << "STOP" ;
                break ;
            }

        } while(!goOn) ;

        if(!F.solverConverged())
            break ;

        Vector strain = F.getAverageField(STRAIN_FIELD, -1, 1.) ;
        Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1.) ;

        double damagePaste = 0. ;
        double damageAggregates = 0. ;
        if(brittle || plastic)
        {
            for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
            {
                bool inPaste = true ;
                for(size_t j = 0 ; j < inInclusions.size() ; j++)
                {
                    if(inInclusions[j]->in(i->getCenter()))
                    {
                        damagePaste += i->area()*i->getBehaviour()->getDamageModel()->getState().max() ;
                    }
                }
                if(inPaste)
                {
                    damageAggregates += i->area()*i->getBehaviour()->getDamageModel()->getState().max() ;
                }
            }
        }
        damagePaste /= areaPaste ;
        damageAggregates /= areaAggregates ;

        std::cout << F.getCurrentTime() << "\t" << state.get("temperature", dummy) << "\t" << state.get("neutron_fluence", dummy) << "\t" << state.get("relative_humidity", dummy) <<  "\t" << strain[0] << "\t" << strain[1] << "\t" << stress[0] << "\t" << stress[1]  << "\t" << damageAggregates << "\t" << damagePaste << std::endl ;
        out << F.getCurrentTime() << "\t" << state.get("temperature", dummy) << "\t" << state.get("neutron_fluence", dummy) << "\t" << state.get("relative_humidity", dummy) <<  "\t" << strain[0] << "\t" << strain[1] << "\t" << stress[0] << "\t" << stress[1] << "\t" << damageAggregates << "\t" << damagePaste << std::endl ;

    }

    TriangleWriter writer(trgFile, &F, 1.) ;
    writer.getField(TWFT_PRINCIPAL_STRAIN) ;
    writer.getField(TWFT_PRINCIPAL_STRESS) ;
    writer.getField(TWFT_STIFFNESS) ;
    writer.getField(TWFT_DAMAGE) ;
    writer.write();

    return 0 ;
}

