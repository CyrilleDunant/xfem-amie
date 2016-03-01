/* this is an auto-generated file created on 29/1/2016 at 12:39  */

#ifndef __OBJECT_TRANSLATOR_H__
#define __OBJECT_TRANSLATOR_H__

#include "../physics/material_laws/material_laws.h"
#include "../elements/integrable_entity.h"
#include "../physics/material_laws/logcreep_accumulator.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../features/microstructuregenerator.h"
#include "../utilities/granulo.h"
#include "../utilities/inclusion_family.h"
#include "../features/features.h"
#include "../geometry/sampler/sampler.h"
#include "../utilities/postprocessor.h"

namespace Amie
{

struct Object
{

    // parsed from header file: ../physics/material_laws/material_laws.h
    static ExternalMaterialLaw * getExternalMaterialLaw(std::string type, std::map<std::string, std::string> & strings, std::map<std::string, std::vector<std::string>> & stringlists, std::map<std::string, double> & values) ;
    static bool isExternalMaterialLaw(std::string type) ;
    static void resetExternalMaterialLaw(ExternalMaterialLaw * target) ;

    // parsed from header file: ../elements/integrable_entity.h
    static Form * getForm(std::string type, std::map<std::string, double> & values, std::map<std::string, LogCreepAccumulator*> & logcreepaccumulators, std::map<std::string, std::string> & strings, std::map<std::string, FractureCriterion*> & fracturecriterions, std::map<std::string, DamageModel*> & damagemodels, std::map<std::string, ExternalMaterialLawList*> & externalmateriallawlists) ;
    static bool isForm(std::string type) ;
    static void resetForm(Form * target) ;

    // parsed from header file: ../physics/damagemodels/damagemodel.h
    static DamageModel * getDamageModel(std::string type, std::map<std::string, double> & values) ;
    static bool isDamageModel(std::string type) ;
    static void resetDamageModel(DamageModel * target) ;

    // parsed from header file: ../physics/fracturecriteria/fracturecriterion.h
    static FractureCriterion * getFractureCriterion(std::string type, std::map<std::string, double> & values, std::map<std::string, std::string> & strings) ;
    static bool isFractureCriterion(std::string type) ;
    static void resetFractureCriterion(FractureCriterion * target, std::map<std::string, double> & values, std::map<std::string, std::string> & strings) ;

    // parsed from header file: ../physics/material_laws/logcreep_accumulator.h
    static LogCreepAccumulator * getLogCreepAccumulator(std::string type, std::map<std::string, double> & values) ;
    static bool isLogCreepAccumulator(std::string type) ;
    static void resetLogCreepAccumulator(LogCreepAccumulator * target) ;

    // parsed from header file: ../features/microstructuregenerator.h
    static InclusionGenerator * getInclusionGenerator(std::string type, std::map<std::string, double> & values) ;
    static bool isInclusionGenerator(std::string type) ;
    static void resetInclusionGenerator(InclusionGenerator * target) ;

    // parsed from header file: ../utilities/granulo.h
    static ParticleSizeDistribution * getParticleSizeDistribution(std::string type, std::map<std::string, double> & values, std::map<std::string, std::string> & strings) ;
    static bool isParticleSizeDistribution(std::string type) ;
    static void resetParticleSizeDistribution(ParticleSizeDistribution * target) ;

    // parsed from header file: ../utilities/inclusion_family.h
    static InclusionFamily * getInclusionFamily(std::string type, std::map<std::string, double> & values, std::map<std::string, ParticleSizeDistribution*> & particlesizedistributions, std::map<std::string, InclusionGenerator*> & inclusiongenerators, std::map<std::string, std::string> & strings) ;
    static bool isInclusionFamily(std::string type) ;
    static void resetInclusionFamily(InclusionFamily * target) ;

    // parsed from header file: ../features/features.h
    static EnrichmentManager * getEnrichmentManager(std::string type, std::map<std::string, FeatureTree*> & featuretrees, std::map<std::string, InclusionFamily*> & inclusionfamilys, std::map<std::string, double> & values, std::map<std::string, std::string> & strings) ;
    static bool isEnrichmentManager(std::string type) ;
    static void resetEnrichmentManager(EnrichmentManager * target) ;

    // parsed from header file: ../geometry/sampler/sampler.h
    static Sampler * getSampler(std::string type, std::map<std::string, Point> & points, std::map<std::string, double> & values, std::map<std::string, std::string> & strings) ;
    static bool isSampler(std::string type) ;
    static void resetSampler(Sampler * target) ;

    // parsed from header file: ../utilities/postprocessor.h
    static PostProcessor * getPostProcessor(std::string type, std::map<std::string, std::string> & strings, std::map<std::string, double> & values) ;
    static bool isPostProcessor(std::string type) ;
    static void resetPostProcessor(PostProcessor * target) ;

} ;

}

#endif // __OBJECT_TRANSLATOR_H__
