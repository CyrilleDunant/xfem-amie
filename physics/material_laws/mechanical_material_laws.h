// C++ Interface: logarithmic_creep
//
// Description: Logarithmic visco-elastic behaviour for the Space-Time Finite Element Method
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2014
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __MECHANICAL_MATERIAL_LAW_H_
#define __MECHANICAL_MATERIAL_LAW_H_

#include "material_laws.h"

namespace Amie
{

/*PARSE BulkShearConversion ExternalMaterialLaw */
struct BulkShearConversionMaterialLaw : public ExternalMaterialLaw
{
	BulkShearConversionMaterialLaw( std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
	virtual ~BulkShearConversionMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE AdjustStrainStressCurve ExternalMaterialLaw */
struct AdjustStrainStressCurveMaterialLaw : public ExternalMaterialLaw
{
    AdjustStrainStressCurveMaterialLaw( std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
    virtual ~AdjustStrainStressCurveMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

/*PARSE BazantLoadNonLinearCreep ExternalMaterialLaw */
struct BazantLoadNonLinearCreepMaterialLaw : public ExternalMaterialLaw
{
    BazantLoadNonLinearCreepMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
    virtual ~BazantLoadNonLinearCreepMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

/*PARSE TensionCompressionCreep ExternalMaterialLaw */
struct TensionCompressionCreepMaterialLaw : public ExternalMaterialLaw
{
    TensionCompressionCreepMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
    virtual ~TensionCompressionCreepMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

struct StrainRateDependentStrengthMaterialLaw : public ExternalMaterialLaw
{
	double strainRateRef ;
	double p ;

	StrainRateDependentStrengthMaterialLaw(double p_ = 0.1, double eps = 1, std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep), strainRateRef(eps), p(p_) { } 
    virtual ~StrainRateDependentStrengthMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;


} 


#endif // __MATERIAL_LAW_H_
