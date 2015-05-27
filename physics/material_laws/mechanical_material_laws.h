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

struct BulkShearConversionExternalMaterialLaw : public ExternalMaterialLaw
{
	BulkShearConversionExternalMaterialLaw( std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
	virtual ~BulkShearConversionExternalMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
} ;

struct LoadNonLinearCreepMaterialLaw : public ExternalMaterialLaw
{
    LoadNonLinearCreepMaterialLaw(std::string args = std::string(), char sep = ',') : ExternalMaterialLaw(args, sep) { } 
    virtual ~LoadNonLinearCreepMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) ;
};

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
