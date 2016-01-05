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
#include "../../utilities/mineral.h"

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

/*PARSE Mineral ExternalMaterialLaw
    @string[filename] // path to the file containing the mineral data
    @string[separators] .- // separator to find the name of the mineral from the filename
    @value[index] -1 // index of the source of data (if several references exist)
    @value[factor] 1 // factor by which the stiffness matrix will be multiplied
    @string<bool>[force] false // ignores discrepancies between the mineral symmetry and its components
    @string<Variable>[cutting_plane] ZETA // cutting plane for the 3D to 2D projection
*/  
struct MineralMaterialLaw : public ExternalMaterialLaw
{
    Mineral mineral ;
    Variable cuttingPlane ;

    MineralMaterialLaw( Mineral mineral, Variable cut = ZETA, std::string args = std::string(), char sep = ',' ) ;
    MineralMaterialLaw( std::string filename, std::string sepstr = ".-", int index = -1, double factor = 1, bool force = false, Variable cut = ZETA, std::string args = std::string(), char sep = ',' ) ;
    virtual ~MineralMaterialLaw() { } ;

    virtual void preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt ) { }  ;
    virtual void preProcess( Matrix & stiffness, Point angle, planeType plane ) ;
} ;


} 


#endif // __MATERIAL_LAW_H_
