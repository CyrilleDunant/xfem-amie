// C++ Interface: generalized_spacetime_viscoelasticity
//
// Description: Generalized visco-elastic behaviour for the Space-Time Finite Element Method
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __VISCOELASTICTY_AND_IMPOSED_DEFORMATION_H_
#define __VISCOELASTICTY_AND_IMPOSED_DEFORMATION_H_

#include "viscoelasticity.h"

namespace Amie
{

/*PARSE ViscoelasticityAndImposedDeformation Form 
    @string<ViscoelasticModel>[model] // rheological assembly
    @value[young_modulus] // value of the Young modulus
    @value[poisson_ratio] // value of the Poisson ratio
    @value[imposed_deformation] // value of the elastic linear imposed deformation
    @value[creep_characteristic_time] 1 // value of the viscosity characteristic time
    @string[file_name] // name of the file containing the modulus of the different branches
    @string<SpaceDimensionality>[dimension] SPACE_TWO_DIMENSIONAL // number of dimensions of the current simulation
    @string<planeType>[plane_type] PLANE_STRESS // 2D hypothesis (plane strain or plane stress)
    @string<IsotropicMaterialParameters>[material_parameters] YOUNG_POISSON // describes how to build the stiffness matrix
 */
struct ViscoelasticityAndImposedDeformation : public Viscoelasticity
{
    Vector imposedStrain ;
    Vector imposedGeneralizedStrain ;
    Vector imposedStress ;

    ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & rig, int additionnalBlocksAfter = 0, double r = 0) ;
    // constructor for pure elasticity or pure viscosity
    ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & rig, Vector & imp, int additionnalBlocksAfter = 0, double r = 0) ;
    // constructor for elementary Kelvin-Voigt or Maxwell
    ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & rig, const Matrix & eta,  Vector & imp, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
    // constructor for Burger (KV + Maxwell in serial)
    ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & c_kv, const Matrix & e_kv, const Matrix & c_mx, const Matrix & e_mx,  Vector & imp, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
    // constructor for generalized KelvinVoigt or Maxwell
    ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & c_0, std::vector<std::pair<Matrix, Matrix> > & branches,  Vector & imp, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
    // constructor for generalized KelvinVoigt or Maxwell with 1 module only
    ViscoelasticityAndImposedDeformation( ViscoelasticModel model, const Matrix & c_0, const Matrix & c_1, const Matrix & e_1, Vector & imp,  int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0, double r = 0) ;
    // constructor for general viscoelasticity (rig and eta are supposed symmetric)
    ViscoelasticityAndImposedDeformation( const Matrix & rig, const Matrix & eta, int blocks, Vector & imp, int additionnalBlocksAfter = 0, double r = 0) ;

    ViscoelasticityAndImposedDeformation( ViscoelasticModel model, double young, double poisson, double alpha, double tau = 1, std::string file = std::string(), SpaceDimensionality dim = SPACE_TWO_DIMENSIONAL, planeType pt = PLANE_STRESS, IsotropicMaterialParameters hooke = YOUNG_POISSON, int additionnalBlocksBefore = 0, int additionnalBlocksAfter = 0) ;

    virtual ~ViscoelasticityAndImposedDeformation() ;

/*    virtual bool isViscous() const {
        return true ;
    }*/

    virtual Form * getCopy() const ;

    virtual void scale( double d )
    {
        param *= d ;
        eta *= d ;
        imposedStrain *= d ;
        imposedGeneralizedStrain *= d ;
        imposedStress *= d ;
    }

    virtual bool hasInducedForces() const {
        return true ;
    }

    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

protected:
    void makeImposedStress() ;

} ;

struct ViscoelasticityAndVariableImposedDeformation : public ViscoelasticityAndImposedDeformation
{
    Function variableImposed ;

    ViscoelasticityAndVariableImposedDeformation( ViscoelasticModel model, const Matrix & rig, const Function & f, int additionnalBlocksAfter = 0, double r = 0) ;

    virtual ~ViscoelasticityAndVariableImposedDeformation() { } ;

    virtual Form * getCopy() const ;

    virtual Vector getImposedStress(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual Vector getImposedStrain(const Point & p, IntegrableEntity * e = nullptr, int g = -1) const ;

    virtual std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

/*    virtual bool isViscous() const {
        return true ;
    }*/

} ;

}


#endif
