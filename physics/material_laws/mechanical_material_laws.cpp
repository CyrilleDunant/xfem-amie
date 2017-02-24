#include "mechanical_material_laws.h"
#include <stdlib.h>
#include <fstream>

namespace Amie
{

void BulkShearConversionMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt)
{
    if(s.has("bulk_modulus") && s.has("shear_modulus"))
    {
        if(! (s.has("young_modulus") && s.has("poisson_ratio") ) )
        {
            double k = s.get("bulk_modulus", defaultValues) ;
            double mu = s.get("shear_modulus", defaultValues) ;
            s.set("young_modulus", 9.*k*mu/(3.*k+mu) ) ;
            s.set("poisson_ratio", (3.*k-2.*mu)/(6.*k+2.*mu) ) ;
        }
    }
    else
    {
        double E = s.get("young_modulus", defaultValues) ;
        double nu = s.get("poisson_ratio", defaultValues) ;
        s.set("bulk_modulus", E/(3.-6.*nu) ) ;
        s.set("shear_modulus", E/(2.+2.*nu) ) ;
    }

    if(s.has("creep_characteristic_time"))
    {
        if(s.has("creep_bulk") && s.has("creep_shear"))
        {
            if(! (s.has("creep_modulus") && s.has("creep_poisson") ) )
            {
                double k = s.get("creep_bulk", defaultValues) ;
                double mu = s.get("creep_shear", defaultValues) ;
                s.set("creep_modulus", 9.*k*mu/(3.*k+mu) ) ;
                s.set("creep_poisson", (3.*k-2.*mu)/(6.*k+2.*mu) ) ;
            }
        }
        else
        {
            double E = s.get("creep_modulus", defaultValues) ;
            double nu = s.get("creep_poisson", defaultValues) ;
            s.set("creep_bulk", E/(3.-6.*nu) ) ;
            s.set("creep_shear", E/(2.+2.*nu) ) ;
        }

        if(s.has("recoverable_bulk") && s.has("recoverable_shear"))
        {
            if(! (s.has("recoverable_modulus") && s.has("recoverable_poisson") ) )
            {
                double k = s.get("recoverable_bulk", defaultValues) ;
                double mu = s.get("recoverable_shear", defaultValues) ;
                s.set("recoverable_modulus", 9.*k*mu/(3.*k+mu) ) ;
                s.set("recoverable_poisson", (3.*k-2.*mu)/(6.*k+2.*mu) ) ;
            }
        }
        else
        {
            double E = s.get("recoverable_modulus", defaultValues) ;
            double nu = s.get("recoverable_poisson", defaultValues) ;
            s.set("recoverable_bulk", E/(3.-6.*nu) ) ;
            s.set("recoverable_shear", E/(2.+2.*nu) ) ;
        }

    }
}

void BazantLoadNonLinearCreepMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("creep_modulus"))
        return ;

    Vector stress(3) ;
    stress=0. ;
    s.getAverageField(REAL_STRESS_FIELD, stress, nullptr, 1.) ;
    double pstress =  stress.max() ;
    double strength = s.get("tensile_strength", defaultValues) ;
    if(stress.min() < 0 && -stress.min() > stress.max())
    {
        pstress = stress.min() ;
        strength = s.get("compressive_strength", defaultValues) ;
    }
    double C = s.get("young_modulus",defaultValues) ; //s.get("creep_modulus",defaultValues) ;
    double alpha = s.get("imposed_deformation", defaultValues) ;
    pstress += C*alpha ;
    double f = pstress/strength ;
    double factor = (1.-std::pow(f, 10.))/(1.+f*f) ;
    s.multiply("creep_modulus", factor) ;
    s.multiply("recoverable_modulus", factor) ;
}

void AdjustStrainStressCurveMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
   double E = s.get("young_modulus", defaultValues) ;
   double sigma = s.get("tensile_strength", defaultValues) ;
   s.set("tensile_strain", sigma/E ) ;
   if(s.has("tensile_fracture_energy"))
   {
       double Gf = s.get("tensile_fracture_energy", defaultValues) ;
       double smax = s.get("tensile_ultimate_strength", defaultValues) ;
       double deltaEpsilon = 2.*Gf/std::abs(sigma-smax) ;
       s.set("tensile_ultimate_strain", sigma/E+deltaEpsilon) ;
   }
   else if(s.has("tensile_softening_modulus"))
   {
       double Esoft = s.get("tensile_softening_modulus", defaultValues) ;
       double smax = s.get("tensile_ultimate_strength", defaultValues) ;
       double deltaEpsilon = std::abs(sigma-smax)/Esoft ;
       s.set("tensile_ultimate_strain", sigma/E+deltaEpsilon) ;
   }
   else if(s.has("tensile_softening_strain"))
   {
       double deltaEpsilon = s.get("tensile_softening_strain", defaultValues) ;
       s.set("tensile_ultimate_strain", sigma/E+deltaEpsilon) ;
   }
}

void TensionCompressionCreepMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("creep_modulus"))
        return ;

    Vector stress(2) ;
    stress=0. ;
    s.getAverageField(PRINCIPAL_REAL_STRESS_FIELD, stress, nullptr, 1.) ;
    double tstress =  stress.max() ;
    double cstress =  stress.min() ;
    if(tstress > 0 && tstress > std::abs(cstress))
    {
        double coef = s.get("tension_compression_creep_coefficient", defaultValues) ;
        s.multiply("creep_modulus", coef) ;
        s.multiply("recoverable_modulus", coef) ;
    }
}

void StrainRateDependentStrengthMaterialLaw::preProcess(GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt)
{
    if(!s.has("tensile_strength") && !s.has("compressive_strength") )
        return ;

    if(!s.has("previous_principal_strain_0") || !s.has("previous_principal_strain_1"))
    {
        s.set("previous_principal_strain_0", 0) ;
        s.set("previous_principal_strain_1", 0) ;
    }

    Vector pstrain(2) ;
    pstrain=0. ;
    s.getAverageField(PRINCIPAL_TOTAL_STRAIN_FIELD, pstrain, nullptr, 1.) ;
    if(std::abs(pstrain.max()) < POINT_TOLERANCE*dt && std::abs(pstrain.min()) < POINT_TOLERANCE*dt)
        return ;

    double strainrate = 0. ;
    if( std::abs(pstrain[0]) > std::abs(pstrain[1])  )
        strainrate = pstrain[0]-s.get("previous_principal_strain_0",defaultValues) ;
    else
        strainrate = pstrain[1]-s.get("previous_principal_strain_1",defaultValues) ;

    if(strainrate < 0)
        strainrate *= -1. ;

    strainrate /= s.getNodalDeltaTime() ;

    if(strainrate < strainRateRef)
    {

        double f = strainrate/strainRateRef ;

        if(s.has("tensile_strength"))
        {
            s.multiply("tensile_strength", 0.6+0.4*std::pow(f, p)) ;
            s.multiply("ultimate_tensile_strain", 1.+std::log10(1./f)) ;
        }

        if(s.has("compressive_strength"))
        {
            s.multiply("compressive_strength", 0.6+0.4*std::pow(f, p)) ;
            s.multiply("ultimate_compressive_strain", 1.+std::log10(1./f)) ;
        }
    }

    s.set("previous_principal_strain_0",pstrain[0]);
    s.set("previous_principal_strain_1",pstrain[1]);

}

void SetInitialDamageMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables &s, double dt )
{
    if(s.has("initial_damage"))
    {
        double d = s.get("initial_damage", defaultValues ) ;
        DamageModel * dfunc = s.getParent()->getBehaviour()->getDamageModel() ;
        if(dfunc && d > 0)
        {
            if( dfunc->getState()[0] < d )
                dfunc->getState(true)[0] = d ;
        }
    }

}

MineralMaterialLaw::MineralMaterialLaw( Mineral min, Variable cut, std::string args, char sep) : ExternalMaterialLaw( args, sep ), mineral(min ), cuttingPlane(cut)
{

}

MineralMaterialLaw::MineralMaterialLaw( std::string filename, std::string sepstr, int index, double factor, bool force, Variable cut, std::string args, char sep) : ExternalMaterialLaw( args, sep ), mineral(filename, sepstr, index, factor, force ), cuttingPlane(cut)
{

}

void MineralMaterialLaw::preProcess( Matrix & stiffness, Point angle, planeType plane ) 
{
    if(mineral.isValid())
    {
        Matrix C = Tensor::rotate4thOrderTensor3D(mineral.stiffness, angle ) ;
        if( stiffness.numCols() == 3 )
        {
            stiffness = Tensor::to2D( C, plane, cuttingPlane ) ;
        }
        else
        {
             stiffness = C ;
        }
    }
} 

void SetStiffnessMatrixMaterialLaw::preProcess( Matrix & stiffness, Point angle, planeType plane ) 
{
    stiffness = localStiffness ;
} 

void BiLinearStiffnessMaterialLaw::preProcess( GeneralizedSpaceTimeViscoElasticElementStateWithInternalVariables & s, double dt )
{
    if(std::abs(factor-1) < POINT_TOLERANCE) { return ; }

    Vector strain(3) ;
    s.getAverageField(TOTAL_STRAIN_FIELD, strain) ;
    double criterion = strain[0]+strain[1] ;
    double cutoff = s.get(limit,defaultValues) ;

    if(cutoff > POINT_TOLERANCE)
    {
        if(criterion > cutoff && !over )
        {
            over = true ;
            if(std::abs(criterion) > POINT_TOLERANCE)
                residual = strain*cutoff/criterion;
            else
                residual = strain ;
            residual *= (1.-1./factor) ;
        }
        else if( criterion < cutoff && over)
        {
            over = false ;
            residual = 0 ;
        }
    }
    else if( cutoff < -POINT_TOLERANCE)
    {
        if(criterion > cutoff && !over)
        {
            over = true ;
            residual = 0 ;
        }
        else if(criterion < cutoff && over)
        {
            over = false ;
            if(std::abs(criterion) > POINT_TOLERANCE)
                residual = strain*cutoff/criterion;
            else
                residual = strain ;
            residual *= (1.-1./factor) ;
        }

    }
    else
    {
        over = criterion > 0 ;
        residual = 0 ;
    }

    s.add("imposed_deformation_xx",residual[0]) ;
    s.add("imposed_deformation_yy",residual[1]) ;
}

void BiLinearStiffnessMaterialLaw::preProcess( Matrix & stiffness, Point angle, planeType plane )
{
    if(over)
        stiffness *= factor ;

}


} 

