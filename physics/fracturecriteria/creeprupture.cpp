#include "creeprupture.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../damagemodels/damagemodel.h"
namespace Mu
{

CreepRupture::CreepRupture(double max, double lim, double strain, MirrorState mirroring, double delta_x, double delta_y, double delta_z): FractureCriterion( mirroring, delta_x, delta_y, delta_z ), maxStress(max), limStress(lim), limStrain(strain)
{
	metInCompression = false ;
	metInTension = false ;
}

CreepRupture::~CreepRupture()
{

}

FractureCriterion * CreepRupture::getCopy() const
{
	return new CreepRupture( *this ) ;
}

double CreepRupture::grade( ElementState &s )
{
	if( s.getParent()->getBehaviour()->fractured())
		return -1 ;

	inIteration = true ;
	metInTension = false ;

	double damage = 0.001 ;
	if(s.getParent()->getBehaviour()->getDamageModel())
		damage = s.getParent()->getBehaviour()->getDamageModel()->getState().max() ;
	if(damage < 0.001)
		damage = 0.001 ;
	double kstrain = limStrain * damage ;

	std::pair<Vector, Vector> stateBefore( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, -1) ) ;
	std::pair<Vector, Vector> stateAfter( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, 1) ) ;

	double maxStrainBefore = stateBefore.second.max() ;
	double maxStressBefore = stateBefore.first.max() ;
	double maxStrainAfter = stateAfter.second.max() ;
	double maxStressAfter = stateAfter.first.max() ;

	double dStress = maxStressAfter - maxStressBefore ;
	double dStrain = maxStrainAfter - maxStrainBefore ;

	double cBefore = (maxStressBefore - (limStress + (maxStress-limStress)*exp( - maxStrainBefore / kstrain ) ))/maxStress ;
	double cAfter = (maxStressAfter - (limStress + (maxStress-limStress)*exp( - maxStrainAfter / kstrain ) ))/maxStress ;

	if(cBefore > 0)
		return 1. ;
	if(cAfter < 0)
		return cAfter ;

	double crit = 1. ;
	double tmin = 0. ;
	double tmax = 1. ;
	double tnext = 0.5 ; 
	int niter = 0 ;
	while(std::abs(crit) > 1e-6 && niter < 10)
	{
		tnext = (tmin + tmax)*0.5 ;
		double sigma = maxStressBefore + (tnext)*dStress ;
		double epsilon = maxStrainBefore +  (tnext)*dStrain ;
		crit = (sigma - (limStress + (maxStress-limStress)*exp( - epsilon / kstrain ) ))/maxStress ;
		if(crit > 0)
			tmax = tnext ;
		if(crit < 0)
			tmin = tnext ;
		niter ++ ;
	}
	metInTension = true ;
	return 1-tnext ;

	
}

Material CreepRupture::toMaterial()
{
	Material mat ;
	return mat ;
}



}
