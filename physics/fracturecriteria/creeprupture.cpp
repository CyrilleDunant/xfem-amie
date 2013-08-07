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
	double kstrain = 1./(limStrain * damage) ;

	std::pair<Vector, Vector> stateBefore( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, -1) ) ;
	std::pair<Vector, Vector> stateAfter( smoothedPrincipalStressAndStrain(s, FROM_STRESS_STRAIN, REAL_STRESS, 1) ) ;

	double maxStrainBefore = stateBefore.second.max() ;
	double maxStressBefore = stateBefore.first.max() ;
	double maxStrainAfter = stateAfter.second.max() ;
	double maxStressAfter = stateAfter.first.max() ;

	double dStress = maxStressAfter - maxStressBefore ;
	double dStrain = maxStrainAfter - maxStrainBefore ;
	double ds = maxStress-limStress ;
	
	double cBefore = maxStressBefore - ds*exp( - maxStrainBefore * kstrain )  ;
	double cAfter = maxStressAfter - ds*exp( - maxStrainAfter * kstrain )  ;

	if(cBefore > limStress)
		return 1. ;
	if(cAfter < limStress)
		return cAfter / maxStress;

	float crit = 1. ;
	float tmin = 0. ;
	float tmax = 1. ;
	float tnext = 0.5 ; 
	int niter = 0 ;
	
	while(std::abs(crit) > 1e-6 && niter++ < 10)
	{
		tnext = tmin*0.66 + tmax*0.34 ;
		float sigma = maxStressBefore + tnext*dStress ;
		float epsilon = -(maxStrainBefore + tnext*dStrain)*kstrain ;
		crit = sigma - ds * expf( epsilon ) ;
		(crit > limStress) ? tmax = tnext : tmin = tnext ;
	}
	metInTension = true ;
	return 1.-tnext ;
	
}

Material CreepRupture::toMaterial()
{
	Material mat ;
	return mat ;
}



}
