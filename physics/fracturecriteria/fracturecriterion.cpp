//
// C++ Implementation: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fracturecriterion.h"
#include "../damagemodels/damagemodel.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../../solvers/assembly.h"
using namespace Mu ;

FractureCriterion::FractureCriterion(MirrorState mirroring, double delta_x, double delta_y, double delta_z) : neighbourhoodradius(.0005), 
neighbourhoodvolume(-1), 
physicalCharacteristicRadius(.008), 
scoreAtState(-1), 
nonLocalScoreAtState(-1),
metAtStep(false), 
mirroring(mirroring), delta_x(delta_x), delta_y(delta_y), delta_z(delta_z), 
deltaScoreAtState(0), 
deltaEnergyAtState(0), 
energyDamageDifferential(0), 
criterionDamageDifferential(0), 
energyIndexed(false), 
noEnergyUpdate(true), 
mesh2d(NULL), mesh3d(NULL), 
stable(true), checkpoint(true), inset(false),
scoreTolerance(.33e-2),
initialScore(1),
cachedInfluenceRatio(1),
currentAngle(0),
minDeltaInNeighbourhood(1),
maxModeInNeighbourhood(-1),
maxAngleShiftInNeighbourhood(0)
{
}

double FractureCriterion::smoothedCrackAngle( ElementState &s) const
{
	double angle = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double area = s.getParent()->area() ;
		double fact =  0 ; //area*s.getParent()->getBehaviour()->getDamageModel()->getState().max() ;
		
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			double dc =  squareDist2D( ci->getCenter(), s.getParent()->getCenter() ) ;
			if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() 
				|| !ci->getBehaviour()->getFractureCriterion() 
				|| ci->getBehaviour()->fractured()
				|| ci->getBehaviour()->type == VOID_BEHAVIOUR
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc >  4.*physicalCharacteristicRadius * physicalCharacteristicRadius)
			{
				continue ;
			}
			
			area = ci->area() ;
			//this is to eliminate scaling effects ;
			double factor = 1.;//-ci->getBehaviour()->getDamageModel()->getState().max() ;
// 			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D && std::abs(ci->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
// 				factor = std::min(std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]),std::abs(s.getParent()->getBehaviour()->param[0][0]/ci->getBehaviour()->param[0][0])) ;
			
			double d = ci->getBehaviour()->getDamageModel()->getState().max() ;
			angle += atan2(ci->getCenter().y-s.getParent()->getCenter().y,ci->getCenter().x-s.getParent()->getCenter().x)*area*d ;
			fact += area*d ;
			
			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				angle += atan2(ci->getCenter().y-s.getParent()->getCenter().y,ci->getCenter().x-s.getParent()->getCenter().x)*area*d ;
				fact += area*d ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				angle += atan2(ci->getCenter().y-s.getParent()->getCenter().y,ci->getCenter().x-s.getParent()->getCenter().x)*area*d ;
				fact += area*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				angle += atan2(ci->getCenter().y-s.getParent()->getCenter().y,ci->getCenter().x-s.getParent()->getCenter().x)*area*d ;
				fact += area*d ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				angle += atan2(ci->getCenter().y-s.getParent()->getCenter().y,ci->getCenter().x-s.getParent()->getCenter().x)*area*d ;
				fact += area*d ;
			}
		}
		
		if(fact > POINT_TOLERANCE_2D)
			angle /= (fact+ s.getParent()->area()*s.getParent()->getBehaviour()->getDamageModel()->getState().max()) ;
		else
			return 0. ;
	}
	return angle ;
}

double FractureCriterion::smoothedPrincipalStressAngle( ElementState &s, StressCalculationMethod m )
{
	double angle = 0 ;
	std::pair<Vector, Vector> ss = smoothedStressAndStrain(s, m) ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		return 0.5 * atan2( ss.first[2], ss.first[0] - ss.first[1] ) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		return  0.5 * atan2(ss.first[3] , ss.first[0] - ss.first[1] ) ;
	}
	
	return angle ;
}

Vector FractureCriterion::smoothedPrincipalStrain(ElementState &s)
{
	return smoothedPrincipalStressAndStrain(s).second ;
}

std::pair< Vector, Vector > FractureCriterion::smoothedPrincipalStressAndStrain(ElementState& s, StressCalculationMethod m , bool useStressLimit )
{

	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{

		std::pair<Vector, Vector> stressAndStrain = smoothedStressAndStrain(s, m) ;
		
		Vector lprincipal( 2 ) ;

		lprincipal[0] = 0.5 * ( stressAndStrain.first[0] + stressAndStrain.first[1] ) +
										 sqrt(
												0.25 *( stressAndStrain.first[0] - stressAndStrain.first[1] ) * ( stressAndStrain.first[0] - stressAndStrain.first[1] ) +
												( stressAndStrain.first[2] * stressAndStrain.first[2] )
										) ;
		
		lprincipal[1] = 0.5 * ( stressAndStrain.first[0] + stressAndStrain.first[1] ) -
										 sqrt(
												0.25 *( stressAndStrain.first[0] - stressAndStrain.first[1] ) * ( stressAndStrain.first[0] - stressAndStrain.first[1] ) +
												( stressAndStrain.first[2] * stressAndStrain.first[2] )
										) ;

		Vector lprincipals( 2 ) ;

		lprincipals[0] = 0.5 * ( stressAndStrain.second[0] + stressAndStrain.second[1] ) +
										 sqrt(
												0.25 *( stressAndStrain.second[0] - stressAndStrain.second[1] ) * ( stressAndStrain.second[0] - stressAndStrain.second[1] ) +
												( stressAndStrain.second[2] * stressAndStrain.second[2] )
										) ;
		
		lprincipals[1] = 0.5 * ( stressAndStrain.second[0] + stressAndStrain.second[1] ) -
										 sqrt(
												0.25 *( stressAndStrain.second[0] - stressAndStrain.second[1] ) * ( stressAndStrain.second[0] - stressAndStrain.second[1] ) +
												( stressAndStrain.second[2] * stressAndStrain.second[2] )
										) ;
		
		return std::make_pair(lprincipal, lprincipals) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		std::pair<Vector, Vector> stressAndStrain = smoothedStressAndStrain(s,m) ;
		
		Vector lprincipal( 3 ) ;
		Matrix stresses(3,3) ;
		stresses[0][0] = stressAndStrain.first[0] ; stresses[0][1] = stressAndStrain.first[3] ; stresses[0][2] = stressAndStrain.first[4] ;
		stresses[1][0] = stressAndStrain.first[3] ; stresses[1][1] = stressAndStrain.first[1] ; stresses[1][2] = stressAndStrain.first[5] ;
		stresses[2][0] = stressAndStrain.first[4] ; stresses[2][1] = stressAndStrain.first[5] ; stresses[2][2] = stressAndStrain.first[2] ;
		Matrix I( 3, 3 ) ;
		I[0][0] = 1 ;
		I[1][1] = 1 ;
		I[2][2] = 1 ;
		double m = ( stresses[0][0] + stresses[1][1] + stresses[2][2] ) / 3. ;
		Matrix Am = stresses - I * m ;
		double q = det( Am ) / 2. ;
		double rr = std::inner_product( &Am.array()[0], &Am.array()[9], &Am.array()[0],  double( 0. ) ) / 6. ;
		double phi = atan2( sqrt( rr * rr * rr - q * q ), q ) / 3. ;

		if( rr * rr * rr - q * q < 1e-12 )
			phi = atan( 0 ) / 3. ;

		if( phi < 0 )
			phi += M_PI ;

		
		lprincipal[0] = m + 2.*sqrt( rr ) * cos( phi ) ;
		lprincipal[1] = m - sqrt( rr ) * ( cos( phi ) + sqrt( 3 ) * sin( phi ) ) ;
		lprincipal[2] = m - sqrt( rr ) * ( cos( phi ) - sqrt( 3 ) * sin( phi ) ) ;
		
		Vector lprincipals( 3 ) ;
		Matrix sstresses(3,3) ;
		sstresses[0][0] = stressAndStrain.second[0] ; sstresses[0][1] = stressAndStrain.second[3] ; sstresses[0][2] = stressAndStrain.second[4] ;
		sstresses[1][0] = stressAndStrain.second[3] ; sstresses[1][1] = stressAndStrain.second[1] ; sstresses[1][2] = stressAndStrain.second[5] ;
		sstresses[2][0] = stressAndStrain.second[4] ; sstresses[2][1] = stressAndStrain.second[5] ; sstresses[2][2] = stressAndStrain.second[2] ;
		m = ( sstresses[0][0] + sstresses[1][1] + sstresses[2][2] ) / 3. ;
		Am = sstresses - I * m ;
		q = det( Am ) / 2. ;
		rr = std::inner_product( &Am.array()[0], &Am.array()[9], &Am.array()[0],  double( 0. ) ) / 6. ;
		phi = atan2( sqrt( rr * rr * rr - q * q ), q ) / 3. ;

		if( rr * rr * rr - q * q < 1e-12 )
			phi = atan( 0 ) / 3. ;

		if( phi < 0 )
			phi += M_PI ;

		
		lprincipals[0] = m + 2.*sqrt( rr ) * cos( phi ) ;
		lprincipals[1] = m - sqrt( rr ) * ( cos( phi ) + sqrt( 3 ) * sin( phi ) ) ;
		lprincipals[2] = m - sqrt( rr ) * ( cos( phi ) - sqrt( 3 ) * sin( phi ) ) ;
		
		return std::make_pair(lprincipal,lprincipals)  ;
	}
	
	return std::make_pair(smoothedPrincipalStress(s, m), smoothedPrincipalStrain(s)) ;
}

double FractureCriterion::smoothedScore(ElementState& s)
{

	double score = 0;
	if(factors.empty())
		initialiseFactors(s) ;
	double total = 0 ;
	auto fiterator = factors.begin() ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		score =scoreAtState*(*fiterator) ;
		total += (*fiterator);
		fiterator++ ;
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			if(ci->getBehaviour()->fractured())
			{
				fiterator++ ;
				continue ;
			}

			score += ci->getBehaviour()->getFractureCriterion()->scoreAtState*(*fiterator) ;
			total += (*fiterator);
			fiterator++ ;
		}
		score /= total ;

		return score ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		score =scoreAtState*(*fiterator) ;
		total += (*fiterator);
		fiterator++ ;
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			
			if( ci->getBehaviour()->fractured() || *fiterator < POINT_TOLERANCE_2D)
			{
				fiterator++ ;
				continue ;
			}
			score += ci->getBehaviour()->getFractureCriterion()->scoreAtState*(*fiterator) ;
			total += (*fiterator);
			fiterator++ ;
		}
		score /= total ;
		return score ;
	}
	
	return 0. ;
}

double FractureCriterion::getSquareInfluenceRatio(ElementState & s, const Point & direction) 
{

// 	return 1. ;
	double dirnorm = direction.norm() ;
	if(dirnorm < POINT_TOLERANCE_2D)
	{
		cachedInfluenceRatio = 1 ;
		return 1. ;
	}
	Point renormdir = direction/-dirnorm ;
	
	if(/*!getInfluenceCalculated()*/true)
	{
		std::vector<Point> u = s.getPrincipalFrame(s.getParent()->getCenter()) ;
		Point tricross = u[2]^(renormdir^u[2]) ;
		double tricrossNorm = tricross.norm() ;
		if(tricrossNorm < POINT_TOLERANCE_2D)
		{
			cachedInfluenceRatio = 1 ;
			return 1. ;
		}
		tricross /= tricrossNorm ;

		double costheta = u[0]*tricross ; costheta *= costheta ;
		double sintheta = u[1]*tricross ; sintheta *= sintheta ;
		double cosphi = u[2]*renormdir ; cosphi *= cosphi ;
		double sinphi = renormdir*tricross ; sinphi *= sinphi ;
		
	
		Vector sigman = s.getPrincipalStressAtNodes() ;
		Vector sigma(0., sigman.size()/s.getParent()->getBoundingPoints().size()) ;
		for(int i = 0 ; i < s.getParent()->getBoundingPoints().size() ; ++i )
		{
			for(int j = 0 ; j <  sigma.size() ; ++j)
			{
				sigma[j] += sigman[i*sigma.size()+j]/s.getParent()->getBoundingPoints().size() ;
			}
		}
		double limit = getTensileLimit(s) ;
		if(std::abs(sigma).min() < POINT_TOLERANCE_2D)
		{
			cachedInfluenceRatio = 0 ;
			return 0. ;
		}
		if(sigma.max() > limit)
		{
			cachedInfluenceRatio = 1 ;
			return 1. ;
		}
		if(sigma.size() == 2)
		{
			double minrho = std::min(s.getParent()->area()/(physicalCharacteristicRadius*physicalCharacteristicRadius*4.), 1.) ;
			double costheta = u[0]*tricross ; costheta *= costheta ;
			double sintheta = u[1]*tricross ; sintheta *= sintheta ;
			double cosphi = u[2]*renormdir ; cosphi *= cosphi ;
			double sinphi = renormdir*tricross ; sinphi *= sinphi ;
			limit =1 ;
			sigma = 1 ;
			cachedInfluenceRatio = std::max(std::min(1./(limit*limit*(sinphi*costheta/(sigma[0]*sigma[0]) + sinphi*sintheta/(sigma[1]*sigma[1]))), 1.), minrho) ;
			
			return cachedInfluenceRatio ;
		}
		if(sigma.size() == 3)
		{
			double minrho = pow(s.getParent()->volume(), .66666666666)/(physicalCharacteristicRadius*physicalCharacteristicRadius*2.) ;
			double costheta = u[0]*tricross ; costheta *= costheta ;
			double sintheta = u[1]*tricross ; sintheta *= sintheta ;
			double cosphi = u[2]*renormdir ; cosphi *= cosphi ;
			double sinphi = renormdir*tricross ; sinphi *= sinphi ;
			limit =1 ;
			sigma = 1 ;
			cachedInfluenceRatio = std::max(std::min(1./(limit*limit*(sinphi*costheta/(sigma[0]*sigma[0]) + sinphi*sintheta/(sigma[1]*sigma[1]) + cosphi/(sigma[2]*sigma[2]))), 1.), minrho) ;
			return cachedInfluenceRatio ;
		}
		cachedInfluenceRatio = 0 ;
		return 0. ;
	}
	else
		return cachedInfluenceRatio ;
	}

void FractureCriterion::initialiseFactors(const ElementState & s)
{
	if(cache.empty())
		initialiseCache(s);
	if(!factors.empty())
		factors.clear() ;
	
	VirtualMachine vm ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		
		Function x = s.getParent()->getXTransform()-s.getParent()->getCenter().x ;
		Function y = s.getParent()->getYTransform()-s.getParent()->getCenter().y ;
		Function rr = x*x+y*y ;
		Function rrn =  rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
		Order order = s.getParent()->getOrder() ;
		s.getParent()->setOrder(CUBIC) ;
		Function smooth =  (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
		double weight = vm.ieval(smooth, s.getParent()) ;
		s.getParent()->setOrder(order) ;
		double fact = weight ;
		factors.push_back(weight);
		
		double selfarea = s.getParent()->area() ;
		double farthest = 0 ;
		if( mirroring == MIRROR_X && std::abs( s.getParent()->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
		{
			x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().x  - delta_x )-s.getParent()->getCenter().x ;
			y = s.getParent()->getYTransform()-s.getParent()->getCenter().y ;
			rr = x*x+y*y ;
			rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth =  (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight ;
			fact += weight ;
		}
		if( mirroring == MIRROR_Y &&  std::abs( s.getParent()->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
		{
			x = s.getParent()->getXTransform() ;
			y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().y  - delta_y )-s.getParent()->getCenter().x ;
			rr = x*x+y*y ;
			rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth =  (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight ;
			fact += weight ;
		}
		if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().x  - delta_x )-s.getParent()->getCenter().x ;
			y = s.getParent()->getYTransform()-s.getParent()->getCenter().y ;
			rr = x*x+y*y ;
			rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth =  (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight ;
			fact += weight ;
		}
		if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform()-s.getParent()->getCenter().x ;
			y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().y  - delta_y )-s.getParent()->getCenter().y ;
			rr = x*x+y*y ;
			rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight ;
			fact += weight ;
		}

		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			if(dynamic_cast<IntegrableEntity *>( ci ) == s.getParent() 
				|| ci->getBehaviour()->type == VOID_BEHAVIOUR
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() )
			{
				factors.push_back(0.);
				continue ;
			}
			double mindist = physicalCharacteristicRadius ;
// 			if( mirroring == NO_MIRROR || 
// 				mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) >= physicalCharacteristicRadius ||
// 				mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) >= physicalCharacteristicRadius ||
// 				mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) >= physicalCharacteristicRadius ||
// 				mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) >= physicalCharacteristicRadius
// 			)
// 			{
// 				Line interceptLine(s.getParent()->getCenter(), ci->getCenter()-s.getParent()->getCenter()) ;
// 				std::vector<Point> intersects = interceptLine.intersection(s.getParent()->getBehaviour()->getSource()) ;
// 				
// 				for(size_t j = 0 ; j < intersects.size() ; j++)
// 				{
// 					double dtest = dist(intersects[j], s.getParent()->getCenter()) ;
// 					if( dtest < mindist && (s.getParent()->getCenter()-ci->getCenter())*(ci->getCenter()-intersects[j]) > 0)
// 					{
// 						mindist = dtest ;
// 					}
// 				}
// 			}
			
			//this is to eliminate scaling effects ;
			double factor = 1.;//-ci->getBehaviour()->getDamageModel()->getState().max() ; ;
// 			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D && std::abs(ci->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
// 				factor = std::min(std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]),std::abs(s.getParent()->getBehaviour()->param[0][0]/ci->getBehaviour()->param[0][0])) ;
			double d = dist(ci->getCenter(), s.getParent()->getCenter()) ; 
			if(d > farthest)
				farthest = d ;
			x = ci->getXTransform()-s.getParent()->getCenter().x ;
			y = ci->getYTransform()-s.getParent()->getCenter().y ;
			rr = x*x+y*y ;
			rrn = rr/(mindist * mindist * -2.) ;
			order = s.getParent()->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.push_back(weight*factor);
			fact += weight*factor ;
			
			if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
			{
				x = ci->getXTransform()*-1 -std::abs( s.getParent()->getCenter().x  - delta_x )-s.getParent()->getCenter().x;
				y = ci->getYTransform()-s.getParent()->getCenter().y ;
				rr = x*x+y*y ;
				rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				order = s.getParent()->getOrder() ;
				ci->setOrder(CUBIC) ;
				smooth = (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				weight = vm.ieval(smooth, ci) ;
				factors.back() += weight ;
				ci->setOrder(order) ;
				fact += weight ;
			}

			if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
			{
				x = ci->getXTransform()-s.getParent()->getCenter().x ;
				y = ci->getYTransform()*-1-std::abs( s.getParent()->getCenter().y  - delta_y )-s.getParent()->getCenter().y ;
				rr = x*x+y*y ;
				rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				order = s.getParent()->getOrder() ;
				ci->setOrder(CUBIC) ;
				smooth = (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				weight = vm.ieval(smooth, ci) ;
				ci->setOrder(order) ;
				factors.back() += weight ;
				fact += weight ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				x = ci->getXTransform()*-1 -std::abs( s.getParent()->getCenter().x  - delta_x )-s.getParent()->getCenter().x;
				y = ci->getYTransform()-s.getParent()->getCenter().y ;
				rr = x*x+y*y ;
				rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				order = s.getParent()->getOrder() ;
				ci->setOrder(CUBIC) ;
				smooth = (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				weight = vm.ieval(smooth, ci) ;
				ci->setOrder(order) ;
				factors.back() += weight ;
				fact += weight ;
			}

			if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
			{
				x = ci->getXTransform()-s.getParent()->getCenter().x ;
				y = ci->getYTransform()*-1-std::abs( s.getParent()->getCenter().y  - delta_y )-s.getParent()->getCenter().y ;
				rr = x*x+y*y ;
				rrn = rr/(physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				order = s.getParent()->getOrder() ;
				ci->setOrder(CUBIC) ;
				smooth = (rrn-1.)*(rrn-1.)*f_negativity(rr-physicalCharacteristicRadius * physicalCharacteristicRadius) ;
				weight = vm.ieval(smooth, ci) ;
				ci->setOrder(order) ;
				factors.back() += weight ;
				fact += weight ;
			}

		}
// 		rrn = f_exp(Function("x 2 ^ y 2 ^ +")/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
// 		double totalAreatarget = 0 ;
// 		double da = farthest*farthest/(10.*10) ;
// 		for(double x = -farthest ; x < farthest ; x += farthest/10.)
// 		{
// 			for(double y = -farthest ; y < farthest ; y += farthest/10.)
// 			{
// 				if(x*x+y*y < farthest*farthest && !s.getParent()->in(Point(x, y)))
// 				{
// 					totalAreatarget += vm.eval(rrn, x, y)*da ;
// 				}
// 			}
// 		}
		
// 		std::cout << totalAreatarget << "  " << fact-factors[0] << "  "<<std::flush ;
// 		double ratio = 1 ; //( fact-factors[0])/totalAreatarget ;
// // 		std::cout << ratio << std::endl ;
// 		double w = factors[0] ;
// 		for(size_t i = 1 ; i < factors.size() ; i++)
// 		{
// 			factors[i] /= ratio ;
// 			w += factors[i] ;
// 		}
// 		factors.push_back(w);
		return ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		Function x = s.getParent()->getXTransform() ;
		Function y = s.getParent()->getYTransform() ;
		Function z = s.getParent()->getZTransform() ;
		Function r(s.getParent()->getCenter(), x, y, z) ;
		Order order = s.getParent()->getOrder() ;
		s.getParent()->setOrder(CUBIC) ;
		Function smooth = f_exp((r*r)/(-2.* physicalCharacteristicRadius * physicalCharacteristicRadius)) ;
		double weight = vm.ieval(smooth, s.getParent()) ;
		s.getParent()->setOrder(order) ;
		factors.push_back(weight);
		double fact = weight ;
		
		if( mirroring == MIRROR_X && std::abs( s.getParent()->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
		{
			x = s.getParent()->getXTransform()*-1.+std::abs( s.getParent()->getCenter().x  - delta_x ) ;
			y = s.getParent()->getYTransform() ;
			z = s.getParent()->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_Y &&  std::abs( s.getParent()->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
		{
			x = s.getParent()->getXTransform() ;
			y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().y  - delta_y ) ;
			z = s.getParent()->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_Z &&  std::abs( s.getParent()->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
		{
			x = s.getParent()->getXTransform() ;
			y = s.getParent()->getYTransform() ;
			z = s.getParent()->getZTransform()*-1-std::abs( s.getParent()->getCenter().z  - delta_z ) ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2. ));
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().x  - delta_x ) ;
			y = s.getParent()->getYTransform() ;
			z = s.getParent()->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(-2.* physicalCharacteristicRadius * physicalCharacteristicRadius)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			factors.back() += weight;
			s.getParent()->setOrder(order) ;
			fact += weight ;
		}

		if( mirroring == MIRROR_XY &&  std::abs( s.getParent()->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform() ;
			y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().y  - delta_y ) ;
			z = s.getParent()->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			factors.back() += weight;
			s.getParent()->setOrder(order) ;
			fact += weight ;
		}

		if( mirroring == MIRROR_XZ &&  std::abs( s.getParent()->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform()*-1-std::abs( s.getParent()->getCenter().x  - delta_x ) ;
			y = s.getParent()->getYTransform() ;
			z = s.getParent()->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_XZ &&  std::abs( s.getParent()->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform() ;
			y = s.getParent()->getYTransform() ;
			z = s.getParent()->getZTransform()*-1-std::abs( s.getParent()->getCenter().z  - delta_z ) ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_YZ &&  std::abs( s.getParent()->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform() ;
			y = s.getParent()->getYTransform()*-1-std::abs( s.getParent()->getCenter().y  - delta_y ) ;
			z = s.getParent()->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_YZ &&  std::abs( s.getParent()->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = s.getParent()->getXTransform() ;
			y = s.getParent()->getYTransform() ;
			z = s.getParent()->getZTransform()*-1-std::abs( s.getParent()->getCenter().z  - delta_z ) ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			s.getParent()->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}
		
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			double dc = squareDist3D( ci->getCenter(), s.getParent()->getCenter() ) ;
			
			if( dynamic_cast<IntegrableEntity *>( ci ) == s.getParent()  
				|| ci->getBehaviour()->getFractureCriterion() 
				|| ci->getBehaviour()->type == VOID_BEHAVIOUR
				|| ci->getBehaviour()->getSource() != s.getParent()->getBehaviour()->getSource() 
				|| dc > 4.* physicalCharacteristicRadius * physicalCharacteristicRadius
				|| ci->getBehaviour()->fractured()
			)
			{
				factors.push_back(0);
				continue ;
			}

			double factor = 1.;//-ci->getBehaviour()->getDamageModel()->getState().max() ; ;
// 			if(std::abs(s.getParent()->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D && std::abs(ci->getBehaviour()->param[0][0]) > POINT_TOLERANCE_3D)
// 				factor = std::min(std::abs(ci->getBehaviour()->param[0][0]/s.getParent()->getBehaviour()->param[0][0]),std::abs(s.getParent()->getBehaviour()->param[0][0]/ci->getBehaviour()->param[0][0])) ;

			x = ci->getXTransform() ;
			y = ci->getYTransform() ;
			z = ci->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.push_back(weight);
			fact += weight ;
			
		if( mirroring == MIRROR_X && std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_X
		{
			x = ci->getXTransform()*-1-std::abs( ci->getCenter().x  - delta_x ) ;
			y = ci->getYTransform() ;
			z = ci->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = s.getParent()->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_Y &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_Y
		{
			x = ci->getXTransform() ;
			y = ci->getYTransform()*-1-std::abs( ci->getCenter().y  - delta_y ) ;
			z = ci->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_Z &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_Y
		{
			x = ci->getXTransform() ;
			y = ci->getYTransform() ;
			z = ci->getZTransform()*-1-std::abs( ci->getCenter().z  - delta_z ) ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = ci->getXTransform()*-1-std::abs( ci->getCenter().x  - delta_x ) ;
			y = ci->getYTransform() ;
			z = ci->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			factors.back() += weight;
			ci->setOrder(order) ;
			fact += weight ;
		}

		if( mirroring == MIRROR_XY &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = ci->getXTransform() ;
			y = ci->getYTransform()*-1-std::abs( ci->getCenter().y  - delta_y ) ;
			z = ci->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			factors.back() += weight;
			ci->setOrder(order) ;
			fact += weight ;
		}

		if( mirroring == MIRROR_XZ &&  std::abs( ci->getCenter().x  - delta_x ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = ci->getXTransform()*-1-std::abs( ci->getCenter().x  - delta_x ) ;
			y = ci->getYTransform() ;
			z = ci->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_XZ &&  std::abs( s.getParent()->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = ci->getXTransform() ;
			y = ci->getYTransform() ;
			z = ci->getZTransform()*-1-std::abs( ci->getCenter().z  - delta_z ) ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, s.getParent()) ;
			ci->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().y  - delta_y ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = ci->getXTransform() ;
			y = ci->getYTransform()*-1-std::abs( ci->getCenter().y  - delta_y ) ;
			z = ci->getZTransform() ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			s.getParent()->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}

		if( mirroring == MIRROR_YZ &&  std::abs( ci->getCenter().z  - delta_z ) < physicalCharacteristicRadius )   // MIRROR_XY
		{
			x = ci->getXTransform() ;
			y = ci->getYTransform() ;
			z = ci->getZTransform()*-1-std::abs( ci->getCenter().z  - delta_z ) ;
			r = Function(s.getParent()->getCenter(), x, y, z) ;
			order = ci->getOrder() ;
			ci->setOrder(CUBIC) ;
			smooth = f_exp((r*r)/(physicalCharacteristicRadius * physicalCharacteristicRadius*-2.)) ;
			weight = vm.ieval(smooth, ci) ;
			ci->setOrder(order) ;
			factors.back() += weight;
			fact += weight ;
		}
		
		}
		factors.push_back(fact);

		return ;
	}
	
}

Vector FractureCriterion::smoothedPrincipalStress( ElementState &s, StressCalculationMethod m)
{
	return smoothedPrincipalStressAndStrain(s, m).first ;
}

std::pair<Vector, Vector> FractureCriterion::smoothedStressAndStrain( ElementState &s, StressCalculationMethod m, bool useStressLimit )
{
// 	useStressLimit = true ;

	size_t vlength = 3 ;
	if(s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
		vlength = 6 ;
	Vector str(0., vlength) ;
	Vector estr(0., vlength) ;
	Vector stra(0., vlength) ;
	Vector tmpstr(vlength) ;
	Vector tmpstra(vlength) ;

	if(factors.empty())
		initialiseFactors(s) ;
	double sumStressFactors = 0 ;
	double sumStrainFactors = 0 ;
	if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		double iteratorValue = factors[0] ;
		s.getStressAndStrainAtCenter(tmpstr, tmpstra, REAL_STRESS) ;
		stra = tmpstra*iteratorValue ;
		str = tmpstr*iteratorValue ;
		if(m == EFFECTIVE_STRESS)
		{
			s.getStressAndStrainAtCenter(tmpstr, tmpstra, EFFECTIVE_STRESS) ;
			estr = tmpstr*iteratorValue ;
		}
		sumStressFactors += iteratorValue ;
		sumStrainFactors += iteratorValue ;

		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTriangle *ci = static_cast<DelaunayTriangle *>( ( *mesh2d )[cache[i]] ) ;
			if(!ci->getBehaviour()->getFractureCriterion())
				factors[i+1] = 0 ;
			
			iteratorValue = factors[i+1] ;
			
			double fractureDistance = 0. ;
			
// 			if(ci->getBehaviour()->getDamageModel() && s.getParent()->getBehaviour()->getDamageModel())
// 			{
// 				Vector dstate = std::abs(ci->getBehaviour()->getDamageModel()->getState() - s.getParent()->getBehaviour()->getDamageModel()->getState()) ;
// 				fractureDistance = sqrt(std::inner_product(&dstate[0], &dstate[dstate.size()], &dstate[0], 0.)) ;
// 			}
			
// 			if(ci->getBehaviour()->getDamageModel() && ci->getBehaviour()->getDamageModel()->fractured())
// 				iteratorValue = 0 ;

// 			iteratorValue *= 1.-fractureDistance ;
			
			
			if(iteratorValue > POINT_TOLERANCE_2D)
			{
				ci->getState().getStressAndStrainAtCenter(tmpstr, tmpstra, REAL_STRESS) ;
				if(useStressLimit && ci->getBehaviour()->getFractureCriterion())
					iteratorValue = pow(iteratorValue, 1./ci->getBehaviour()->getFractureCriterion()->getSquareInfluenceRatio(ci->getState(),ci->getCenter()-s.getParent()->getCenter())) ;
				
				if(!ci->getBehaviour()->fractured())
				{
					stra += tmpstra*iteratorValue ;
					str += tmpstr*iteratorValue ;
					
					if(m == EFFECTIVE_STRESS)
					{
						ci->getState().getStressAndStrainAtCenter(tmpstr, tmpstra, EFFECTIVE_STRESS) ;
						estr += tmpstr*iteratorValue ;
					}
					
				}
// 				else
// 				{
// 				  stra += tmpstra*iteratorValue ;
// 				}
				
			}

			if(ci->getBehaviour()->getDamageModel())
			{
				
				if(!ci->getBehaviour()->fractured())
				{
					sumStrainFactors += iteratorValue ;
				  sumStressFactors += iteratorValue ;
				}
			}
		}

		str /= sumStressFactors ;
		estr /= sumStressFactors ;
		stra /= sumStrainFactors ;
		if(std::abs(stra[0]-stra[1]) > POINT_TOLERANCE_2D)
		{
			currentAngle =  0.5*atan2( stra[2], stra[0] - stra[1] ) ;
			if(currentAngle < 0)
				currentAngle += M_PI ;
		}

		if(m == REAL_STRESS)
		{
			return std::make_pair(str, stra) ;
		}
		return std::make_pair(estr, stra) ;
	}
	else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
	{
		double iteratorValue = factors[0];
		s.getStressAndStrainAtCenter(tmpstr, tmpstra, m) ;
		str = tmpstr*iteratorValue ;
		stra =tmpstra*iteratorValue ;
		sumStressFactors += iteratorValue ;
		sumStrainFactors += iteratorValue ;
		for( size_t i = 0 ; i < cache.size() ; i++ )
		{
			DelaunayTetrahedron *ci = static_cast<DelaunayTetrahedron *>( ( *mesh3d )[cache[i]] ) ;
			ci->getState().getStressAndStrainAtCenter(tmpstr, tmpstra, m) ;
			iteratorValue = factors[i+1] ;
			Point direction =  ci->getCenter()-s.getParent()->getCenter() ; 
			if(useStressLimit && ci->getBehaviour()->getFractureCriterion())
				iteratorValue = pow(iteratorValue, 1./ci->getBehaviour()->getFractureCriterion()->getSquareInfluenceRatio(ci->getState(),direction)) ;
			if( ci->getBehaviour()->fractured())
			{
				stra += tmpstra*iteratorValue ;
				sumStrainFactors += iteratorValue ;
				continue ;
			}
			str += tmpstr*iteratorValue ;
			stra += tmpstra*iteratorValue ;
			sumStressFactors += iteratorValue ;
			sumStrainFactors += iteratorValue ;
		}
		str /= sumStressFactors ;
		stra /= sumStrainFactors ;
		currentAngle = 0.5 * atan2(str[3] , str[0] - str[1] ) ;
		return std::make_pair(str,stra)  ;
	}
	
	return std::make_pair(smoothedPrincipalStress(s, m), smoothedPrincipalStrain(s)) ;
}

double FractureCriterion::getDeltaEnergy(const ElementState & s, double delta_d)
{
	Assembly K ;
	double originalscore = 0 ;
	double originalenergy = 0 ;
	double score = 0 ;
	double energy = 0 ;
	double volume = 0 ;

	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = static_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				originalenergy += elements[i]->getState().elasticEnergy()*a ;
			}
		}
		
		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
		
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron *>((*mesh3d)[ cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*3];
									double ey = tri->getState().getDisplacements()[l*3+1];
									double ez = tri->getState().getDisplacements()[l*3+2];
									K.setPoint(ex, ey, ez ,id);
									break ;
								}
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double v = elements[i]->volume() ;
				originalenergy += elements[i]->getState().elasticEnergy()*v ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle *>((*mesh2d)[ cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;

			if(tric == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;

		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				energy += elements[i]->getState().elasticEnergy()*a ;
			}
		}

		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * >((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			if(tet == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								int id = tri->getBoundingPoint(k).id ;
								double ex = tri->getState().getDisplacements()[k*3];
								double ey = tri->getState().getDisplacements()[k*3+1];
								double ez = tri->getState().getDisplacements()[k*3+2];
								K.setPoint(ex, ey, ez ,id);
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	return (originalenergy-energy)/(delta_d) ;
}

void FractureCriterion::initialiseCache(const ElementState & s)
{
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	if(testedTri)
	{
		if(!cache.empty())
		{
			cache.clear();
		}
		Circle epsilon(std::max(std::max(neighbourhoodradius, testedTri->getRadius()*2.1), physicalCharacteristicRadius*2.1),testedTri->getCenter()) ;
		if(!testedTri->tree)
			return ;
		mesh2d = &testedTri->tree->getTree() ;
		std::vector<DelaunayTriangle *> tempcache = testedTri->tree->getNeighbouringElementsInGeometry(testedTri, &epsilon);
		std::vector<DelaunayTriangle *> neighbourhood ;
		for(size_t i = 0 ; i < testedTri->neighbourhood.size() ; i++)
		{
			if(testedTri->getNeighbourhood(i)->getBehaviour() && 
				testedTri->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR && 
				testedTri->getNeighbourhood(i)->getBehaviour()->getFractureCriterion())
			{
				neighbourhood.push_back(testedTri->getNeighbourhood(i));
				cache.push_back(testedTri->getNeighbourhood(i)->index);
			}
		}
		
		
		for(size_t i = 0 ; i < tempcache.size() ; i++)
		{
			if(tempcache[i]->getBehaviour() && 
				tempcache[i]->getBehaviour()->type != VOID_BEHAVIOUR && 
				tempcache[i]->getBehaviour()->getFractureCriterion())
			{
				bool inNeighbourhood = false ;
				for(size_t j = 0 ; j < neighbourhood.size() ; j++)
				{
					if(neighbourhood[j] == tempcache[i])
					{
						inNeighbourhood = true ;
						break ;
					}
				}
				if(!inNeighbourhood)
				{
					cache.push_back(tempcache[i]->index);
				}
			}
		}
		
		if(cache.empty())
			cache.push_back(testedTri->index);
	}
	else if(testedTet)
	{
		if(!cache.empty())
			cache.clear();
		
		Sphere epsilon(neighbourhoodradius,testedTet->getCenter()) ;
		if(!testedTet->tree)
			return ;
		mesh3d = &testedTet->tree->getTree() ;
		std::vector<DelaunayTetrahedron *> tempcache3d = testedTet->tree->getConflictingElements(&epsilon);
		std::vector<DelaunayTetrahedron *> neighbourhood ;
		for(size_t i = 0 ; i < testedTet->neighbourhood.size() ; i++)
		{
			if(testedTet->getNeighbourhood(i)->getBehaviour()
				&& testedTet->getNeighbourhood(i)->getBehaviour()->type != VOID_BEHAVIOUR  
				&& testedTet->getNeighbourhood(i)->getBehaviour()->getFractureCriterion())
			{
				neighbourhood.push_back(testedTet->getNeighbourhood(i));
				cache.push_back(testedTet->getNeighbourhood(i)->index);
			}
		}
		
		for(size_t i = 0 ; i < tempcache3d.size() ; i++)
		{
			if(tempcache3d[i]->getBehaviour()
				&& tempcache3d[i]->getBehaviour()->type != VOID_BEHAVIOUR 
				&& tempcache3d[i]->getBehaviour()->getFractureCriterion())
			{
				bool inNeighbourhood = false ;
				for(size_t j = 0 ; j < neighbourhood.size() ; j++)
				{
					if(neighbourhood[j] == tempcache3d[i])
					{
						inNeighbourhood = true ;
						break ;
					}
				}
				if(!inNeighbourhood)
				{
					cache.push_back(tempcache3d[i]->index);
				}
			}
		}
		if(cache.empty())
			cache.push_back(testedTet->index);
	}
	
	std::sort(cache.begin(), cache.end()) ;
}

std::pair<double, double> FractureCriterion::getDeltaEnergyDeltaCriterion(const ElementState & s, double delta_d) const
{
	Assembly K ;
	double originalscore = 0 ;
	double originalenergy = 0 ;
	double score = 0 ;
	double energy = 0 ;
	double volume = 0 ;

	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,  tric->second,   tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j < tric->neighbour.size() ; j++)
			{
				if(tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>(tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
					originalscore += 1./(1.-sc) ;
				originalenergy += elements[i]->getState().elasticEnergy()*a ;
			}
		}
		
		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
		
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = dynamic_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*3];
									double ey = tri->getState().getDisplacements()[l*3+1];
									double ez = tri->getState().getDisplacements()[l*3+2];
									K.setPoint(ex, ey, ez ,id);
									break ;
								}
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double v = elements[i]->volume() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
					originalscore += 1./(1.-sc) ;
				originalenergy += elements[i]->getState().elasticEnergy()*v ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	
	if(s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;

		std::vector<DelaunayTriangle *> elements ;
		std::vector<LinearForm *> behaviours ;
		TriElement father(LINEAR) ;

		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTriangle * tric = static_cast<DelaunayTriangle * > ((*mesh2d)[cache[i]]) ;
			elements.push_back(new DelaunayTriangle(NULL, NULL,   tric->first,   tric->second,    tric->third,  NULL) );
			elements.back()->setBehaviour(tric->getBehaviour()->getCopy()) ;

			if(tric == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;

			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume += elements.back()->area() ;
			for(size_t j = 0 ; j <  tric->neighbour.size() ; j++)
			{
				if( tric->getNeighbour(j)->isTriangle )
				{
					DelaunayTriangle * tri = dynamic_cast<DelaunayTriangle *>( tric->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								for(size_t l = 0 ; l <  tri->getBoundingPoints().size() ; l++)
								{
									int id = tri->getBoundingPoint(k).id ;
									double ex = tri->getState().getDisplacements()[l*2];
									double ey = tri->getState().getDisplacements()[l*2+1];
									K.setPoint(ex, ey ,id);
									break ;
								}
							}
						}
					}
				}
			}
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;

		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double a = elements[i]->area() ;
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->grade(elements[i]->getState()) ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_2D)
					score += 1./(1.-sc);
				energy += elements[i]->getState().elasticEnergy()*a ;
			}
		}

		std::valarray<Point *> nularray(0) ;
		
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			elements[i]->setBoundingPoints(nularray) ;
			delete elements[i] ;
		}
	}
	else
	{
		Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
		
		std::vector<DelaunayTetrahedron *> elements ;
		std::vector<LinearForm *> behaviours ;
		TetrahedralElement father ;
		for(size_t i = 0 ; i < cache.size() ; i++)
		{
			DelaunayTetrahedron * tet = static_cast<DelaunayTetrahedron * > ((*mesh3d)[cache[i]]) ;
			elements.push_back(new DelaunayTetrahedron(NULL, NULL,   tet->first,  tet->second,   tet->third, tet->fourth,  NULL) );
			elements.back()->setBehaviour(tet->getBehaviour()->getCopy()) ;
			if(tet == s.getParent())
				elements.back()->getBehaviour()->setTensor(elements.back()->getBehaviour()->getTensor(elements.back()->getCenter())*(1.-delta_d)) ;
			elements.back()->getBehaviour()->getFractureCriterion()->setEnergyIndexed(false) ;
			elements.back()->refresh(&father);
			elements.back()->getState().initialize() ;
			K.add(elements.back());
			volume+= elements.back()->volume() ;
			
			for(size_t j = 0 ; j < tet->neighbour.size() ; j++)
			{
				if(tet->getNeighbour(j)->isTetrahedron() )
				{
					DelaunayTetrahedron * tri = dynamic_cast<DelaunayTetrahedron *>(tet->getNeighbour(j)) ;
					if(tri->getBehaviour()->type != VOID_BEHAVIOUR)
					{
						for(size_t k = 0 ; k <  tri->getBoundingPoints().size() ; k++)
						{
							if(!c.in(tri->getBoundingPoint(k)))
							{
								int id = tri->getBoundingPoint(k).id ;
								double ex = tri->getState().getDisplacements()[k*3];
								double ey = tri->getState().getDisplacements()[k*3+1];
								double ez = tri->getState().getDisplacements()[k*3+2];
								K.setPoint(ex, ey, ez ,id);
							}
						}
					}
				}
			}
		
			
		}
	
		K.cgsolve() ;
		for(size_t i = 0 ; i < elements.size() ; i++)
			elements[i]->step(0., &K.getDisplacements()) ;
		
		double energy = 0 ;
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			if(elements[i]->getBehaviour()->getFractureCriterion())
			{
				double sc = elements[i]->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
				if(dist(elements[i]->getCenter(), s.getParent()->getCenter()) < POINT_TOLERANCE_3D)
					score += 1./(1.-sc) ;
				energy += elements[i]->getState().elasticEnergy()*elements[i]->volume() ;
			}
		}
		for(size_t i = 0 ; i < elements.size() ;i++)
		{
			delete elements[i] ;
		}
	}
	
	
	return std::make_pair( (energy-originalenergy)/(delta_d),(score-originalscore)/(delta_d)) ;
}

std::pair<double, double> FractureCriterion::setChange(const ElementState &s) 
{
	stable = true ;
	
	if( !s.getParent()->getBehaviour()->getDamageModel())
		return std::make_pair(0.,0.) ;
	if(cache.size() == 0)
	{
		initialiseCache(s);
	}
	
	if(mesh2d)
	{
		
		// outside of the checkpoints, we only care about the order of the elements in 
		// term of their score. At the checkpoint, we consider the elements which
		// have met their criterion
		if(checkpoint) //new iteration
		{
			inset = false ;
			damagingSet.clear();
			proximitySet.clear() ;
			
			std::vector<unsigned int> newSet ;
			std::set<unsigned int> newProximity ;
			std::multimap<double, DelaunayTriangle *> sortedElements ;
			double thresholdScore = 0 ;
			
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;
				if(ci->getBehaviour()->getFractureCriterion() && !ci->getBehaviour()->fractured())
				{
					double renormScore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
					sortedElements.insert( std::make_pair(-renormScore, ci)) ;
				}
			}
			if(!sortedElements.empty())
			{
				thresholdScore = -sortedElements.begin()->first ;
			}
			if(thresholdScore > 0 && s.getParent()->getState().getDeltaTime() > POINT_TOLERANCE_2D)
				initialScore = std::max(1.+thresholdScore, POINT_TOLERANCE_2D) ;
			double minscore = thresholdScore ;
			double maxscore = 0 ;
			bool foundmaxscore = false ;
			minDeltaInNeighbourhood = 1 ;
			maxModeInNeighbourhood = -1 ;
			maxAngleShiftInNeighbourhood = 0 ;
			if(!sortedElements.empty() && thresholdScore > 0 )
			{
				for(auto i = sortedElements.begin() ; i != sortedElements.end() ; i++ )
				{
					if(std::abs(-i->first-thresholdScore) <= scoreTolerance*initialScore && -i->first > 0)
					{
						if(i->second == s.getParent() && met())
							inset = true ;
						if(i->second->getBehaviour()->getDamageModel()->getDelta() > POINT_TOLERANCE_2D)
							minDeltaInNeighbourhood = std::min(minDeltaInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getDelta()) ;
						maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getMode()) ;
						maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, i->second->getBehaviour()->getDamageModel()->getAngleShift()) ;
						newSet.push_back(i->second->index);
						minscore = std::min(-i->first, thresholdScore) ;
					}
					else
					{
						newProximity.insert(i->second->index) ;
						if(!foundmaxscore)
						{
							maxscore = -i->first ;
							foundmaxscore = true ;
						}
					}
				}
			}
			maxScoreInNeighbourhood = minscore ;
// 			std::cout << newSet.size() << std::endl ;
			
			if(!inset)
			{
				damagingSet.clear();
				proximitySet.clear() ;
				return std::make_pair(0.,0.) ;
			}
			
			if(!newSet.empty())
				std::stable_sort(newSet.begin(), newSet.end());
			
			damagingSet = newSet ;
			proximitySet.insert(proximitySet.end(), newProximity.begin(), newProximity.end()) ;
			
			
			return std::make_pair(minscore - maxscore /*- scoreTolerance*2.*initialScore*/, thresholdScore - minscore /*+ scoreTolerance*initialScore*/) ;
		}
		else if (inset)
		{

			checkpoint = false ;
			DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[damagingSet[0]]) ;
			double maxscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
			maxModeInNeighbourhood = ci->getBehaviour()->getDamageModel()->getMode() ;
			maxAngleShiftInNeighbourhood = ci->getBehaviour()->getDamageModel()->getAngleShift() ;
			
			double thresholdScore = maxscore ;
			for(size_t i = 1 ; i < damagingSet.size() ; i++)
			{
				ci = static_cast<DelaunayTriangle *>((*mesh2d)[damagingSet[i]]) ;
				double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
				maxModeInNeighbourhood = std::max(maxModeInNeighbourhood, ci->getBehaviour()->getDamageModel()->getMode()) ;
				maxAngleShiftInNeighbourhood = std::max(maxAngleShiftInNeighbourhood, ci->getBehaviour()->getDamageModel()->getAngleShift()) ;
				maxscore = std::min(maxscore,nls) ;
				thresholdScore = std::max(thresholdScore,nls) ;

			}
			maxScoreInNeighbourhood = maxscore ;

			double minscore = 0 ;
			if(!proximitySet.empty())
			{
				ci = static_cast<DelaunayTriangle *>((*mesh2d)[proximitySet[0]]) ;
				minscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
				
				for(size_t i = 1 ; i < proximitySet.size() ; i++)
				{
					ci = static_cast<DelaunayTriangle *>((*mesh2d)[proximitySet[i]]) ;
					double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
					
					minscore = std::max(nls, minscore) ;
				}
			}
			return std::make_pair(maxscore - minscore /*- scoreTolerance*2.*initialScore*/, thresholdScore - maxscore /*+ scoreTolerance*initialScore*/) ;
		}
	}
	else
	{
		
		// outside of the checkpoints, we only care about the order of the elements in 
		// term of their score. At the checkpoint, we consider the elements which
		// have met their criterion
		if(checkpoint) //new iteration
		{
			inset = false ;
			damagingSet.clear();
			proximitySet.clear() ;
			
			std::vector<unsigned int> newSet ;
			std::set<unsigned int> newProximity ;
			std::multimap<double, DelaunayTetrahedron *> sortedElements ;
			double thresholdScore = 0 ;
			
			for(size_t i = 0 ; i< cache.size() ; i++)
			{
				DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
				if(ci->getBehaviour()->getFractureCriterion())
				{
					double renormScore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
					sortedElements.insert( std::make_pair(-renormScore, ci)) ;
				}
			}
			
			
			if(!sortedElements.empty())
				thresholdScore = -sortedElements.begin()->first ;
			if(thresholdScore < 0 && s.getParent()->getState().getDeltaTime() > POINT_TOLERANCE_2D)
				initialScore = std::max(1.+thresholdScore, POINT_TOLERANCE_2D) ;
			
			double minscore = thresholdScore ;
			double maxscore = 0 ;
			bool foundmaxscore = false ;
			if(!sortedElements.empty() && thresholdScore > 0 )
			{
				for(auto i = sortedElements.begin() ; i != sortedElements.end() ; i++ )
				{
					if(-i->first >= thresholdScore - scoreTolerance*initialScore)
					{
						if(i->second == s.getParent() && met())
							inset = true ;
						newSet.push_back(i->second->index);
						minscore = -i->first ;
					}
					else
					{
						newProximity.insert(i->second->index) ;
						if(!foundmaxscore)
						{
							maxscore =  -i->first ;
							foundmaxscore = true ;
						}
					}
				}
			}

// 			std::cout << newSet.size() << std::endl ;
			
			if(!inset)
			{
				damagingSet.clear();
				proximitySet.clear() ;
				return std::make_pair(0.,0.) ;
			}
			
			if(!newSet.empty())
				std::stable_sort(newSet.begin(), newSet.end());
			
			damagingSet = newSet ;
			proximitySet.insert(proximitySet.end(), newProximity.begin(), newProximity.end()) ;
			
			
			return std::make_pair(minscore - maxscore - scoreTolerance*2.*initialScore, thresholdScore - minscore + scoreTolerance*initialScore) ;
		}
		else if (inset)
		{

			checkpoint = false ;
			DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[damagingSet[0]]) ;
			double maxscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

			double thresholdScore = maxscore ;
			for(size_t i = 1 ; i < damagingSet.size() ; i++)
			{
				ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[damagingSet[i]]) ;
				double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

				if(nls < maxscore)
				{
					maxscore = nls ;
				}
				if(nls > thresholdScore)
				{
					thresholdScore = nls ;
				}
			}
			double minscore = 0 ;
			if(!proximitySet.empty())
			{
				ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[proximitySet[0]]) ;
				minscore = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;
				
				for(size_t i = 1 ; i < proximitySet.size() ; i++)
				{
					ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[proximitySet[i]]) ;
					double nls = ci->getBehaviour()->getFractureCriterion()->nonLocalScoreAtState ;

					if(nls > minscore)
						minscore = nls ;
				}
			}
			return std::make_pair(maxscore - minscore - scoreTolerance*2.*initialScore, thresholdScore - maxscore + scoreTolerance*initialScore) ;
		}
	}
	
	//shut up the compiler
	return std::make_pair(0., 0.) ;
}

void FractureCriterion::step(ElementState &s)
{
	if(cache.empty() )
				initialiseCache(s) ;

	if(energyIndexed && s.getDeltaTime() > POINT_TOLERANCE_2D )
		noEnergyUpdate = true ;
	
	if(energyIndexed && noEnergyUpdate /*&& scoreAtState > -.5*/)
	{
		noEnergyUpdate = false ;
		currentEnergy = 0 ;
		bool buildvolume = false ;
		if(neighbourhoodvolume < 0)
		{
			neighbourhoodvolume = 0 ;
			buildvolume = true ;
		}
	
		if(mesh2d)
		{
			Circle c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			for(size_t i = 0 ; i < getCache().size() ; i++)
			{
				DelaunayTriangle *tri = static_cast<DelaunayTriangle *>( (*mesh2d)[getCache()[i]] ) ;
				if(c.in(*tri->first) || c.in(*tri->second) || c.in(*tri->third))
				{
					currentEnergy += tri->getState().elasticEnergy()*tri->area() ;
					if(buildvolume)
						neighbourhoodvolume += tri->area() ;
				}
			}
		}
		else if(mesh3d)
		{
			Sphere c(getMaterialCharacteristicRadius(), s.getParent()->getCenter()) ;
			
			for(size_t i = 0 ; i < getCache().size() ; i++)
			{
				DelaunayTetrahedron *tet = static_cast<DelaunayTetrahedron *>( (*mesh3d)[getCache()[i]] ) ;
				if(c.in(*tet->first) || c.in(*tet->second) || c.in(*tet->third) || c.in(*tet->fourth))
				{
					currentEnergy += tet->getState().elasticEnergy()*tet->volume() ;
					if(buildvolume)
						neighbourhoodvolume += tet->volume() ;
				}
			}
		}
		deltaEnergyAtState = (currentEnergy-previousEnergy) ;
		
		previousEnergy = currentEnergy ;
		
	}
	
	if(energyIndexed /*&& scoreAtState > -.5*/)
	{
		std::pair<double, double> dedc = getDeltaEnergyDeltaCriterion(s, 0.0001) ;
		energyDamageDifferential = dedc.first ;
		criterionDamageDifferential = dedc.second ;
	}
	
	if(s.getParent()->getBehaviour()->fractured())
	{
		scoreAtState = -1 ;
		return ;
	}
	scoreAtState = grade(s) ;
// 	if(!directionMet(0) && !directionMet(1) && !directionMet(2))
// 		scoreAtState = -1 ;
}

void FractureCriterion::computeNonLocalState(ElementState &s, NonLocalSmoothingType st)
{
	metAtStep = false ;

// 	if( s.getParent()->getBehaviour()->getDamageModel() && s.getParent()->getBehaviour()->getDamageModel()->fractured())
// 	{
// 		metAtStep = 0 ; scoreAtState > 0 ;
// 		nonLocalScoreAtState = -1 ; scoreAtState ;
// 		return  ;
// 	}
	
	DelaunayTriangle * testedTri = dynamic_cast<DelaunayTriangle *>(s.getParent()) ;
	DelaunayTetrahedron * testedTet = dynamic_cast<DelaunayTetrahedron *>(s.getParent()) ;
	HexahedralElement * testedHex = dynamic_cast<HexahedralElement *>(s.getParent()) ;
	
	switch (st)
	{
		case NULL_SMOOTH :
		{
			metAtStep = scoreAtState > 0 ;
			nonLocalScoreAtState = scoreAtState ;
			return ;
		}
		case MAX_PROXIMITY_SMOOTH :
		{
			if(testedTri)
			{
				if(scoreAtState <= 2.*scoreTolerance)
				{
					nonLocalScoreAtState = scoreAtState ;
					metAtStep = false ;
					return ;
				}
				
				if(!cache.empty())
				{
					nonLocalScoreAtState = scoreAtState ;
					for(size_t i = 0 ; i< cache.size() ; i++)
					{
						DelaunayTriangle * ci = static_cast<DelaunayTriangle *>((*mesh2d)[cache[i]]) ;

						if(squareDist2D(ci->getCenter(), s.getParent()->getCenter()) > physicalCharacteristicRadius*physicalCharacteristicRadius
							|| ci->getBehaviour()->getSource() !=  s.getParent()->getBehaviour()->getSource()
							|| ci->getBehaviour()->fractured()
						)
							continue ;
							
							nonLocalScoreAtState = std::max(nonLocalScoreAtState, ci->getBehaviour()->getFractureCriterion()->scoreAtState) ;
						}
					}
					metAtStep = true ;
					return ;

			}
			
			if(testedTet)
			{
				nonLocalScoreAtState = scoreAtState ;
				if(testedTet->visited())
				{
					metAtStep = false ;
					return  ;
				}
						
				if (scoreAtState <= 0)
				{
					metAtStep = false ;
					return  ;
				}

				double maxNeighbourhoodScore = 0 ;
// 				double matchedArea = 0 ;
// 				std::map<double, DelaunayTetrahedron *> scores ;
// 				std::map<DelaunayTetrahedron *, double> areatemp ;
// 				DelaunayTetrahedron * maxLocus = NULL;
				
				if(!cache.empty())
				{
					for(size_t i = 0 ; i< cache.size() ; i++)
					{
						DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;

						if( !ci->getBehaviour()->fractured())
						{
							double s = ci->getBehaviour()->getFractureCriterion()->getScoreAtState() ;
// 							scores[-s] =  ci;
							if(s > maxNeighbourhoodScore)
							{
								maxNeighbourhoodScore = s ;
// 								maxLocus = ci ;
							}
						}
// 						else if(ci->getBehaviour()->fractured())
// 						{
// 							double s = POINT_TOLERANCE_2D ;
// // 							scores[-s] =  ci;
// 						}
// 						areatemp[ci] = area[i] ;

					}
				}
				
				if(maxNeighbourhoodScore < 0)
				{
					metAtStep = false ;
					return  ;
				}
				
				std::vector<DelaunayTetrahedron *> maxloci ;
				
				for(size_t i = 0 ; i< cache.size() ; i++)
				{
					DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
					if(ci->getBehaviour()->getFractureCriterion())
						if(std::abs(ci->getBehaviour()->getFractureCriterion()->getScoreAtState()-maxNeighbourhoodScore) < scoreTolerance)
							maxloci.push_back(ci) ;
				}


				for(size_t i = 0 ; i < maxloci.size() ; i++)
					if(squareDist3D(maxloci[i]->getCenter(), s.getParent()->getCenter()) < physicalCharacteristicRadius*physicalCharacteristicRadius)
				{
					metAtStep = true ;
					return  ;
				}
				metAtStep = false ;
				return  ;
			}
			else if(testedHex)
			{
				std::set<HexahedralElement *> neighbourhood ;
				std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
				for(size_t i = 0 ; i < neighbours.size() ; i++)
				{
					for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
					{
						if(neighbours[i]->neighbourhood[j] != testedHex 
							&& !neighbours[i]->neighbourhood[j]->getBehaviour()->fractured())
							neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
					}
				}
				double score = grade(s) ;
				double maxNeighbourhoodScore = 0 ;
				if(!neighbourhood.empty())
				{
					for(auto i= neighbourhood.begin() ; i != neighbourhood.end() ; ++i)
					{
						if((*i)->getBehaviour()->getFractureCriterion() 
							&& !(*i)->getBehaviour()->fractured())
							maxNeighbourhoodScore = std::max(maxNeighbourhoodScore,
																							(*i)->getBehaviour()->getFractureCriterion()->grade((*i)->getState())) ;
						if((*i)->getBehaviour()->changed())
						{
							maxNeighbourhoodScore = 10.*score ;
							break ;
						}
						
						if (maxNeighbourhoodScore > score)
							break ;
						
					}
				}
				
				if( score > 0 )
				{
					if(score > maxNeighbourhoodScore)
					{
						metAtStep = true ;
						return  ;
					}
				}
				metAtStep = false ;
				return  ;
			}
			else
			{
				std::cout << " criterion not implemented for this kind of element" << std::endl ;
				metAtStep = false ;
				return  ;
			}
		}
		case GAUSSIAN_SMOOTH :
		{
			if(testedTri)
			{
				if(testedTri->visited)
				{
					metAtStep = false ;
					return  ;
				}
				
				double smoothscore = smoothedScore(s) ;
				
				metAtStep =  (smoothscore > scoreTolerance) ;
				if(scoreAtState >= 0)
				  std::cout << scoreAtState << "  " << smoothscore << std::endl ;
				nonLocalScoreAtState = smoothscore ;
				return ;

			}
			if(testedTet)
			{

				if(testedTet->visited())
				{
					metAtStep = false ;
					return  ;
				}
				
				double str = 0 ;
				double fact = 0 ;
				
				for(size_t i = 0 ; i< cache.size() ; i++)
				{
					DelaunayTetrahedron * ci = static_cast<DelaunayTetrahedron *>((*mesh3d)[cache[i]]) ;
					double dc = squareDist3D(s.getParent()->getCenter(), ci->getCenter()) ;
					if(ci->getBehaviour()->getFractureCriterion())
					{
						double d =  exp(-dc/(physicalCharacteristicRadius*physicalCharacteristicRadius) );
						double a = ci->volume() ;
						str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
						fact+=d*a ;
						if(mirroring == MIRROR_X && std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Y &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Z &&  std::abs(ci->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs(ci->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs(ci->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs(ci->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs(ci->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=ci->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
					}
				}
				
				double smoothscore = str/fact ;
				metAtStep =  smoothscore > 0 ;
				nonLocalScoreAtState = smoothscore ;
				return ;
			}
			else if(testedHex)
			{
				std::set<HexahedralElement *> neighbourhood ;
				std::vector<HexahedralElement *> neighbours = testedHex->neighbourhood ;
				for(size_t i = 0 ; i < neighbours.size() ; i++)
				{
					for(size_t j = 0 ; j <  neighbours[i]->neighbourhood.size() ; j++)
					{
							neighbourhood.insert(neighbours[i]->neighbourhood[j]) ;
					}
				}
				double str = 0 ;
				double fact = 0 ;
				
				for(auto ci = neighbourhood.begin() ; ci != neighbourhood.end() ; ++ci)
				{
					double dc = squareDist3D(s.getParent()->getCenter(), (*ci)->getCenter()) ;
					if((*ci)->getBehaviour()->getFractureCriterion())
					{
						double d =  exp(-dc/(physicalCharacteristicRadius*physicalCharacteristicRadius) );
						double a = (*ci)->volume() ;
						str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
						fact+=d*a ;
						if(mirroring == MIRROR_X && std::abs((*ci)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_X
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Y &&  std::abs((*ci)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_Z &&  std::abs((*ci)->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_Y
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs((*ci)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XY &&  std::abs((*ci)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs((*ci)->getCenter().x  - delta_x) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_XZ &&  std::abs((*ci)->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs((*ci)->getCenter().y  - delta_y) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
						if(mirroring == MIRROR_YZ &&  std::abs((*ci)->getCenter().z  - delta_z) < physicalCharacteristicRadius) // MIRROR_XY
						{
							str +=(*ci)->getBehaviour()->getFractureCriterion()->getScoreAtState() *a*d ;
							fact+=d*a ;
						}
					}
				}
				
				double smoothscore = str/fact ;
				metAtStep =  smoothscore > 0 ;
				nonLocalScoreAtState = smoothscore ;
				return ;
			}
			else
			{
				std::cout << " criterion not implemented for this kind of element" << std::endl ;
				metAtStep = false ;
				return  ;
			}
		}
	}
	//shut up the compiler
		metAtStep = false ;
		return  ;
}

FractureCriterion::~FractureCriterion()
{
}

void FractureCriterion::setNeighbourhoodRadius(double r)
{
	neighbourhoodradius = r ;
	cache.clear() ;
	factors.clear();
}

void FractureCriterion::setMaterialCharacteristicRadius(double r)
{
	physicalCharacteristicRadius = r ;
	neighbourhoodradius = std::max(r*2.5, neighbourhoodradius) ;
	cache.clear() ;
	factors.clear();
}

bool FractureCriterion::met() const
{
	return metAtStep ;
}

Material FractureCriterion::toMaterial()
{
	Material mat ;
	return mat ;
}


