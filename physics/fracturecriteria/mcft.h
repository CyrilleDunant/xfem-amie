//
// C++ Interface: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MCFT_H
#define MCFT_H

#include "fracturecriterion.h"

namespace Mu {

typedef enum
{
	UPPER_BOUND,
	LOWER_BOUND,
	AVERAGE
} RedistributionType ;

class NonLocalMCFT : public FractureCriterion
{

	double getBareConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress) ;
	double getRebarConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress) ;
	double getConcreteTensileCriterion(const ElementState & s, double pseudoYoung, double tstrain, double tstress) ;
	double getConcreteCompressiveCriterion(const ElementState & s, double pseudoYoung, double cstrain, double tstress, double cstress) ;
	void initialise() ;
	RedistributionType rtype ;
public:
	bool initialised ;
	bool firstTension ;
	bool firstCompression ;
	bool secondTension ;
	bool secondCompression ;
	bool firstMet ;
	bool secondMet ;
	double upVal ;
	double downVal ;
	double tensionCritStrain ;
	double critStrain ;
	double youngModulus ;
	double k ;
	double strain_ch ;
	double strain_te ;

	double scaleFactor ;
	virtual bool directionInTension(size_t direction) 
	{
		if(direction == 0)
			return firstTension ;
		if(direction == 1)
			return secondTension ;
		return false ;
	}
	virtual bool directionInCompression(size_t direction) 
	{
		if(direction == 0)
			return firstCompression ;
		if(direction == 1)
			return secondCompression ;
		return false ;
	}
	virtual bool directionMet(size_t direction) 
	{
		if(direction == 0)
			return firstMet ;
		if(direction == 1)
			return secondMet ;
		return false ;
	}
	std::vector<std::pair<double, double> > rebarLocationsAndDiameters ;
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
NonLocalMCFT(double up, double down, double youngModulus, double charDistance, RedistributionType r = UPPER_BOUND, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;

	virtual ~NonLocalMCFT();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual Material toMaterial() ;
	
	virtual void scale(double d) 
	{
		scaleFactor = d ;
// 		upVal *= d ; 
// 		downVal *= d ;
// 		youngModulus *=d ;
		initialised = false ;
// 		tensionCritStrain /= d ;
// 		critStrain /= d ;
	};
	
	virtual double getTensileLimit(const ElementState & s) const ;
};

}

#endif
