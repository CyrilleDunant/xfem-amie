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
#ifndef MUMOHRCOULOMB_H
#define MUMOHRCOULOMB_H

#include "fracturecriterion.h"

namespace Mu {

/** \brief Mohr - Coulomb fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The Mohr-Coulomb is met when the maximum principal stresses is below or above the prescribed limits 
	
*/
class MohrCoulomb : public FractureCriterion
{
protected:
	PointArray testPoints ;
public:
	double upVal ;
	double downVal ;
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	MohrCoulomb(double up, double down, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

	virtual ~MohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(const ElementState &s)  ;

	virtual Material toMaterial() ;
	
	virtual void multiply(double f) { upVal *=f ; downVal *= f ; } ;
};

class NonLocalMohrCoulomb : public FractureCriterion
{
protected:
	PointArray testPoints ;
public:
	double upVal ;
	double downVal ;
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	NonLocalMohrCoulomb(double up, double down, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

	virtual ~NonLocalMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(const ElementState &s)  ;

	virtual Material toMaterial() ;

	virtual void multiply(double f) { upVal *=f ; downVal *= f ; } ;
};

}

#endif
