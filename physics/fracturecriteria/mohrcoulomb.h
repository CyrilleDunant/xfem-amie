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
public:
	double upVal ;
	double downVal ;
	bool metInCompression  ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction) {return metInTension ;}
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
	virtual double grade(ElementState &s)  ;

	virtual Material toMaterial() ;
	
	virtual void scale(double f) { upVal *=f ; downVal *= f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};

class NonLocalMohrCoulomb : public FractureCriterion
{
protected:
	PointArray testPoints ;
public:
	double upVal ;
	double downVal ;
	double stiffness ;
	bool metInCompression  ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction) {return true ;}
	virtual bool directionInCompression(size_t direction) {return true ;}
	virtual bool directionMet(size_t direction) 
	{
		if(direction == 0)
			return metInTension ;
		if(direction == 1)
			return metInCompression ;
		
		return false ;
	}
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	NonLocalMohrCoulomb(double up, double down, double E, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

	virtual ~NonLocalMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual Material toMaterial() ;

	virtual void scale(double f) { upVal *=f ; downVal *= f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};

class NonLocalLinearlyDecreasingMohrCoulomb : public FractureCriterion
{
public:
	double upVal ;
	double downVal ;
	double stiffness ;
	double limittstrain ;
	double limitcstrain ;
	bool metInCompression  ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction) {return !metInCompression ;}
	virtual bool directionInCompression(size_t direction) {return !metInTension ;}
	virtual bool directionMet(size_t direction) 
	{
		if(direction == 0)
			return metInTension ;
		if(direction == 1)
			return metInCompression ;
		
		return false ;
	}
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	NonLocalLinearlyDecreasingMohrCoulomb(double up, double down, double limittstrain, double limitcstrain,  double E, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

	virtual ~NonLocalLinearlyDecreasingMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual Material toMaterial() ;

	virtual void scale(double f) { upVal *=f ; downVal *= f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};

class NonLocalExponentiallyDecreasingMohrCoulomb : public FractureCriterion
{
public:
	double upVal ;
	double downVal ;
	double stiffness ;
	double limittstrain ;
	double limitcstrain ;
	bool metInCompression  ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction) 
	{
		if(direction == 1)
			return false ;
		return true ;
		
	}
	virtual bool directionInCompression(size_t direction) 
	{
		return false ;
	}
	virtual bool directionMet(size_t direction) 
	{
		if(direction == 0)
			return metInTension ;
		if(direction == 1)
			return metInCompression ;
		
		return false ;
	}
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	NonLocalExponentiallyDecreasingMohrCoulomb(double up, double down, double limittstrain, double limitcstrain,  double E, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

	virtual ~NonLocalExponentiallyDecreasingMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual Material toMaterial() ;

	virtual void scale(double f) { upVal *=f ; downVal *= f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};

class NonLocalInverseRootMohrCoulomb : public FractureCriterion
{
public:
	double upVal ;
	double stiffness ;
	double limitstrain ;
	double limitystrain ;
	double c ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction) {return true ;}
	virtual bool directionInCompression(size_t direction) {return false ;}
	virtual bool directionMet(size_t direction) 
	{
		if(direction == 0)
			return true ;
		
		return false ;
	}
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	NonLocalInverseRootMohrCoulomb(double limitstrain, double limitystrain,  double E, double c, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

	virtual ~NonLocalInverseRootMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual Material toMaterial() ;

	virtual void scale(double f) { upVal *=f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};


}

#endif
