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

namespace Amie {

/** \brief Mohr - Coulomb fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The Mohr-Coulomb is met when the maximum principal stresses is below or above the prescribed limits 
	
*/

/*PARSE MohrCoulomb FractureCriterion 
    @value[tensile_strength] // maximum stress in tension (positive)
    @value[compressive_strength] // maximum stress in compression (negative)
*/
class MohrCoulomb : public FractureCriterion
{
public:
	double upVal ;
	double downVal ;
	bool metInCompression  ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction, double t = 0) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction, double t = 0) {return metInTension ;}
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	MohrCoulomb(double up, double down);

	virtual ~MohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;
		
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};

/*PARSE NonLocalMohrCoulomb FractureCriterion 
    @value[tensile_strength] // maximum stress in tension (positive)
    @value[compressive_strength] // maximum stress in compression (negative)
    @value[young_modulus] // Young modulus of the material
*/
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
	
	virtual bool directionInTension(size_t direction, double t = 0) {return true ;}
	virtual bool directionInCompression(size_t direction, double t = 0) {return true ;}
	virtual bool directionMet(size_t direction, double t = 0) 
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
	NonLocalMohrCoulomb(double up, double down, double E);

	virtual ~NonLocalMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual void scale(double f) { upVal *=f ; downVal *= f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};


/*PARSE SpaceTimeNonLocalMohrCoulomb FractureCriterion 
    @value[tensile_strength] // maximum stress in tension (positive)
    @value[compressive_strength] // maximum stress in compression (negative)
    @value[young_modulus] // Young modulus of the material
*/
class SpaceTimeNonLocalMohrCoulomb : public NonLocalMohrCoulomb
{
protected:
	PointArray testPoints ;
public:
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
 * @param E young's Modulus
*/
	SpaceTimeNonLocalMohrCoulomb(double up, double down, double E) : NonLocalMohrCoulomb(up, down, E) { } ;

	virtual ~SpaceTimeNonLocalMohrCoulomb() { } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalMohrCoulomb(*this) ; }

	virtual double grade(ElementState &s)  ;
    virtual double gradeAtTime(ElementState &s, double t)  ;

};


/*PARSE NonLocalLinearlyDecreasingMohrCoulomb FractureCriterion 
    @value[tensile_strength] // maximum stress in tension (positive)
    @value[compressive_strength] // maximum stress in compression (negative)
    @value[tensile_ultimate_strain] // strain in tension at the end of the softening (positive)
    @value[compressive_ultimate_strain] // strain in compression at the end of the softening (negative)
    @value[young_modulus] // Young modulus of the material
*/
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
	
	virtual bool directionInTension(size_t direction, double t = 0) {return !metInCompression ;}
	virtual bool directionInCompression(size_t direction, double t = 0) {return !metInTension ;}
	virtual bool directionMet(size_t direction, double t = 0) 
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
	NonLocalLinearlyDecreasingMohrCoulomb(double up, double down, double limittstrain, double limitcstrain,  double E);

	virtual ~NonLocalLinearlyDecreasingMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual void scale(double f) { upVal *=f ; downVal *= f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};

/*PARSE NonLocalExponentiallyDecreasingMohrCoulomb FractureCriterion 
    @value[tensile_strength] // maximum stress in tension (positive)
    @value[compressive_strength] // maximum stress in compression (negative)
    @value[tensile_ultimate_strain] // strain in tension at the end of the softening (positive)
    @value[compressive_ultimate_strain] // strain in compression at the end of the softening (negative)
    @value[young_modulus] // Young modulus of the material
*/
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
	
	virtual bool directionInTension(size_t direction, double t = 0) 
	{
		if(direction == 1)
			return false ;
		return true ;
		
	}
	virtual bool directionInCompression(size_t direction, double t = 0) 
	{
		return false ;
	}
	virtual bool directionMet(size_t direction, double t = 0) 
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
	NonLocalExponentiallyDecreasingMohrCoulomb(double up, double down, double limittstrain, double limitcstrain,  double E);

	virtual ~NonLocalExponentiallyDecreasingMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual void scale(double f) { upVal *=f ; downVal *= f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};

/*PARSE NonLocalInverseRootMohrCoulomb FractureCriterion 
    @value[tensile_strain] // maximum strain in tension (positive)
    @value[tensile_ultimate_strain] // strain in tension at the end of the softening (positive)
    @value[young_modulus] // Young modulus of the material
    @value[tensile_strain_coefficient] // indicates how fast the softening occurs
*/
class NonLocalInverseRootMohrCoulomb : public FractureCriterion
{
public:
	double upVal ;
	double stiffness ;
	double limitstrain ;
	double limitystrain ;
	double c ;
	bool metInTension  ;
	
	virtual bool directionInTension(size_t direction, double t = 0) {return true ;}
	virtual bool directionInCompression(size_t direction, double t = 0) {return false ;}
	virtual bool directionMet(size_t direction, double t = 0) 
	{
		if(direction == 0)
			return true ;
		
		return false ;
	}
	
/** \brief Constructor, set the maximum and minimum strain
 * @param limitstrain Maximum strain (tension) for elastic behaviour
 * @param limitystrain maximum strain before failure
 * @param E initial Young's modulus
 * @param c in 1/(sqrt(c epsilon))
*/
	NonLocalInverseRootMohrCoulomb(double limitstrain, double limitystrain,  double E, double c);

	virtual ~NonLocalInverseRootMohrCoulomb();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(ElementState &s)  ;

	virtual void scale(double f) { upVal *=f ; } ;
	
	virtual double getTensileLimit(const ElementState & s) const {return upVal ;};
};


}

#endif
