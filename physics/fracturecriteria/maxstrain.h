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
#ifndef MAX_STRAIN_H__
#define MAX_STRAIN_H__

#include "fracturecriterion.h"
#include "../../mesher/delaunay_3d.h"

namespace Amie {

/** \brief Maximum strain fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The maximum (tensile) strain criterion is met when a strain limit is reached.
	
*/
/*PARSE MaximumStrain FractureCriterion
    @value[tensile_strain] // maximum strain in tension
*/
class MaximumStrain : public FractureCriterion
{
protected:
	bool metInCompression  ;
	bool metInTension  ;
public:
	double upVal ;
	
	virtual bool directionInTension(size_t direction, double t = 0) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction, double t = 0) {return metInTension ;}

	/** \brief Constructor 
	 * @param up Set the maximum strain. 
	 */
	MaximumStrain(double up) ;

	virtual ~MaximumStrain();

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; strain}{max\; Mises; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
	virtual double grade(ElementState &s)  ;
 	
	/** \brief Return a copy of this criterion
	 */
	virtual FractureCriterion * getCopy() const;

	virtual double getTensileLimit(const ElementState & s) const {return upVal*20e9 ; } ;
};

/*PARSE SpaceTimeNonLocalMaximumStrain FractureCriterion
    @value[tensile_strain] // maximum strain in tension
*/
class SpaceTimeNonLocalMaximumStrain : public MaximumStrain
{
protected:
	PointArray testPoints ;
public:
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	SpaceTimeNonLocalMaximumStrain(double up) : MaximumStrain(up) { } ;

	virtual ~SpaceTimeNonLocalMaximumStrain() { } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalMaximumStrain(*this) ; }

	virtual double grade(ElementState &s)  ;

};

/*PARSE SpaceTimeNonLocalLinearSofteningMaximumStrain FractureCriterion --has-reset
    @value[tensile_strain] // strain at the peak in tension
    @value[tensile_strength] // stress at the peak in tension
    @value[tensile_ultimate_strain] // maximum strain in tension
*/
class SpaceTimeNonLocalLinearSofteningMaximumStrain : public SpaceTimeNonLocalMaximumStrain
{
protected:
	PointArray testPoints ;
public:
	double yieldstrain ;
	double maxstress ;
	
	SpaceTimeNonLocalLinearSofteningMaximumStrain(double up, double mstr, double lim) : SpaceTimeNonLocalMaximumStrain(up),yieldstrain(lim),maxstress(mstr) { } ;

	void reset(double up, double mstr, double lim) { upVal = up; maxstress = mstr; yieldstrain = lim ; }

	virtual ~SpaceTimeNonLocalLinearSofteningMaximumStrain() { } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalLinearSofteningMaximumStrain(*this) ; }

	virtual double grade(ElementState &s)  ;
        virtual double gradeAtTime(ElementState &s, double t)  ;

};

/*PARSE SpaceTimeNonLocalMaximumStress FractureCriterion
    @value[tensile_strength] // stress at the peak in tension
*/
class SpaceTimeNonLocalMaximumStress : public MaximumStrain
{
protected:
	PointArray testPoints ;
public:
	double maxstress ;
	
/** \brief Constructor, set the maximum strain
 * @param mstr Maximum stress (tension)
*/
	SpaceTimeNonLocalMaximumStress(double mstr) : MaximumStrain(mstr),maxstress(mstr) { } ;

	virtual ~SpaceTimeNonLocalMaximumStress() { } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalMaximumStress(*this) ; }

	virtual double grade(ElementState &s)  ;

};

class SpaceTimeNonLocalEllipsoidalMixedCriterion : public MaximumStrain
{
protected:
	PointArray testPoints ;
public:
	double maxstress ;
	double E_inst ;
	double E_relaxed ;
	Ellipse* surface ;
	double base ;
	double renormStrain ;
	double renormStress ;
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	SpaceTimeNonLocalEllipsoidalMixedCriterion(double up, double mstr, double E0, double Einf) ;

	virtual ~SpaceTimeNonLocalEllipsoidalMixedCriterion() { delete surface ; } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalEllipsoidalMixedCriterion(upVal, maxstress, E_inst, E_relaxed) ; }

	virtual double grade(ElementState &s)  ;

};


}

#endif
