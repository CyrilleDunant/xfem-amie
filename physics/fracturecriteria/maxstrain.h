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
class MaximumStrain : public FractureCriterion
{
protected:
	bool metInCompression  ;
	bool metInTension  ;
public:
	double upVal ;
	
	virtual bool directionInTension(size_t direction) {return metInCompression ;}
	virtual bool directionInCompression(size_t direction) {return metInTension ;}

	/** \brief Constructor 
	 * @param up Set the maximum strain. 
	 */
	MaximumStrain(double up, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;

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

	virtual Material toMaterial() ;
		
	virtual double getTensileLimit(const ElementState & s) const {return upVal*20e9 ; } ;
};

class SpaceTimeNonLocalMaximumStrain : public MaximumStrain
{
protected:
	PointArray testPoints ;
public:
	double maxstress ;
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	SpaceTimeNonLocalMaximumStrain(double up, double mstr, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) : MaximumStrain(up, mirroring, delta_x, delta_y, delta_z),maxstress(mstr) { } ;

	virtual ~SpaceTimeNonLocalMaximumStrain() { } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalMaximumStrain(*this) ; }

	virtual double grade(ElementState &s)  ;

};

class SpaceTimeNonLocalLinearSofteningMaximumStrain : public SpaceTimeNonLocalMaximumStrain
{
protected:
	PointArray testPoints ;
public:
	double yieldstrain ;
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	SpaceTimeNonLocalLinearSofteningMaximumStrain(double up, double mstr, double lim, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) : SpaceTimeNonLocalMaximumStrain(up, mstr, mirroring, delta_x, delta_y, delta_z),yieldstrain(lim) { } ;

	virtual ~SpaceTimeNonLocalLinearSofteningMaximumStrain() { } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalLinearSofteningMaximumStrain(*this) ; }

	virtual double grade(ElementState &s)  ;

};

class SpaceTimeNonLocalMaximumStress : public MaximumStrain
{
protected:
	PointArray testPoints ;
public:
	double maxstress ;
	
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
*/
	SpaceTimeNonLocalMaximumStress(double up, double mstr, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) : MaximumStrain(up, mirroring, delta_x, delta_y, delta_z),maxstress(mstr) { } ;

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
	SpaceTimeNonLocalEllipsoidalMixedCriterion(double up, double mstr, double E0, double Einf, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0) ;

	virtual ~SpaceTimeNonLocalEllipsoidalMixedCriterion() { delete surface ; } ;

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const { return new SpaceTimeNonLocalEllipsoidalMixedCriterion(upVal, maxstress, E_inst, E_relaxed) ; }

	virtual double grade(ElementState &s)  ;

};


}

#endif
