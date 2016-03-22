//
// C++ Interface: mazars
//
// Description: 
//
//
// Author: Adrien Hilaire <adrien.hilairet@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MAZARS_H
#define MAZARS_H

#include "fracturecriterion.h"

namespace Amie {

/*PARSE NonLocalMazars FractureCriterion --has-reset
    @value[tensile_strain] // maximum strain in tension (positive)
    @value[young_modulus] // Young modulus of the material
    @value[poisson_ratio] // Poisson ratio of the material
    @value[tensile_fracture_energy] // fracture energy in tension
    @value[compressive_strain] 1 // maximum strain in compression (negative; a positive value disables the compressive part of the criterion)
    @value[compressive_stress] 1 // maximum stress in compression (negative; a positive value disables the compressive part of the criterion)
    @value[material_characteristic_radius] // characteristic length of the non-local approximation
    @string<planeType>[plane_type] PLANE_STRESS // planar approximation in 2 dimensions 
*/

	/** \brief The Mazars fracture criterion is met when the equivalent Mazars strain reaches a threshold level
	      @author Adrien Hilaire <adrien.hilaire@epfl.ch>
	*/
	class NonLocalMazars : public FractureCriterion
	{
    protected:
		bool ismet ;
		bool tensionOnly ;
		double B_t;
		double B_c ;
	public:
	/** \brief  
	 *   @param threshold Maximum strain (tension)
	 *   @param E Young modulus
	 *   @param Gf Fracture energy (J.m^-2)
	 *   @param nu Poisson coefficient
	 *   @param cstrain strain at stress peak (compression)
	 *   @param cstress peak stress  (compression)
	 */ 
		double threshold ;
		double E ;
		double Gf ;
		double nu ; 
		double A_c ; 
		double cstrain; 
		double cstress;
		planeType pt;
		virtual bool directionInTension(size_t direction, double t = 0) {return ismet ;}
		virtual bool directionInCompression(size_t direction, double t = 0) {return ismet ;}
		virtual bool directionMet(size_t direction, double t = 0) {return ismet;}
	public:
	/** \brief Constructor 
	 * @param thres Set the maximum equivalent Mazars strain. 
	 */
		NonLocalMazars(double thres, double E, double nu, double Gf, double cstress, double cstrain, double radius, planeType pt, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

		NonLocalMazars(double thres, double E, double nu, double Gf, double radius, planeType pt, MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);

		virtual ~NonLocalMazars();

	/** \brief Return a copy of this criterion
	 */
		virtual FractureCriterion * getCopy() const;

	/** \brief Return normalised distance to the fracture surface
	 *
	 * The distance is computed as: \f$ 1.-|\frac{Limit\; stress}{max\; principal\; strain\; in\; element}|  \f$
	 * @param s ElementState to consider
	*/
		virtual double grade(ElementState &s)  ;
		double gradeAtTime(ElementState &s, double t) ;
		
		virtual double getTensileLimit(const ElementState & s) const {return threshold ;};

		virtual void reset( double thr, double E, double nu, double Gf, double cstress, double cstrain, double r = -1, planeType = PLANE_STRESS) ;

		virtual void reset( double thr, double E, double nu, double Gf) ;
	

	};
	

/*PARSE NonLocalSpaceTimeMazars FractureCriterion --has-reset
    @value[tensile_strain] // maximum strain in tension (positive)
    @value[young_modulus] // Young modulus of the material
    @value[poisson_ratio] // Poisson ratio of the material
    @value[tensile_fracture_energy] // fracture energy in tension
    @value[compressive_strain] 1 // maximum strain in compression (negative; a positive value disables the compressive part of the criterion)
    @value[compressive_stress] 1 // maximum stress in compression (negative; a positive value disables the compressive part of the criterion)
    @value[material_characteristic_radius] // characteristic length of the non-local approximation
    @string<planeType>[plane_type] PLANE_STRESS // planar approximation in 2 dimensions 
*/
class NonLocalSpaceTimeMazars : public NonLocalMazars
{
public:
	NonLocalSpaceTimeMazars(double thres, double E, double nu, double Gf, double cstress, double cstrain, double radius, planeType pt,  MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);
	NonLocalSpaceTimeMazars(double thres, double E, double nu, double Gf, double radius, planeType pt,  MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);
	virtual ~NonLocalSpaceTimeMazars();
	virtual double grade(ElementState &s)  ;
	virtual FractureCriterion * getCopy() const;
} ;

} 

#endif
