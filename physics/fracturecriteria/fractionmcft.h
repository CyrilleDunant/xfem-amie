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
#ifndef FRAC_MCFT_H
#define FRAC_MCFT_H

#include "fracturecriterion.h"

namespace Amie {

/** \brief Mohr - Coulomb fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The MCFT fraction is designed to model failing concrete around steel 
	
*/
// class FractionMCFT : public FractureCriterion
// {
// 
// public:
// 	double upVal ;
// 	double downVal ;
// 	Matrix steelCGTensor ;
// 	double phi ;
// 	double tensionCritStrain ;
// 	bool metInCompression  ;
// 	bool metInTension  ;
// 	
// 	virtual bool directionInTension(size_t direction) {return metInCompression ;}
// 	virtual bool directionInCompression(size_t direction) {return metInTension ;}
// 	
// /** \brief Constructor, set the maximum and minimum strain
//  * @param up Maximum stress (tension)
//  * @param down Minimum stress (compression)
//  * @param concreteCGTensor stiffness tensor of the concrete
// */
// 	FractionMCFT(double up, double down, Matrix concreteCGTensor, double youngModulus, double phi,  MirrorState mirroring = NO_MIRROR, double delta_x = 0, double delta_y = 0, double delta_z = 0);
// 
// 	virtual ~FractionMCFT();
// 
// /** \brief Return a copy of this fracture criterion*/
// 	virtual FractureCriterion * getCopy() const;
// 
// /** \brief return the normalised distance to the fracture surface
//  *
//  * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
//  * @param s ElementState to consider
// */
// 	virtual double grade(ElementState &s)  ;
// 	
// 	virtual void scale(double d) {upVal *= d ; downVal *=d ;};
// 	
// 	virtual double getTensileLimit(const ElementState & s) const {return upVal ; } ;
// 
// 	virtual Material toMaterial() ;
// };

}

#endif
