//
// C++ Interface: mohrcoulomb
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef FRAC_MCFT_H
#define FRAC_MCFT_H

#include "fracturecriterion.h"

namespace Mu {

/** \brief Mohr - Coulomb fracture criterion
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	The MCFT fraction is designed to model failing concrete around steel 
	
*/
class FractionMCFT : public FractureCriterion
{

public:
	double upVal ;
	double downVal ;
	Matrix concreteCGTensor ;
/** \brief Constructor, set the maximum and minimum strain
 * @param up Maximum stress (tension)
 * @param down Minimum stress (compression)
 * @param concreteCGTensor stiffness tensor of the concrete
*/
	FractionMCFT(double up, double down, Matrix concreteCGTensor);

	virtual ~FractionMCFT();

/** \brief Return a copy of this fracture criterion*/
	virtual FractureCriterion * getCopy() const;

/** \brief return the normalised distance to the fracture surface
 *
 * The distance is computed as: \f$ 1.-|\frac{max\; principal\; strain\; in\; element}{Limit\; strain}|  \f$
 * @param s ElementState to consider
*/
	virtual double grade(const ElementState &s)  ;

	virtual Material toMaterial() ;
};

}

#endif
