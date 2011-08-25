//
// C++ Interface: kelvinvoight
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __KELVIN_VOIGHT_H_
#define __KELVIN_VOIGHT_H_

#include "physics_base.h"

namespace Mu
{


/** \brief A Kelvin-Voight Law
* The field param is the Cauchy-Green Strain Tensor
* The viscosity tensor is also stored
*/

struct KelvinVoight : public LinearForm
{
	Matrix eta ;
	std::vector<Variable> v ;
	/** \brief Constructor
	*
	* @param rig Complete expression of the Cauchy-Green Strain Tensor
	* @param eta Complete expression of the viscosity Tensor
	*/
	KelvinVoight( const Matrix & rig, const Matrix & eta ) ;

	virtual ~KelvinVoight() ;

	/** \brief Apply the law.
	 *
	 * The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j + \dot{\nabla}^T h_i E \dot{\nabla} h_j\f$
	 * @param p_i first basis polynomial.
	 * @param p_j second basis polynomial.
	 * @param gp Gauss Points used for the quadrature
	 * @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
	 * @param ret Matrix to store the result
	 * @param vm virtualMachine to use to compute the result
	 */
	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	/** \brief is Element fractured
	 *
	 * @return false
	 */
	virtual bool fractured() const ;

	/** \brief return a copy of the behaviour
	 *
	 * @return a new KelvinVoight
	 */
	virtual Form * getCopy() const ;

	/** \brief return true*/
	virtual bool changed() const ;

	virtual void scale( double d )
	{
		param *= d ;
		eta *= d ;
	}
} ;

struct IncrementalKelvinVoight : public LinearForm
{
	Matrix stiff ;
	Matrix eta ;
	std::vector<Variable> v ;
	Matrix N ;
	Vector phi ;
	double tau ;
	std::valarray<bool> up ;

	IncrementalKelvinVoight( const Matrix & rig, const Matrix & eta, double dt ) ;

	virtual ~IncrementalKelvinVoight() ;

	virtual void apply( const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm ) const ;

	virtual bool fractured() const
	{
		return false ;
	}

	virtual bool changed() const
	{
		return false ;
	}

	virtual Form * getCopy() const ;

	virtual void scale( double d )
	{
		param *= d ;
		eta *= d ;
		stiff *= d ;
	}

	virtual Vector getImposedStress( const Point &p ) const ;

	virtual std::vector<BoundaryCondition * > getBoundaryConditions( const ElementState &s, size_t id, const Function &p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv ) const ;

	virtual void step( double timestep, ElementState &s ) ;

	void resize( size_t num_points ) ;
};

} ;

#endif
