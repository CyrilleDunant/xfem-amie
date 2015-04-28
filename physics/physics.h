// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __PHYSICS_H_
#define __PHYSICS_H_

#include "physics_base.h"
#include "void_form.h"
#include "stiffness.h"
#include "diffusion.h"
#include "weibull_distributed_stiffness.h"
#include "stiffness_and_fracture.h"
#include "stiffness_with_imposed_deformation.h"

namespace Amie
{



// class Diffusion2D :public LinearForm
// {
// public:
// 
// 	double c ;
// 	
// 	double alpha ;
// 	double tau   ;
// 	
// 	Diffusion2D(double capacity, double conductivityX, double conductivityY) ;
// 	
// 	virtual ~Diffusion2D() ;
// 	
// 	
// 	virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const;
// 	
// 	virtual Matrix apply(const Function & p_i, const Function & p_j, const std::valarray< std::pair<Point,double> > &gp, const std::valarray<Matrix> &Jinv, Matrix & ret) const;
// 	
// 	virtual void step(double timestep, ElementState * currentState) ;
// 	
// 	/** Check for fracture state
// 	 *
// 	 * @return true if the element is fractured
// 	 */
// 	virtual bool fractured() const;
// 	
// 	
// 	/** get Copy of the behaviour
// 	 *
// 	 * @return pointer to the copy. Caller is responsible for cleaning memory
// 	 */
// 	virtual Form * getCopy() const ;
// 
// 	virtual bool hasInducedForces() const ;
// 	
// 	virtual void getForces(const ElementState & s, const Function & p_i, const std::valarray< std::pair<Point, double> > &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
// 	
// } ;

class TwoDCohesiveForces : public NonLinearForm
{
public:

	std::vector<Point> normals ;
	const IntegrableEntity * source ;
	const IntegrableEntity * target ;
	double startArea ;
	
	bool active ;
	
	TwoDCohesiveForces(const IntegrableEntity *s, const IntegrableEntity *t, const SegmentedLine * sl) ;

	virtual Form * getCopy() const ;
	
	virtual ~TwoDCohesiveForces() ;
		
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
	
	virtual bool hasInducedForces() const ;
	
	virtual bool hasInducedMatrix() const ;
	
	virtual void getForces( ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v)  ;

	virtual void step(double timestep, ElementState & currentState) ;

	virtual bool isActive() const ;
} ;

/** ViscoElasticity law.
 */
struct ViscoElasticity: public LinearForm
{
	Vector tau_g ;
	Vector tau_k ;
	Vector g ;
	Vector k ;
	Matrix a_g ;
	Vector a_k ;
	Vector average_delta_sigma ;

	ViscoElasticity( double _tau_k, double _tau_g, Vector g, Vector k);
	
	virtual ~ViscoElasticity() ;
		
	virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const;
	/** \todo remove usage of previousState. complement state instead*/
	virtual void step(double timestep, ElementState & currentState);
	virtual void getForces(ElementState& s, const Function& p_i, const Amie::GaussPointArray& gp, const std::valarray< Amie::Matrix >& Jinv, Vector& f) ;
	
	virtual bool hasInducedForces();
	
	virtual Form * getCopy() const ;
} ;

} ;

#endif // __ PHYSICS_H_



