/*
    implementation of MAD for AMIE
    Copyright (C) 2010  Cyrille Dunant

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#ifndef MAD_H
#define MAD_H

#include "physics_base.h"
#include "fracturecriteria/fracturecriterion.h"
#include "damagemodels/damagemodel.h"
#include "../features/features.h"
#include "homogenization/elastic_homogenization.h"
#include "homogenization/converter.h"
#include "../utilities/optimizer.h"
#include "../features/crack.h"

namespace Mu
{
	struct MultipleAggregatingDiscontinuities : public LinearForm
	{
		void homogenize() ;
		
		std::vector<Variable> v ;
		FractureCriterion * crit ;
		FeatureTree * featureTree ;
		DelaunayTriangle * self2d ;
		DelaunayTetrahedron * self3d ;
		GeneralizedSelfConsistent scheme ;
		ElementDefinedBoundaryCondition * bc ;
		LeastSquaresApproximation * leastSquares ;
		BranchedCrack * equivalentCrack ;
		Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2D ;
		
		bool change ;
		/** \brief Constructor
		* 
		*/
		MultipleAggregatingDiscontinuities(FeatureTree * featureTree, DelaunayTriangle * self ,  FractureCriterion * crit) ;
		
		MultipleAggregatingDiscontinuities(FeatureTree * featureTree, DelaunayTetrahedron * self ,  FractureCriterion * crit) ;
		
		virtual ~MultipleAggregatingDiscontinuities() 
		{ 
			delete bc ; 
			delete featureTree ;
		};

		virtual XMLTree * toXML() { return new XMLTree("MultipleAggregatingDiscontinuities", param) ; } ;
		
		/** \brief Apply the law.
		* The apparent stifness tensor is computed from a self-consistent scheme applied on the submesh.
		*/
		virtual Matrix apply(const Function & p_i, const Function & p_j, const IntegrableEntity *e) const; 
		
		/** \brief Apply the law.
			*
			* The matrix is computed as: \f$ \nabla^T h_i K \nabla h_j \f$
			* @param p_i first basis polynomial.
			* @param p_j second basis polynomial.
			* @param gp Gauss Points used for the quadrature
			* @param Jinv Inverse Jacobian Matrices corresponding to the gauss points
			* @param ret Matrix to store the result
			* @param vm virtualMachine to use to compute the result
			*/
		virtual void apply(const Function & p_i, const Function & p_j, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Matrix & ret, VirtualMachine * vm) const ;
		
		/** \brief update behaviour
		* Apparent cracks are computed and enrich the element.
		* 
		* @param timestep elapsed time
		* @param currentState state of the element
		* 
		*/
		virtual void step(double timestep, ElementState & currentState) ;
		
		virtual bool changed() const ;
		virtual Matrix getTensor(const Point & p) const ;
		
		/** \brief Return false.*/
		virtual bool fractured() const ;
		
		/** \brief Return a copy of the behaviour*/
		virtual Form * getCopy() const ;
		
		/** \brief Return a 0-length Vector*/
		virtual void getForces(const ElementState & s, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv, Vector &v) const ;
	};
} ;
#endif // MAD_H
