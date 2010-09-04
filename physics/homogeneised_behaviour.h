//
// C++ Interface: stiffness
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __HOMOGENEISED_H_
#define __HOMOGENEISED_H_


#include "physics_base.h"
#include "homogenization/elastic_homogenization.h"
#include "../mesher/mesh.h"
#include "../mesher/delaunay.h"
#include "../mesher/delaunay_3d.h"

namespace Mu
{
	class Feature ;
	class FeatureTree ;

	/** \brief A linear Elastic Law
	* The field param is the Cauchy-Green Strain Tensor
	*/
	struct HomogeneisedBehaviour : public LinearForm
	{
		Form * equivalent ;
		Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2d ;
		DelaunayTriangle * self2d ;
		Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * mesh3d ;
		DelaunayTetrahedron * self3d ;
		std::vector<Variable> v ;
		std::vector<DelaunayTriangle *> source ;
		std::vector<DelaunayTetrahedron *> source3d ;
                std::vector<Feature *> ft ;

		/** \brief Constructor
		* 
		* @param mesh2d The 2D mesh for 2D homogeneisation
		* @param self The element in which the homogeneisation will occur
		*/
		HomogeneisedBehaviour(Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh2d, DelaunayTriangle * self) ;
		/** \brief Constructor
		* 
		* @param mesh3d The 3D mesh for 3D homogeneisation
		* @param self The element in which the homogeneisation will occur
		*/
		HomogeneisedBehaviour(Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * mesh3d, DelaunayTetrahedron * self) ;

		HomogeneisedBehaviour(std::vector<Feature *> feats, DelaunayTriangle * self) ;


		virtual ~HomogeneisedBehaviour() ;
		
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
		
		/** \brief Return false.*/
		virtual bool fractured() const ;
		
		/** \brief Return a copy of the behaviour*/
		virtual Form * getCopy() const ;

		/** \brief Homogenizes the elastic behaviour*/
		void homogenize() ;

		void refresh() ;

//		void setFeatureTree(FeatureTree ft) {subTree = ft ; } ;

		virtual void step(double timestep, ElementState & currentState) ;
		virtual void stepBack() ;

		Material toMaterial() {return equivalent->toMaterial() ; } ;
		Material homogenize(Material mat) ;

                virtual Vector getImposedStress(const Point & p) const ;

                std::vector<BoundaryCondition * > getBoundaryConditions(const ElementState & s,  size_t id, const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const ;

		Form * getEquivalentBehaviour(Material mat) ;	

	} ;


} ;

#endif //__HOMOGENEISED_H_
