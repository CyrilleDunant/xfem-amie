 //
// C++ Interface: multi-grid solver
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "solver.h"
#include "conjugategradient.h"
#include "gausseidell.h"
#include "inversediagonal.h"

namespace Amie 
{
/** \brief MultiGrid for symmetric systems*/
	template<class MESH_T, class ETYPE>
	struct MultiGrid : public LinearSolver
	{
		std::vector<Assembly *> A1 ;
		MESH_T * mesh0 ;
		std::vector<MESH_T *> mesh1 ;
		Vector r0 ;
		Vector r1 ;
		ConjugateGradient * cg0 ;
		GaussSeidel * gs0 ;
		LinearSolver * subsolver ;
		bool coarseSolved ;
		std::vector<ETYPE *> elements0 ;
		std::vector<ETYPE *> elements1 ;
		virtual ~MultiGrid() 
		{ 
			delete cg0 ;
			delete gs0 ;
			delete subsolver ;
		} ;
		MultiGrid(Assembly * a, const std::vector<Assembly *> & A1, MESH_T * mesh0, const std::vector<MESH_T *> & mesh1, Vector & f) : LinearSolver(a), A1(A1), mesh0(mesh0), mesh1(mesh1), r0(0.,a->getMatrix().stride()*a->getMaxDofID()), r1(A1.back()->getMatrix().row_size.size()*A1.back()->getMatrix().stride)  
		{ 
			coarseSolved = false ;
			cg0 = new ConjugateGradient(a) ;
			gs0 = new GaussSeidel(a) ;
			if(mesh1.size() == 1)
			{
				subsolver = new ConjugateGradient(A1.back()) ;
			}
			else if (mesh1.size() > 1)
			{
				std::vector<Assembly *> A1copy(A1.begin(), A1.end()-1) ;
				std::vector<MESH_T *> mesh1copy(mesh1.begin(), mesh1.end()-1) ;
				subsolver = new MultiGrid<MESH_T, ETYPE>(*(A1.back()), A1copy, mesh1.back(), mesh1copy, r1) ;
			}
			elements0 = mesh0->getElements() ;
			elements1 = mesh1.back()->getElements() ;

		};
		
		virtual bool solve(const Vector &x0, Preconditionner * precond = nullptr, const double eps = 5e-8, const int maxit = -1, bool verbose = true)
		{
			
			if(mesh1.empty())
			{
				std::cout << "No Coarse grids, falling back to CG" << std::endl ;
				bool solve = cg0->solve(x0, nullptr, eps, maxit, verbose) ;
				x = cg0->x ;
				return solve ;
				
			}

			double smoothingFactor = 0 ;
			int smoothingSteps = 10 ;
			
			int Maxit = maxit ;
			int nit = 0 ;
			if(maxit == -1)
				Maxit = x.size() ;

			
			if(x0.size() == assembly->getForces().size())
				x = x0 ;
			else
				x = 0. ;
			
			assign(r0, assembly->getMatrix()*x-assembly->getForces()) ;
			r0 = -r0 ;			
			if(std::abs(r0).max() < eps)
			{
				if(verbose)
				{
					std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< eps << std::endl  ;
				}
				return true ;
			}
			

			
			
			while(nit < Maxit )
			{

				nit++ ;
				if(!coarseSolved)
				{
					cg0->solve(x, nullptr, std::max(std::abs(r0).max()*smoothingFactor, eps), smoothingSteps, verbose) ;
					x = cg0->x ;
					assign(r0, assembly->getMatrix()*x-assembly->getForces()) ;
					r0 = -r0 ;
					if(verbose && nit%10 == 0)
						std::cout << A1.size() << "  "<< std::abs(r0).max() << "  " << eps  << std::endl ;
					if(std::abs(r0).max() < eps)
					{

						if(true)
							std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< std::abs(r0).max() << std::endl  ;

						return true ;
					}
					

					//restrict
					mesh1.back()->project(mesh0, r1, r0, false) ;
					
					//V iteration
					subsolver->assembly->getForces() = r1 ;
					subsolver->x = 0 ;
					subsolver->solve(subsolver->x, nullptr, eps, -1, false) ;
						
					//extend
					mesh0->project(mesh1.back(), r0, subsolver->x, false) ;
					x -= r0 ;

					if(std::abs(subsolver->x).max() < eps)
					{
						coarseSolved = true ;
						std::cout << "Grid " << A1.size()-1 << " converged." << std::endl ;
					}
					else
					{
						std::cout << A1.size()-1 << " not converged... " << std::flush ;
					}
					
					cg0->solve(x, nullptr, std::max(std::abs(r0).max()*smoothingFactor, eps), smoothingSteps, verbose) ;
					x = cg0->x ;
				}
				else
				{
					// finish computation on fine grid ;
					bool solved = cg0->solve(x, nullptr, eps, -1, verbose) ;
					x = cg0->x ;
// 					coarseSolved = false ;
					if(solved)
					{
						if(true)
							std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< eps << std::endl  ;
						
						return true ;
					}
					else
					{
						if(true)
							std::cout << "Grid " << A1.size() << " not solved in " << nit << " iterations, err = "<< std::abs((Vector)(assembly->getMatrix()*x-assembly->getForces())).max() << std::endl  ;
						
						return false ;
					}
				}
			}
			
			if(verbose)
				std::cout << "Grid " << A1.size() << " not solved in " << nit << " iterations, err = "<< std::abs(r0).max() << std::endl  ;
			
			return false ;
		} ;
	} ;

} 

#endif
