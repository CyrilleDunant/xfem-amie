 //
// C++ Interface: multi-grid solver
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
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

namespace Mu 
{
/** \brief MultiGrid for symmetric systems*/
	template<class MESH_T, class ETYPE>
	struct MultiGrid : public LinearSolver
	{
		std::vector<const CoordinateIndexedSparseMatrix *> A1 ;
		MESH_T * mesh0 ;
		std::vector<MESH_T *> mesh1 ;
		Vector r0 ;
		Vector r1 ;
		ConjugateGradient * cg0 ;
		LinearSolver * subsolver ;
		virtual ~MultiGrid() 
		{ 
			delete cg0 ;
			delete subsolver ;
		} ;
		MultiGrid(const CoordinateIndexedSparseMatrix & A0, const std::vector<const CoordinateIndexedSparseMatrix *> & A1, MESH_T * mesh0, const std::vector<MESH_T *> & mesh1, Vector & f) : LinearSolver(A0, f), A1(A1), mesh0(mesh0), mesh1(mesh1), r0(b.size()), r1(A1.back()->row_size.size()*A1.back()->stride)  
		{ 
			cg0 = new ConjugateGradient(A0, f) ;
			
			if(mesh1.size() == 1)
			{
				subsolver = new ConjugateGradient(*A1.back(), r1) ;
			}
			else if (mesh1.size() > 1)
			{
				std::vector<const CoordinateIndexedSparseMatrix *> A1copy(A1.begin(), A1.end()-1) ;
				std::vector<MESH_T *> mesh1copy(mesh1.begin(), mesh1.end()-1) ;
				subsolver = new MultiGrid<MESH_T, ETYPE>(*(A1.back()), A1copy, mesh1.back(), mesh1copy, r1) ;
			}
		};
		
		virtual bool solve(const Vector &x0, const Preconditionner * precond = NULL, const double eps = 1e-8, const int maxit = -1, bool verbose = true)
		{
			if(mesh1.empty())
			{
				std::cout << "No Coarse grids, falling back to CG" << std::endl ;
				bool solve = cg0->solve(x0, NULL, eps, maxit, verbose) ;
				x = cg0->x ;
				return solve ;
				
			}
// 			GaussSeidel cg0(A, b) ;
			double smoothingFactor = 1.1 ;
			int smoothingSteps = -1 ;
			
			int Maxit = maxit ;
			int nit = 0 ;
			if(maxit == -1)
				Maxit = x.size() ;
			
// 			
			std::vector<ETYPE *> elements0 = mesh0->getElements() ;
			std::vector<ETYPE *> elements1 = mesh1.back()->getElements() ;

			
			if(x0.size() == b.size())
				x = x0 ;
			else
				x = 0. ;
			
			assign(r0, A*x-b) ;
			
			if(std::abs(r0).max() < eps)
			{
				if(verbose)
				{
					std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< eps << std::endl  ;
				}
				return true ;
			}
			

			
			bool coarseSolved = false ;
			while(nit < Maxit )
			{
				int toto ;

				nit++ ;
				if(!coarseSolved)
				{
					if(nit%10 == 0)
						std::cout << A1.size() << " " << nit << ", err = " << std::abs(r0).max() << std::endl ;
					
					cg0->solve(x, NULL, std::max(std::abs(r0).max()*smoothingFactor, eps), smoothingSteps, verbose) ;
					x = cg0->x ;
					assign(r0, A*x-b) ;
					
					if(std::abs(r0).max() < eps)
					{

						if(verbose)
							std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< std::abs(r0).max() << std::endl  ;

						return true ;
					}
					
// 					if(verbose)
						

					//restrict
					mesh1.back()->project(mesh0, r1, r0, true) ;
					
					//V iteration
					subsolver->b = r1 ;
					subsolver->solve(subsolver->x, NULL, eps, -1, false) ;
						
					//extend
					mesh0->project(mesh1.back(), r0, subsolver->x, true) ;
					x -= r0 ;					
					
					if(std::abs(subsolver->x).max() < eps)
						coarseSolved = true ;
				}
				else
				{
					// finish computation on fine grid ;
					bool solved = cg0->solve(x, NULL, eps, Maxit-nit, verbose) ;
					x = cg0->x ;
					
					if(solved)
					{
						if(verbose)
							std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< eps << std::endl  ;
						
						return true ;
					}
					else
					{
						if(verbose)
							std::cout << "Grid " << A1.size() << " not solved in " << nit << " iterations, err = "<< std::abs((Vector)(A*x-b)).max() << std::endl  ;
						
						return false ;
					}
				}
			}
			
			if(verbose)
				std::cout << "Grid " << A1.size() << " not solved in " << nit << " iterations, err = "<< std::abs(r0).max() << std::endl  ;
			
			return false ;
		} ;
	} ;

} ;

#endif
