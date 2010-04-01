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
		
		virtual ~MultiGrid() { } ;
		MultiGrid(const CoordinateIndexedSparseMatrix & A0, const std::vector<const CoordinateIndexedSparseMatrix *> & A1, MESH_T * mesh0, const std::vector<MESH_T *> & mesh1, Vector & f) : LinearSolver(A0, f), A1(A1), mesh0(mesh0), mesh1(mesh1)  { };
		
		virtual bool solve(const Vector &x0, const Preconditionner * precond = NULL, const double eps = 1e-8, const int maxit = -1, bool verbose = true)
		{
			LinearSolver * subsolver = NULL;
			int Maxit = maxit ;
			Preconditionner * p1 = NULL ;
			InverseDiagonal * p0 = new InverseDiagonal(A) ;
			ConjugateGradient cg0(A, b) ;
// 			GaussSeidel cg0(A, b) ;
			int nit = 0 ;
			if(x0.size() == b.size())
			{
				x.resize(x0.size()) ;
				x = x0 ;
			}
			else
			{
				if(x0.size())
				{
					std::cout << "ouch" << std::endl ;
					exit(0) ;
				}
				x.resize(b.size()) ;
				x = 0. ;
			}
			
			
			if(maxit == -1)
				Maxit = x.size() ;
			bool solve = cg0.solve(x, p0, std::max(std::abs(b).max()*.5, eps), -1, true) ;
// 			bool solve = cg0.solve(x, p0, eps, 10, false) ;
			x = cg0.x ;
			Vector r0 = A*x-b ;
			
			if(std::abs(r0).max() < eps)
			{
				delete p0 ;
// 				if(verbose)
// 				{
					std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< eps << std::endl  ;
// 				}
				return true ;
			}
			
			
			//mesh0 refresh elements with r0 ;
			std::vector<ETYPE *> elements0 = mesh0->getElements() ;
			std::vector<ETYPE *> elements1 = mesh1.back()->getElements() ;
			
			for(size_t i = 0 ; i < elements0.size() ; i++)
				elements0[i]->step(1., &r0) ;
			Vector r1(A1.back()->row_size.size()*A1.back()->stride) ;
			mesh1.back()->project(mesh0, r1, r0) ;
			
			if(mesh1.size() == 1)
			{
				p1 = new InverseDiagonal(*A1.back()) ;
				subsolver = new ConjugateGradient(*A1.back(), r1) ;
			}
			else
			{
				std::vector<const CoordinateIndexedSparseMatrix *> A1copy(A1.begin(), A1.end()-1) ;
				std::vector<MESH_T *> mesh1copy(mesh1.begin(), mesh1.end()-1) ;
				subsolver = new MultiGrid<MESH_T, ETYPE>(*(A1.back()), A1copy, mesh1.back(), mesh1copy, r1) ;
			}

			subsolver->solve(r1, p1, eps, -1, false) ;

			//mesh1 refresh elements with r1
			for(size_t i = 0 ; i < elements1.size() ; i++)
				elements1[i]->step(1., &subsolver->x) ;
			mesh0->project(mesh1.back(), r0, subsolver->x) ;
			x -= r0 ;
			
			bool coarseSolved = false ;
			while(nit < Maxit )
			{
				if(!coarseSolved)
				{
					solve = cg0.solve(x, p0, std::max(std::abs(r0).max()*.5, eps), -1, true) ;
	//  				solve = cg0.solve(x, p0, eps, 10, false) ;
					x = cg0.x ;
					r0 = A*x-b ;
					if(std::abs(r0).max() < eps)
					{
						delete subsolver ;
						delete p0 ;
						delete p1 ;
	// 					if(verbose)
	// 					{
							std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< std::abs(r0).max() << std::endl  ;
	// 					}
						return true ;
					}
					
					nit++ ;
					if(verbose)
						std::cout << A1.size() << " " << nit << ", err = " << std::abs(r0).max() << std::endl ;
					//mesh0 refresh elements with r ;

					for(size_t i = 0 ; i < elements0.size() ; i++)
						elements0[i]->step(1., &r0) ;
					mesh1.back()->project(mesh0, r1, r0) ;
					
					subsolver->b = r1 ;
					subsolver->solve(subsolver->x, p1, eps, -1, false) ;
						
					//mesh1 refresh elements with r1
					for(size_t i = 0 ; i < elements1.size() ; i++)
						elements1[i]->step(1., &subsolver->x) ;
					
					if(std::abs(subsolver->x).max() < eps)
						coarseSolved = true ;
					
					mesh0->project(mesh1.back(), r0, subsolver->x) ;
					x -= r0 ;
				}
				else
				{
					bool solved = cg0.solve(x, p0, eps, -1, true) ;
					x = cg0.x ;
					if(solved)
					{
						std::cout << "Grid " << A1.size() << " solved in " << nit << " iterations, err = "<< eps << std::endl  ;
						return true ;
					}
					else
					{
						std::cout << "Grid " << A1.size() << " not solved in " << nit << " iterations, err = "<< std::abs((Vector)(A*x-b)).max() << std::endl  ;
						return false ;
					}
				}
			}
// 			if(verbose)
				std::cout << "Grid " << A1.size() << " not solved in " << nit << " iterations, err = "<< std::abs(r0).max() << std::endl  ;
			delete p0 ;
			delete p1 ;
			delete subsolver ;
			return false ;
		} ;
	} ;

} ;

#endif
