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
#include "inversediagonal.h"

namespace Mu 
{
/** \brief preconditionned Conjugate Gradient for symmetric systems*/
	template<class MESH_T, class ETYPE>
	struct MultiGrid : public LinearSolver
	{
		const CoordinateIndexedSparseMatrix & A1 ;
		MESH_T * mesh0 ;
		MESH_T * mesh1 ;
		
		virtual ~MultiGrid() { } ;
		MultiGrid(const CoordinateIndexedSparseMatrix & A0, const CoordinateIndexedSparseMatrix & A1, MESH_T * mesh0, MESH_T * mesh1, Vector & f) : LinearSolver(A0, f), A1(A1), mesh0(mesh0), mesh1(mesh1)  { };
		
		virtual bool solve(const Vector &x0, const Preconditionner * precond = NULL, const double eps = 1e-9, const int maxit = -1, bool verbose = false)
		{
			LinearSolver * subsolver = NULL;
			int Maxit = maxit ;

			InverseDiagonal * p0 =new InverseDiagonal(A) ;
			ConjugateGradient cg0(A, b) ;
			if(x0.size() == b.size())
			{
				x.resize(x0.size()) ;
				x = x0 ;
			}
			else
			{
				x.resize(b.size()) ;
				x = 0. ;
			}
			
			if(maxit == -1)
				Maxit = x.size() ;
			cg0.solve(x, p0, 1e-10, 2, false) ;
			x = cg0.x ;
			Vector r0 = A*x-b ;
			int nit = 0 ;

			//mesh0 refresh elements with r0 ;
			std::vector<ETYPE *> elements0 = mesh0->getElements() ;
			for(size_t i = 0 ; i < elements0.size() ; i++)
				elements0[i]->step(1., &r0) ;
			Vector r1(A1.row_size.size()*A1.stride) ;
			mesh1->project(mesh0, r1) ;

			
			//if last grid
			subsolver = new ConjugateGradient(A1, r1) ;

			InverseDiagonal * p1 =new InverseDiagonal(A1) ;
			subsolver->solve(r1, p1, 1e-10, -1, true) ;

			//mesh1 refresh elements with r1
			std::vector<ETYPE *> elements1 = mesh1->getElements() ;
			for(size_t i = 0 ; i < elements1.size() ; i++)
				elements1[i]->step(1., &subsolver->x) ;
			mesh0->project(mesh1, r0) ;
			x -= r0 ;

			//iterate once
			bool coarseConverged = false ;
			while(nit < Maxit)
			{
				
				bool solve = cg0.solve(x, p0, 1e-10, 2, false) ;
				x = cg0.x ;
				if(solve)
				{
					delete subsolver ;
					delete p0 ;
					delete p1 ;
					return true ;
				}
				if( !coarseConverged)
				{
					r0 = A*x-b ;
					std::cout << nit++<< ", err = " << std::abs(r0).max()  << std::endl ;
					//mesh0 refresh elements with r ;

					for(size_t i = 0 ; i < elements0.size() ; i++)
						elements0[i]->step(1., &r0) ;
					mesh1->project(mesh0, r1) ;
					subsolver->b = r1 ;
					subsolver->solve(subsolver->x, p1, 1e-10, -1, true) ;
					if(std::abs(subsolver->x).max() < 1e-8)
						coarseConverged = true ;
					//mesh1 refresh elements with r1

					for(size_t i = 0 ; i < elements1.size() ; i++)
						elements1[i]->step(1., &subsolver->x) ;
					mesh0->project(mesh1, r0) ;
					std::cout << "err = " << std::abs(r0).max()  << std::endl ;
					x -= r0 ;
				}
				else
				{
					std::cout << nit++ << std::endl;
					bool solve = cg0.solve(x, p0, 1e-10, -1, true) ;
					x = cg0.x ;
					if(solve)
					{
						delete subsolver ;
						delete p0 ;
						delete p1 ;
						return true ;
					}
				}
			}
			delete subsolver ;
			return false ;
		} ;
	} ;

} ;

#endif
