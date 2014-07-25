//
// C++ Interface: multi-grid step for preconditionning purposes
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MULTIGRIDSTEP_H
#define MULTIGRIDSTEP_H

#include "preconditionners.h"
#include "gaussseidellstep.h"
#include "inversediagonal.h"
#include "multigrid.h"


namespace Mu {

/** \brief Preconditionner, perform a GS step
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
template<class MESH_T, class ETYPE>
struct MultiGridStep : public Preconditionner
{
	ConjugateGradient * subsolver ;
	std::vector< ETYPE *> elements0 ;
	std::vector< ETYPE *> elements1 ;
	MESH_T * mesh0; 
	MESH_T * mesh1 ;
	Vector v1 ;
	Vector b ;
	double averageRadius ;
	std::vector< Point *> points ;
	std::vector<std::vector<Point *> > smoothingPoints ;
	InverseDiagonal ivd ;
	int count ;
	virtual ~MultiGridStep()  
	{ 
		delete subsolver ;
	};
	
	MultiGridStep(MESH_T * mesh0, MESH_T * mesh1, const CoordinateIndexedSparseMatrix * A0,  const CoordinateIndexedSparseMatrix * A1, Vector &b) : elements0(mesh0->getElements()), elements1(mesh1->getElements()), mesh0(mesh0), mesh1(mesh1), v1(b), b(b), averageRadius(0), ivd(*A0)
	{ 
		count = 0 ;
		subsolver = new ConjugateGradient(*A1, b) ;
// 		subsolver->solve(b, nullptr, 1e-9) ;
// 		b = subsolver->getX() ;
		
// 		for(size_t i = 0 ; i < elements0.size() ; i++)
// 		{
// 			averageRadius += elements0[i]->getRadius()/elements0.size() ;
// 			for(size_t j = 0 ; j < elements0[i]->getBoundingPoints().size() ; j++)
// 			{
// 				points.push_back(&elements0[i]->getBoundingPoint(j)) ;
// 			}
// 			std::sort(points.begin(), points.end()) ;
// 			points.erase(std::unique(points.begin(), points.end()), points.end()) ;
// 		}
		
// 		for(size_t i = 0 ; i < points.size() ; i++)
// 		{
// 			Circle c(averageRadius*0.1, *points[i]) ;
// 			std::vector< ETYPE *> elems = mesh0->getConflictingElements(&c) ;
// 			std::vector<Point*> pts ;
// 			for(size_t j = 0 ; j < elems.size() ; j++)
// 			{
// 				for(size_t k = 0 ; k < elems[j]->getBoundingPoints().size() ; k++)
// 				{
// 					pts.push_back(&elems[j]->getBoundingPoint(k)) ;
// 				}
// 			}
// 			std::sort(pts.begin(), pts.end()) ;
// 			pts.erase(std::unique(pts.begin(), pts.end()), pts.end()) ;
// 			smoothingPoints.push_back(pts) ;
// 		}
		
	};
	
	virtual void precondition(const Vector &v,Vector & t)  
	{
		Vector vorig(v) ;
		
		count++ ;
		if(true)
		{
			//restrict
// 			for(size_t i = 0 ; i < elements0.size() ; i++)
// 				elements0[i]->step(1., &v) ;
			mesh1->project(mesh0, v1, v, false) ;
			
			//V iteration
			subsolver->b = -v1;
			subsolver->x = 0 ;
			subsolver->solve(subsolver->x, nullptr, 1e-9, -1, false) ;
			subsolver->nit = 0 ;
			
			//extend
// 			for(size_t i = 0 ; i < elements1.size() ; i++)
// 				elements1[i]->step(1., &subsolver->getX()) ;
			mesh0->project(mesh1, t, subsolver->x, false) ;

			t = vorig-t ;
		}
		else
		{
			ivd.precondition(v, t) ;
		}
			
// 		Vector torig(t) ;
// 		for(size_t i = 0 ; i < points.size() ; i++)
// 		{
// 			for(size_t j = 0 ; j < smoothingPoints[i].size() ; j++)
// 			{
// 				t[points[i]->getId()*2] += torig[smoothingPoints[i][j]->getId()*2]/smoothingPoints[i].size() ;
// 				t[points[i]->getId()*2+1] += torig[smoothingPoints[i][j]->getId()*2+1]/smoothingPoints[i].size() ;
// 			}
// 		}
// 		t = t*.25+v*0.75 ;

// 		ivd.precondition(t,t) ;
// 		t = v ;
	}

};

}

#endif
