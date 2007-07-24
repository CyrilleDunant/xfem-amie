//
// C++ Interface: samplingcriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SAMPLINGCRITERION_H
#define SAMPLINGCRITERION_H

#include "delaunay.h"
#include "delaunay_3d.h"
namespace Mu
{

class DelaunayTriangle ;

/** \todo ajouter un SamplingCriterion InGeometry*/

/**
This class defines a sampling criterion for mesh refinement.

@author Cyrille Dunant
*/

class SamplingCriterion
{
protected:
	bool toAdd ;
public:
	SamplingCriterion() ;

	virtual ~SamplingCriterion() ;
	
	/*! returns true if triangle is valid for the given criterion*/
	virtual bool meetsCriterion(const DelaunayTriangle * t) ;
	virtual void reset()  = 0;
	virtual bool add() const ;
	virtual std::vector<Point> suggest(const DelaunayTriangle * t) const  = 0;
};

class MinimumAngle : public SamplingCriterion
{
protected:
	double angle ;
public:
	MinimumAngle(double a) ;
	
	virtual ~MinimumAngle() ;
	
	virtual bool meetsCriterion(const DelaunayTriangle * t) ;
	virtual std::vector<Point> suggest(const DelaunayTriangle * t) const ;
	virtual void reset()  ;
};

class MaximumLength : public SamplingCriterion
{
protected:
	double length ;
public:
	MaximumLength(double l) ;
	
	virtual ~MaximumLength() ;
	
	virtual bool meetsCriterion(const DelaunayTriangle * t) ;
	virtual void reset()  ;
};

class Counter : public SamplingCriterion
{
protected:
	size_t c ;
	size_t s ;
public:
	Counter(size_t nit) ;
	
	virtual ~Counter() ;
	
	virtual bool meetsCriterion(const DelaunayTriangle * t)  ;
	virtual void reset()  ;
};

};

#endif


