//
// C++ Interface: scalarfield
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2008-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUSCALARFIELD_H
#define MUSCALARFIELD_H
#include "../utilities/matrixops.h"
#include "../geometry/geometry_base.h"

namespace Amie {

/**
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	This scalar field returns the value interpolated at
	x of a field. The values in the field are given as 
	regularly sampled on a line.
*/
class LineScalarField{
protected:
	Vector values ;
	Point offset ;
	double scale ;
public:
	LineScalarField(const Vector & vals, Point offset, double scale);

	~LineScalarField();
	
	double operator() (double x) const ;
	
	const Vector & getValues() const ;
	const Point & getOffset() const ;
	const double & getScale() const ;
	

};

/**
@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	This scalar field returns the value interpolated at
	x, y of a field. The values in the field are given 
	as regularly sampled on a square surface. If the
	sampling point number is not the same in the x and
	y direction, the results will be distorted.
*/
class SurfaceScalarField
{
protected:
	std::valarray<LineScalarField> values;
public:
	SurfaceScalarField(const std::valarray<LineScalarField> & vals);
	
	~SurfaceScalarField();
	
	double operator() (double x, double y) const ;
	
	const std::valarray<LineScalarField> & getValues() const ;
	const Point & getOffset() const ;
	const double & getScale() const ;
	
};

/**
@author Cyrille Dunant <cyrille.dunant@epfl.ch>
	This scalar field returns the value interpolated at
	x, y, z of a field. The values in the field are given 
	as regularly sampled on a cubic surface. If the
	sampling point number is not the same in the x and
	y direction, the results will be distorted.
*/
class VolumeScalarField
{
protected:
	std::valarray<SurfaceScalarField> values;
public:
	VolumeScalarField(const std::valarray<SurfaceScalarField> & vals);
	
	~VolumeScalarField();
	
	double operator() (double x, double y, double z) const ;
	
	const std::valarray<SurfaceScalarField> & getValues() const ;
	const Point & getOffset() const ;
	const double & getScale() const ;
	
};

}

#endif
