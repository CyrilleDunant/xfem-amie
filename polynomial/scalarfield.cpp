//
// C++ Implementation: scalarfield
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "scalarfield.h"

namespace Mu {

LineScalarField::LineScalarField(const Vector & vals, Point o, double s): values(vals), offset(o), scale(s)
{
}


LineScalarField::~LineScalarField()
{
}


const Vector & LineScalarField::getValues() const
{
	return values ;
}

const Point & LineScalarField::getOffset() const
{
	return offset ;
}

const double & LineScalarField::getScale() const
{
	return scale ;
}

double LineScalarField::operator() (double x) const
{
	double indexLow = floor((x-offset.x)*(values.size()-1)/scale) ;
	double indexHigh = ceil((x-offset.x)*(values.size()-1)/scale) ;
	if(indexLow < 0)
		return 0. ;
	if(indexHigh >= values.size())
		return 0. ;

	double xlow = indexLow/(values.size()-1)*scale + offset.x;
	double xhigh = indexHigh/(values.size()-1)*scale + offset.x;
	
	return ((x-xlow)*values[indexHigh] + (xhigh-x)*values[indexLow])/(xhigh-xlow) ;
}


SurfaceScalarField::SurfaceScalarField(const std::valarray<LineScalarField> & vals) : values(vals)
{
	
}

SurfaceScalarField::~SurfaceScalarField() { }

const std::valarray<LineScalarField> & SurfaceScalarField::getValues() const
{
	return values ;
}

const Point & SurfaceScalarField::getOffset() const
{
	return values[0].getOffset() ;
}

const double & SurfaceScalarField::getScale() const
{
	return values[0].getScale() ;
}


double SurfaceScalarField::operator() (double x, double y) const
{
	if(!values.size())
		return 0 ;
	
	double scale = values[0].getScale() ;
	Point offset = values[0].getOffset() ;
	double indexLow = floor((y-offset.y)*(values.size()-1)/scale) ;
	double indexHigh = ceil((y-offset.y)*(values.size()-1)/scale) ;
	
	if(indexLow < 0)
		return 0. ;
	if(indexHigh >= values.size())
		return 0. ;
	
	double xlow = values[indexLow](x) ;
	double xhigh = values[indexHigh](x) ;
	
	double ylow = indexLow/(values.size()-1)*scale + offset.y;
	double yhigh = indexHigh/(values.size()-1)*scale + offset.y;
	
	return ((y-ylow)*xhigh + (yhigh-y)*xlow)/(yhigh-ylow) ;
}

VolumeScalarField::VolumeScalarField(const std::valarray<SurfaceScalarField> & vals) : values(vals)
{
	
}

VolumeScalarField::~VolumeScalarField() { }

const std::valarray<SurfaceScalarField> & VolumeScalarField::getValues() const
{
	return values ;
}

const Point & VolumeScalarField::getOffset() const
{
	return values[0].getOffset() ;
}

const double & VolumeScalarField::getScale() const
{
	return values[0].getScale() ;
}


double VolumeScalarField::operator() (double x, double y, double z) const
{
	if(!values.size())
		return 0 ;
	
	double scale = values[0].getScale() ;
	Point offset = values[0].getOffset() ;
	double indexLow = floor((z-offset.z)*(values.size()-1)/scale) ;
	double indexHigh = ceil((z-offset.z)*(values.size()-1)/scale) ;
	
	if(indexLow < 0)
		return 0. ;
	if(indexHigh >= values.size())
		return 0. ;
	
	double xylow = values[indexLow](x,y) ;
	double xyhigh = values[indexHigh](x,y) ;
	
	double zlow = indexLow/(values.size()-1)*scale + offset.z;
	double zhigh = indexHigh/(values.size()-1)*scale + offset.z;
	
	return ((z-zlow)*xyhigh + (zhigh-z)*xylow)/(zhigh-zlow) ;
}

}
