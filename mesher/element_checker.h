/*
    mesh abstract implementation for AMIE
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


#ifndef ELEMENT_CHECKER_H
#define ELEMENT_CHECKER_H

#include "../elements/integrable_entity.h"

namespace Amie
{

struct ElementChecker
{
	ElementChecker() { } ;
	virtual bool checkElement(IntegrableEntity * e) const { return (bool) e ; }
} ;

struct StiffnessHigherThanElementChecker : public ElementChecker
{
	double threshold ;
	size_t row ;
	size_t col ;
	StiffnessHigherThanElementChecker(double t, size_t r = 0, size_t c = 0) : ElementChecker(), threshold(t), row(r), col(c) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		if(e->getBehaviour() )
		{ 
			Point p ;
			return e->getBehaviour()->getTensor(p)[row][col] > threshold ;
		}
		return false ;
	}	
} ;

struct StiffnessLowerThanElementChecker : public StiffnessHigherThanElementChecker
{
	StiffnessLowerThanElementChecker(double t, size_t r = 0, size_t c = 0) : StiffnessHigherThanElementChecker(t,r,c) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		if(e->getBehaviour() )
		{ 
			Point p ;
			return e->getBehaviour()->getTensor(p)[row][col] < threshold ;
		}
		return false ;
	}	
} ;

struct ViscosityHigherThanElementChecker : public StiffnessHigherThanElementChecker
{
	ViscosityHigherThanElementChecker(double t, size_t r = 0, size_t c = 0) : StiffnessHigherThanElementChecker(t,r,c) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		if(e->getBehaviour() )
		{ 
			Point p ;
			return e->getBehaviour()->getViscousTensor(p)[row][col] > threshold ;
		}
		return false ;
	}	
} ;

struct ViscosityLowerThanElementChecker : public StiffnessHigherThanElementChecker
{
	ViscosityLowerThanElementChecker(double t, size_t r = 0, size_t c = 0) : StiffnessHigherThanElementChecker(t,r,c) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		if(e->getBehaviour() )
		{ 
			Point p ;
			return e->getBehaviour()->getViscousTensor(p)[row][col] < threshold ;
		}
		return false ;
	}	
} ;

struct FieldHigherThanElementChecker : public ElementChecker
{
	FieldType field ;
	double threshold ;
	double instant ;
	FieldHigherThanElementChecker(FieldType f, double t, double i = 1.) : ElementChecker(), field(f), threshold(t), instant(i) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		if(e->getBehaviour())
		{ 
			size_t blocks = e->getBehaviour()->getNumberOfDegreesOfFreedom() / e->spaceDimensions() ;
			Vector values ( 0., fieldTypeElementarySize ( field, e->spaceDimensions(), blocks ) ) ;
			e->getState().getAverageField( field, values, nullptr, -1., 1.) ;
			return values.max() > threshold ;
		}
		return false ;
	}	
} ;

struct FieldLowerThanElementChecker : public ElementChecker
{
	FieldType field ;
	double threshold ;
	double instant ;
	FieldLowerThanElementChecker(FieldType f, double t, double i = 1.) : ElementChecker(), field(f), threshold(t), instant(i) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		if(e->getBehaviour())
		{ 
			size_t blocks = e->getBehaviour()->getNumberOfDegreesOfFreedom() / e->spaceDimensions() ;
			Vector values ( 0., fieldTypeElementarySize ( field, e->spaceDimensions(), blocks ) ) ;
			e->getState().getAverageField( field, values, nullptr, -1., 1.) ;
			return values.min() < threshold ;
		}
		return false ;
	}	
} ;

/*struct InFeatureElementChecker : public ElementChecker
{
	std::vector<Feature *> inclusions ;

	InFeatureElementChecker(std::vector<Feature *> inc) : ElementChecker(), inclusions(inc) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		for(size_t i = 0 ; i < inclusions.size() ; i++)
		{
			if(inclusions[i]->in(e->getCenter()))
				return true ;
		}
		return false ;
	}	
	

} ;

struct NotInFeatureElementChecker : public InFeatureElementChecker
{
	NotInFeatureElementChecker(std::vector<Feature *> inc) : InFeatureElementChecker(inc) { } ;
	virtual bool checkElement(IntegrableEntity * e) const 
	{ 
		for(size_t i = 0 ; i < inclusions.size() ; i++)
		{
			if(inclusions[i]->in(e->getCenter()))
				return false ;
		}
		return true
 ;
	}	
	

} ;*/

}


#endif // ELEMENT_CHECKER_H
