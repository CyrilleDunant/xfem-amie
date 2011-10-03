
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

#ifndef MICROSTRUCTURE_GENERATOR_H
#define MICROSTRUCTURE_GENERATOR_H

#include "feature_base.h"
#include "inclusion.h"
#include "../geometry/geometry_base.h"
#include "../utilities/random.h"

namespace Mu
{
	struct MicrostructureGenerator
	{
		
		MicrostructureGenerator() { }; 
		virtual std::vector<Feature *> getFeatures(Geometry * sample) = 0;
		virtual double score() = 0;
		virtual void print() const = 0 ;
	} ;
	
	
	struct AggregateDistribution2DGenerator : public MicrostructureGenerator
	{
		
		double area ;
		double dmax ;
		double fill ;
		double minMaxRatio ;
		double itzSize ;
		double inclusionNumber ;
		double massOfAggregates ;
		Geometry * sample ;
		
		AggregateDistribution2DGenerator(double area, double dmax, double itzSize, double fill, double minMaxRatio) ;
		
		virtual std::vector<Feature *> getFeatures(Geometry * sample) ;
		virtual double score() ;
		virtual void print() const ;
	} ;
	
	/** \brief Utility class to convert circles to any distribution*/
	struct InclusionConverter
	{
		GeometryType geom ;
		RandomDistribution * area ;
		RandomDistribution * aspectRatio ;
		RandomDistribution * orientation ;
		
		InclusionConverter(GeometryType type, RandomDistribution * a = new ConstantDistribution(1.), RandomDistribution * ar = new ConstantDistribution(1.), RandomDistribution * o = new ConstantDistribution(0.)) ;
		
		void setArea(RandomDistribution * a) ;
		void setAspectRatio(RandomDistribution * ar) ;
		void setOrientation(RandomDistribution * o) ;
		
		void setArea(double a) ;
		void setAspectRatio(double ar) ;
		void setOrientation(double o) ;
		
		Feature * convert(Inclusion * inc) const ;
		std::vector<Feature *> convert(std::vector<Inclusion *> inc) const ;
		
		
	} ;
	
	
} ;





#endif // MICROSTRUCTURE_GENERATOR_H