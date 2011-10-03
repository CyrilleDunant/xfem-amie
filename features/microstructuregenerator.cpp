// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "microstructuregenerator.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/optimizer.h"
#include <iostream>

namespace Mu
{
	
	AggregateDistribution2DGenerator::AggregateDistribution2DGenerator(double area, double dmax, double itzSize, double fill, double minMaxRatio) : area(area), dmax(dmax), fill(fill), minMaxRatio(minMaxRatio), itzSize(itzSize)
	{
		
	}
	
	double AggregateDistribution2DGenerator::score()
	{
//		srand(0);
		std::vector<Inclusion *> inclusions ;
		
		if(dmax == 0.004)
			inclusions = GranuloBolome(exp(massOfAggregates), 1, BOLOME_D)(dmax*.5, .0001, inclusionNumber, itzSize, false);
		else
			inclusions = GranuloBolome(exp(massOfAggregates), 1, BOLOME_A)(dmax*.5, .0001, inclusionNumber, itzSize, false);

		std::vector<Feature *> feats ;
		for(size_t i = 0; i < inclusions.size() ; i++)
			feats.push_back(inclusions[i]) ;

		int nAgg = 1 ;
		feats = placement(sample, feats, &nAgg, 6400, false);
		
		double volume = 0 ;
		for(size_t i = 0 ; i < feats.size() ; i++)
		{
			inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
			volume += feats[i]->area() ;
		}
		if(!feats.empty())
		{
			std::cout << "mass = " << exp(massOfAggregates) << " n = " << feats.size() << ", largest r = " << feats.front()->getRadius() 
			<< ", smallest r =" << feats.back()->getRadius() << ", ratio = " << feats.front()->getRadius()/feats.back()->getRadius()
			<< ", filling = " << volume/sample->area()*100.<< "%"<< std::endl ; 
		

			double currentMinToMax = feats.front()->getRadius()/feats.back()->getRadius() ;
			double currentFill = volume/sample->area() ;
			
			for(int i = 0  ; i < inclusions.size() ; i++)
				delete inclusions[i] ;
			return 1./(.01*currentMinToMax/minMaxRatio + currentFill/fill) ;
		}
		else
		{
			for(int i = 0  ; i < inclusions.size() ; i++)
				delete inclusions[i] ;
			return 1000 ;
		}
	}

	void AggregateDistribution2DGenerator::print() const
	{
		std::vector<Inclusion *> inclusions ;
		
		if(dmax == 0.004)
			inclusions = GranuloBolome(massOfAggregates, 1, BOLOME_D)(dmax*.5, .00001, inclusionNumber, itzSize, false);
		else
			inclusions = GranuloBolome(massOfAggregates, 1, BOLOME_A)(dmax*.5, .00001, inclusionNumber, itzSize, false);

		std::vector<Feature *> feats ;
		for(size_t i = 0; i < inclusions.size() ; i++)
			feats.push_back(inclusions[i]) ;

		int nAgg = 1 ;
		feats = placement(sample, feats, &nAgg, 6400, false);
		
		double volume = 0 ;
		for(size_t i = 0 ; i < feats.size() ; i++)
			volume += feats[i]->area() ;
		if(!feats.empty())
			std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() 
			<< ", smallest r =" << feats.back()->getRadius() << ", ratio = " << feats.front()->getRadius()/feats.back()->getRadius()
			<< ", filling = " << volume/sample->area()*100.<< "%"<< std::endl ; 
		
	}
	
	std::vector<Feature *> AggregateDistribution2DGenerator::getFeatures(Geometry * sample)
	{
		this->sample = sample ;
		
		
		inclusionNumber = 20000 ;
		if(dmax == 0.004)
			inclusionNumber = 4096 ;
// 		for(int i = 2 ; i < 1000 ; i++)
// 		{
//			srand(0) ;
// 			massOfAggregates = 10.*area*fill*(double)i/1000 ; // .00000743*.75 ;
			if(dmax == 0.004)
				massOfAggregates = log(0.0000416) ;

			std::vector<double * > val ;
// 			val.push_back(&inclusionNumber);
			val.push_back(&massOfAggregates);
			
			std::vector<std::pair<double, double> > bounds ;
			
// 			bounds.push_back(std::make_pair(5000, 15000)) ;
			bounds.push_back(std::make_pair(log(.1e-5), log(1e-4))) ;
			GeneticAlgorithmOptimizer ga(val, bounds, this) ;
			ga.optimize(1e-4, 200, 20,  .1, .1) ;

			std::vector<Inclusion *> inclusions ;
			
			if(dmax == 0.004)
				inclusions = GranuloBolome(exp(massOfAggregates), 1, BOLOME_D)(dmax*.5, .000001, inclusionNumber, itzSize);
			else
				inclusions = GranuloBolome(exp(massOfAggregates), 1, BOLOME_A)(dmax*.5, .000001, inclusionNumber, itzSize);

			std::vector<Feature *> feats ;
			for(size_t i = 0; i < inclusions.size() ; i++)
				feats.push_back(inclusions[i]) ;

			int nAgg = 1 ;
			
			feats = placement(sample, feats, &nAgg, 6400, false);
			
			double volume = 0 ;
			for(size_t i = 0 ; i < feats.size() ; i++)
				volume += feats[i]->area() ;
			if(!feats.empty())
				std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() 
				<< ", smallest r =" << feats.back()->getRadius() 
				<< ", filling = " << volume/sample->area()*100.<< "%"<< std::endl ; 

			for(size_t i = 0 ; i < feats.size() ; i++)
				inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;

// 			if(!inclusions.empty())
// 			{
// 				std::cout << "largest inclusion with r = " << (*inclusions.begin())->getRadius() << std::endl ;
// 				std::cout << "smallest inclusion with r = " << (*inclusions.rbegin())->getRadius() << std::endl ;
// 			}
// 		}
		return feats ;
	}



	InclusionConverter::InclusionConverter(GeometryType type, RandomDistribution * a, RandomDistribution * ar, RandomDistribution * o) 
	{
		geom = type ;
		area = a ;
		aspectRatio = ar ;
		orientation = o ;
	}

	void InclusionConverter::setArea(RandomDistribution * a) 
	{
		area = a ;	
	}

	void InclusionConverter::setAspectRatio(RandomDistribution * ar) 
	{
		aspectRatio = ar ;
	}

	void InclusionConverter::setOrientation(RandomDistribution * o) 
	{
		orientation = o ;
	}

	void InclusionConverter::setArea(double a) 
	{
		area = new ConstantDistribution(a) ;
	}

	void InclusionConverter::setAspectRatio(double ar)
	{
		aspectRatio = new ConstantDistribution(ar) ;
	}

	void InclusionConverter::setOrientation(double o)
	{
		orientation = new ConstantDistribution(o) ;
	}

	Feature * InclusionConverter::convert(Inclusion * inc) const 
	{
		double newr, ar, a, b, o, ax, ay, bx, by, cx, cy ;
		Point center, A, B, C, D ;
			switch (geom) 
			{
				case CIRCLE:
					newr = inc->getRadius()*std::sqrt(area->draw()) ;
					return new Inclusion(inc->getFather(), newr, inc->getCenter()) ;
					
				case ELLIPSE:
					newr = inc->getRadius()*std::sqrt(area->draw()) ;
					
					ar = aspectRatio->draw() ;
					if(ar > 1)
					{
						a = newr*std::sqrt(ar) ;
						b = newr/std::sqrt(ar) ;
					}
					else
					{
						a = newr/std::sqrt(ar) ;
						b = newr*std::sqrt(ar) ;
					}
					
					o = orientation->draw() ;
					ax = a*std::cos(o) ;
					ay = a*std::sin(o) ;
					bx = b*std::sin(o) ;
					by = -b*std::cos(o) ;
					
					return new EllipsoidalInclusion(inc->getFather(), inc->getCenter(), Point(ax,ay), Point(bx,by)) ;
					
				case TRIANGLE:
					// produces rectangle triangle only !
					newr = inc->getRadius()*std::sqrt(2.*area->draw()*M_PI) ;
					
					ar = aspectRatio->draw() ;
					a = newr*std::sqrt(ar) ;
					b = newr/std::sqrt(ar) ;
					
					o = orientation->draw() ;
					bx = a*std::cos(o) ;
					by = a*std::sin(o) ;
					cx = b*std::sin(o) ;
					cy = -b*std::cos(o) ;
					
					A = Point(0,0) ;
					B = Point(bx,by) ;
					C = Point(cx,cy) ;
					
					center = (A+B+C)*0.3333333333333 ;
					center += inc->getCenter() ;
					
					return new TriangularInclusion(inc->getFather(), A-center, B-center, C-center) ;
					
				case RECTANGLE:
					newr = inc->getRadius()*std::sqrt(area->draw()*M_PI) ;
					
					ar = aspectRatio->draw() ;
					a = newr*std::sqrt(ar) ;
					b = newr/std::sqrt(ar) ;
					
					o = orientation->draw() ;
					bx = a*std::cos(o) ;
					by = a*std::sin(o) ;
					cx = b*std::sin(o) ;
					cy = -b*std::cos(o) ;
					
					A = Point(0,0) ;
					B = Point(bx,by) ;
					C = Point(cx,cy) ;
					D = Point(bx,cy) ;
					
					center = (A+B+C+D)*0.25 ;
					center += inc->getCenter() ;
					
					return new RectangularInclusion(inc->getFather(), A-center, B-center, C-center, D-center) ;
					
				case SPHERE:
					newr = inc->getRadius()*std::sqrt(area->draw()) ;
					return new Inclusion3D(inc->getFather(), newr, inc->getCenter()) ;
					
					
					
					
			}
			std::cout << "geometry type unsupported for inclusion translation" << std::endl ;
			return NULL ;
	}

	std::vector<Feature *> InclusionConverter::convert(std::vector<Inclusion *> inc) const 
	{
			std::vector<Feature *> ret ;
			for(size_t i = 0 ; i < inc.size() ; i++)
				ret.push_back(this->convert(inc[i])) ;
			return ret ;
	}

}