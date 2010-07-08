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
		srand(0);
		std::vector<Inclusion *> inclusions ;
		
		if(dmax == 0.004)
			inclusions = GranuloBolome(massOfAggregates, 1, BOLOME_D)(dmax*.5, .0001, inclusionNumber, itzSize, false);
		else
			inclusions = GranuloBolome(massOfAggregates, 1, BOLOME_A)(dmax*.5, .0001, inclusionNumber, itzSize, false);

		std::vector<Feature *> feats ;
		for(size_t i = 0; i < inclusions.size() ; i++)
			feats.push_back(inclusions[i]) ;

		int nAgg = 1 ;
		feats = placement(sample, feats, &nAgg, 6400, false);
		
		double volume = 0 ;
		for(size_t i = 0 ; i < feats.size() ; i++)
			volume += feats[i]->area() ;
		if(!feats.empty())
		{
			std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() 
			<< ", smallest r =" << feats.back()->getRadius() << ", ratio = " << feats.front()->getRadius()/feats.back()->getRadius()
			<< ", filling = " << volume/sample->area()*100.<< "%"<< std::endl ; 
		

			double currentMinToMax = feats.front()->getRadius()/feats.back()->getRadius() ;
			double currentFill = volume/sample->area() ;
			
			for(int i = 0  ; i < inclusions.size() ; i++)
				delete inclusions[i] ;
			return 1./(currentMinToMax/minMaxRatio + currentFill/fill) ;
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
		
		
		inclusionNumber = 8000 ;
		if(dmax == 0.004)
			inclusionNumber = 4096 ;
// 		for(int i = 2 ; i < 1000 ; i++)
// 		{
			srand(0) ;
// 			massOfAggregates = 10.*area*fill*(double)i/1000 ; // .00000743*.75 ;
			if(dmax == 0.004)
				massOfAggregates = 0.0000416 ;

			std::vector<double * > val ;
			val.push_back(&inclusionNumber);
			val.push_back(&massOfAggregates);
			
			std::vector<std::pair<double, double> > bounds ;
			
			bounds.push_back(std::make_pair(5000, 15000)) ;
			bounds.push_back(std::make_pair(6.9151e-06/16, 6.9151e-06*16)) ;
			GeneticAlgorithmOptimizer ga(val, bounds, this) ;
			ga.optimize(1e-4, 20, 10,  .1, .1) ;

			std::vector<Inclusion *> inclusions ;
			
			if(dmax == 0.004)
				inclusions = GranuloBolome(massOfAggregates, 1, BOLOME_D)(dmax*.5, .000001, inclusionNumber, itzSize);
			else
				inclusions = GranuloBolome(massOfAggregates, 1, BOLOME_A)(dmax*.5, .000001, inclusionNumber, itzSize);

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

}