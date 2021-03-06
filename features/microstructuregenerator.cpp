// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "microstructuregenerator.h"
#include "../physics/stiffness.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/optimizer.h"
#include <iostream>

namespace Amie
{
	
	AggregateDistribution2DGenerator::AggregateDistribution2DGenerator(double area, double dmax, double itzSize, double fill, double minMaxRatio) : area(area), dmax(dmax), fill(fill), minMaxRatio(minMaxRatio), itzSize(itzSize)
	{
		
	}
	
	double AggregateDistribution2DGenerator::score()
	{
//		srand(0);
		std::vector<Inclusion *> inclusions ;
		
		if(dmax == 0.004)
			inclusions = PSDGenerator::get2DInclusions(dmax*0.5, massOfAggregates, new PSDBolomeD(), PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;
		else
			inclusions = PSDGenerator::get2DInclusions(dmax*0.5, massOfAggregates, new PSDBolomeA(), PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;

		std::vector<Feature *> feats ;
		for(size_t i = 0; i < inclusions.size() ; i++)
			feats.push_back(inclusions[i]) ;

		feats = placement2D(sample, feats, itzSize, 0, 6400);

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
			
			for(size_t i = 0  ; i < inclusions.size() ; i++)
				delete inclusions[i] ;
			return 1./(.01*currentMinToMax/minMaxRatio + currentFill/fill) ;
		}
		else
		{
			for(size_t i = 0  ; i < inclusions.size() ; i++)
				delete inclusions[i] ;
			return 1000 ;
		}
	}

	void AggregateDistribution2DGenerator::print() const
	{
		std::vector<Inclusion *> inclusions ;
		
		if(dmax == 0.004)
			inclusions = PSDGenerator::get2DInclusions(dmax*0.5, massOfAggregates, new PSDBolomeD(), PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;
		else
			inclusions = PSDGenerator::get2DInclusions(dmax*0.5, massOfAggregates, new PSDBolomeA(), PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;

		std::vector<Feature *> feats ;
		for(size_t i = 0; i < inclusions.size() ; i++)
			feats.push_back(inclusions[i]) ;

		feats = placement2D(sample, feats, itzSize, 0, 6400);
		
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

		if(dmax == 0.004)
			massOfAggregates = log(0.0000416) ;

		std::vector<double * > val ;
		val.push_back(&massOfAggregates);
		
		std::vector<std::pair<double, double> > bounds ;
		
		bounds.push_back(std::make_pair(log(.1e-5), log(1e-4))) ;
		GeneticAlgorithmOptimizer ga(val, bounds, this) ;
		ga.optimize(1e-4, 200, 20,  .1, .1) ;

		std::vector<Inclusion *> inclusions ;
		
		if(dmax == 0.004)
			inclusions = PSDGenerator::get2DInclusions(dmax*0.5, massOfAggregates, new PSDBolomeD(), PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;
		else
			inclusions = PSDGenerator::get2DInclusions(dmax*0.5, massOfAggregates, new PSDBolomeD(), PSDEndCriteria(-1, 0.0001, inclusionNumber)) ;

		std::vector<Feature *> feats ;
		for(size_t i = 0; i < inclusions.size() ; i++)
			feats.push_back(inclusions[i]) ;

		feats = placement2D(sample, feats, itzSize, 0, 6400);
		
		double volume = 0 ;
		for(size_t i = 0 ; i < feats.size() ; i++)
			volume += feats[i]->area() ;
		if(!feats.empty())
			std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() 
			<< ", smallest r =" << feats.back()->getRadius() 
			<< ", filling = " << volume/sample->area()*100.<< "%"<< std::endl ; 

		for(size_t i = 0 ; i < feats.size() ; i++)
			inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;

		return feats ;
	}



std::vector<Feature *> InclusionGenerator::convert(std::vector<Inclusion *> inc)  
{
	std::vector<Feature *> ret ;
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		ret.push_back( this->convert( inc[i] ) ) ;
		if(i%100 == 0)
			std::cerr << "\rsetting up inclusion geometry: feature " << i << "/" << inc.size() << std::flush ;
	}
	std::cerr << "\rsetting up inclusion geometry: feature " << inc.size() << "/" << inc.size() << std::endl ;
	return ret ;
}

Feature * EllipsoidalInclusionGenerator::convert(Inclusion * inc)  
{
	double r = inc->getRadius() ;
	double aspect = shape(rng) ;
	if(aspect < 0.1)
		aspect = 0.1 ;
	if(aspect > 1-POINT_TOLERANCE)
	{
		EllipsoidalInclusion * ret = new EllipsoidalInclusion( inc->getFather(), inc->getCenter(), Point(r,0), Point(0,r) ) ;
		ret->setBehaviour( inc->getBehaviour()) ;
		return ret ;
	}
	double a = r/std::sqrt(aspect) ;
	double b = r*std::sqrt(aspect) ;

	double phase = orientation(rng) ;
	EllipsoidalInclusion * ret = new EllipsoidalInclusion(inc->getFather(), inc->getCenter(), Point(a*cos( phase ), a*sin(phase)), Point( b*(-sin(phase)), b*cos(phase)) ) ;
	ret->setBehaviour( inc->getBehaviour()) ;
	return ret ;
}

Feature * RectangularInclusionGenerator::convert(Inclusion * inc)  
{
	double r = inc->getRadius() ;

	double aspect = shape(rng) ;
	if(aspect < 0.1)
		aspect = 0.1 ;

	double a = r/std::sqrt(aspect) ;
	double b = r*std::sqrt(aspect) ;

	double direction = orientation(rng) ;

	Point A( 0,0 ) ;
	Point B( a*cos(direction), a*sin(direction) ) ;
	Point C( b*(-sin(direction)), b*cos(direction) ) ;
	Point D = B+C ;

	Point center = (A+B+C+D)*0.25 ;
	center -= inc->getCenter() ;
					
	RectangularInclusion * ret = new RectangularInclusion(inc->getFather(), A-center, C-center, D-center, B-center) ;
	ret->setBehaviour( inc->getBehaviour()) ;
	return ret ;
}

Feature * XFEMInclusionGenerator::convert(Inclusion * inc)  
{
	return new ExpansiveZone( inc->getFather(), inc->getRadius(), inc->getCenter().getX(), inc->getCenter().getY(), inc->getBehaviour() ) ;
}

Feature * SpaceTimeXFEMInclusionGenerator::convert(Inclusion * inc)  
{
	Function r( inc->getRadius() ) ;
	return new GrowingExpansiveZone( inc->getFather(), r, inc->getCenter().getX(), inc->getCenter().getY() ) ;
}

PolygonalSample * PolygonalInclusionGenerator::generatePolygon(double radius) 
{
	size_t npoints = vertex(rng) ;
	double phase = orientation(rng) ;
	std::valarray< Point *> points ; points.resize( npoints) ;
	for(size_t i = 0 ; i < npoints ; i++)
	{
		double theta = phase + 2*M_PI*(double) i/(double) npoints ;
		points[i] = new Point( radius*cos(theta+phase), radius*sin(theta+phase) ) ;
	}

	return new PolygonalSample( nullptr, points ) ;
}

Feature * PolygonalInclusionGenerator::convert( Inclusion * inc)  
{
	PolygonalSample * ret = this->generatePolygon( inc->getRadius() ) ;
	double area = ret->area() ;
	double target = inc->area() ;
	double r = ret->getRadius() * sqrt( target / area ) ;
	transform( dynamic_cast<Polygon *>(ret), SCALE, Point( r/ret->getRadius(), r/ret->getRadius()) ) ;
	if(forceOrientation)
	{
		Point p = ret->getOrientation() ;
		double phase = p.angle() - orientation(rng)  ;
		if(std::abs(phase) > POINT_TOLERANCE)
			transform(  dynamic_cast<Polygon *>(ret), ROTATE, Point( 0,0, phase ) ) ;
	}
	ret->setCenter(inc->getCenter()) ;
	ret->setBehaviour( inc->getBehaviour() ) ;
	ret->setFather( inc->getFather() ) ;
	return ret ;
}

VoronoiPolygonalInclusionGenerator::VoronoiPolygonalInclusionGenerator( double box, size_t seed, double minDist, double s, double o, double ov, double sv, double rot, bool force) : PolygonalInclusionGenerator(5, o, ov, 0, rot, force)
{
	Rectangle rect( box, box, 0,0 ) ;
	std::vector<VoronoiGrain> morphology ;
	morphology.push_back( VoronoiGrain( nullptr, minDist, 1., seed, 1.) ) ;
	source = PSDGenerator::get2DSourceVoronoiPolygons(&rect, morphology, seed, minDist)[0] ;
        index = std::uniform_int_distribution< size_t >(0,source.size()-1) ;
	if(source.size() == 0)
	{
		std::cout << "no polygons available after Voronoi tesselation, exiting now" << std::endl ;
		exit(0) ;
	}
	std::uniform_real_distribution< double > shape(s-sv,s+sv) ;
        for(size_t i = 0 ; i < source.size() ; i++)
	{
		std::valarray<Point> pts = source[i]->getOriginalPoints() ;
		for(size_t j = 0 ; j < pts.size() ; j++)
			pts[j].getX() *= shape( rng ) ;
		source[i]->setOriginalPoints( pts, true ) ;
	}
}

PolygonalSample * VoronoiPolygonalInclusionGenerator::generatePolygon(double radius) 
{
	size_t i = index( rng ) ;
	std::valarray<Point> opts = source[i]->getOriginalPoints() ;
	std::valarray<Point *> pts( opts.size() ) ;
	for(size_t j = 0 ; j < opts.size() ; j++)
		pts[j] = new Point( opts[j] ) ;
	PolygonalSample * ret = new PolygonalSample( nullptr, pts ) ;
	ret->setCenter( Point(0,0) ) ;
	return ret ;
}

PolygonalSample * GravelPolygonalInclusionGenerator::generatePolygon(double radius) 
{
	size_t npoints = vertex(rng) ;
	double phase = orientation( rng ) ;
	std::valarray< Point *> points ; points.resize( npoints) ;
	Vector A ; A.resize(m) ; A = 0. ;
	Vector alpha ; alpha.resize(m) ; alpha = 0. ;
	for(size_t i = 0 ; i < m ; i++)
	{
		A[i] = radius*exp(-p*log(i+1)-b) ;
		alpha[i] = phi(rng) ;
	}
	for(size_t i = 0 ; i < npoints ; i++)
	{
		double theta = phase + 2*M_PI*(double) i/(double) npoints ;
		double r = radius ;
		for(size_t j = 0 ; j < m ; j++)
			r += A[j]*cos( (j+1)*theta + alpha[j]) ;
		points[i] = new Point( r*cos(theta), r*sin(theta) ) ;
	}

	return new PolygonalSample( nullptr, points ) ;
}

PolygonalSample * CrushedPolygonalInclusionGenerator::generatePolygon(double radius) 
{
	size_t npoints = vertex(rng) ;
	double phase = orientation(rng) ;
	std::valarray< Point *> points ;points.resize( npoints) ;
	std::vector<double> theta ;
	for(size_t i = 0 ; i < npoints ; i++)
		theta.push_back( phi(rng) ) ;
	std::sort( theta.begin(), theta.end() ) ;
	for(size_t i = 0 ; i < npoints ; i++)
	{
		double r = radius  * (1.+shape(rng)) ;
		points[i] = new Point( r*cos(theta[i]+phase), r*sin(theta[i]+phase) ) ;
	}

	return new PolygonalSample( nullptr, points ) ;
}

PolygonalSample * CrushedSubtendedPolygonalInclusionGenerator::generatePolygon(double radius) 
{
	size_t npoints = vertex(rng) ;
	double phase = orientation(rng) ;
	std::valarray< Point *> points ;points.resize( npoints) ;
	std::vector<double> phi ;
	double beta = 2.*M_PI/(double) npoints ;
	double sumphi = 0. ;
	for(size_t i = 0 ; i < npoints ; i++)
	{
		phi.push_back( beta * ( 1.+ delta(rng) ) ) ;
		sumphi += phi[i] ;
	}
	for(size_t i = 0 ; i < phi.size() ; i++)
	{
		phi[i] *= 2.*M_PI/sumphi ;
	}
	double theta = phase ;
	for(size_t i = 0 ; i < npoints ; i++)
	{
		double r = radius * (1.+shape(rng)) ;
		points[i] = new Point( r*cos(theta), r*sin(theta) ) ;
		theta += phi[i] ;
	}

	return new PolygonalSample( nullptr, points ) ;
}

}


