
// Author: Jérôme Krebs <jerome.krebs@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "placement.h"
#include "random.h"
#include "../geometry/geometry_base.h"


using namespace Amie ;


double Amie::chiffreAleatoire(double longueur) // fonction qui retourne une valeur aléatoire
{
	double chiffreAleatoire = longueur*(double)rand()/(double)RAND_MAX;
	return chiffreAleatoire;
}



bool bord(double r, double longueurX, double longueurY, double x, double y)//fonction du problème de bord
{
	if(x+r > longueurX/2. 
	|| x-r < -longueurX/2.
	|| y+r > longueurY/2.
	|| y-r < -longueurY/2.
	)
	{			
		return true;	
	}
	return false;
}

bool intersections( Feature * feat, const std::vector<Geometry *> & exclusionZones)
{
	for(size_t i = 0 ; i < exclusionZones.size() ; i++)
	{
		if(feat->intersects(exclusionZones[i]))
			return true ;
	}
	return false ;
}

bool isInside( Feature * feat, const std::vector<Geometry *> & base)
{
	for(size_t j = 0 ; j < base.size() ; j++)
	{
		if(feat->intersects(base[j]))
		{
			Point p = feat->getCenter() ;
			base[j]->project(&p) ;
			return !feat->in(p) ;
		}
	}
	return false ;
}


void transform2D( Feature * inc, RandomDistribution & xDistribution, RandomDistribution & yDistribution, RandomDistribution & rDistribution)
{
	Point c( xDistribution.draw(), yDistribution.draw() ) ;
	Point theta( 0,0, rDistribution.draw() ) ;
	inc->setCenter( c ) ;
	inc->transform(ROTATE, theta) ;
}

std::vector<Feature *> Amie::placement2D(const Geometry* box, std::vector<Feature *> inclusions, double minDist, int placedAggregates, int triesMax, double orientation,  std::vector<Geometry *> exclusionZones) 
{
	std::vector<Feature *> ret ;
	int tries = 0 ;
	
	std::vector<Point> boundingBox = box->getBoundingBox() ;
	UniformDistribution xDistribution( boundingBox[0].getX(), boundingBox[2].getX() ) ;
	UniformDistribution yDistribution( boundingBox[0].getY(), boundingBox[2].getY() ) ;
	UniformDistribution rDistribution( -orientation, orientation ) ;
	Grid grid(boundingBox[2].getX()-boundingBox[0].getX(), boundingBox[0].getY()-boundingBox[2].getY(), 10, box->getCenter()) ;

	for(size_t i = 0 ; i < placedAggregates ; i++)
	{
		ret.push_back(inclusions[i]);
		grid.add(inclusions[i]) ;
	}

	for(size_t i = placedAggregates ; i < inclusions.size() && tries < triesMax ; i++)
	{
		tries++ ;
		
		double scale = 1. ;
		if(minDist > POINT_TOLERANCE_2D)
		{
			scale = (inclusions[i]->getRadius()+minDist)/inclusions[i]->getRadius() ;
			Point s( scale, scale ) ;
			inclusions[i]->transform(SCALE, s) ;
		}
		
		transform2D( inclusions[i], xDistribution, yDistribution, rDistribution); 
		std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
		while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || intersections(inclusions[i], exclusionZones) && tries < triesMax)
		{
			tries++ ;
			transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);  
			bbox = inclusions[i]->getBoundingBox() ;
		}
		
		while(!grid.add(inclusions[i]) && tries < triesMax)
		{
			tries++ ;
			
			transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);  
			std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
			while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))|| intersections(inclusions[i], exclusionZones) && tries < triesMax )
			{
				tries++ ;
				transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);  
				bbox = inclusions[i]->getBoundingBox() ;
			}
			
		}
		
		if(tries < triesMax)
		{
			if(i%100 == 0)
				std::cout << "\rplaced " << i << " particles (tries " << tries << "/" << triesMax << ")" << std::flush ;
			if(scale > 1.)
			{
				Point s(1./scale, 1./scale) ;
				inclusions[i]->transform( SCALE , s) ;
			}
			ret.push_back(inclusions[i]);
		}
		  
	}
	
	std::cout << "\n" << ret.size() << " inclusions placed after " << tries << " tries" << std::endl ;
	
	return ret ;
}

std::vector<Feature *> Amie::placement2DInInclusions(const Geometry* box, std::vector<Geometry *> base, std::vector<Feature *> inclusions, double minDist, int placedAggregates, int triesMax, double orientation,  std::vector<Geometry *> exclusionZones) 
{
	std::vector<Feature *> ret ;
	int tries = 0 ;
	
	std::vector<Point> boundingBox = box->getBoundingBox() ;
	UniformDistribution xDistribution( boundingBox[0].getX(), boundingBox[2].getX() ) ;
	UniformDistribution yDistribution( boundingBox[0].getY(), boundingBox[2].getY() ) ;
	UniformDistribution rDistribution( -orientation, orientation ) ;
	Grid grid(boundingBox[2].getX()-boundingBox[0].getX(), boundingBox[0].getY()-boundingBox[2].getY(), 10, box->getCenter()) ;

	for(size_t i = 0 ; i < placedAggregates ; i++)
	{
		ret.push_back(inclusions[i]);
		grid.add(inclusions[i]) ;
	}

	for(size_t i = placedAggregates ; i < inclusions.size() && tries < triesMax ; i++)
	{
		tries++ ;
		
		double scale = 1. ;
		if(minDist > POINT_TOLERANCE_2D)
		{
			scale = (inclusions[i]->getRadius()+minDist)/inclusions[i]->getRadius() ;
			Point s( scale, scale ) ;
			inclusions[i]->transform(SCALE, s) ;
		}
		
		transform2D( inclusions[i], xDistribution, yDistribution, rDistribution); 
		std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
		while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || !isInside(inclusions[i], base) || intersections(inclusions[i], exclusionZones) && tries < triesMax)
		{
			tries++ ;
			transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);  
			bbox = inclusions[i]->getBoundingBox() ;
		}
		
		while(!grid.add(inclusions[i]) && tries < triesMax)
		{
			tries++ ;
			
			transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);  
			std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
			while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3]))  || !isInside(inclusions[i], base) || intersections(inclusions[i], exclusionZones) && tries < triesMax)
			{
				tries++ ;
				transform2D( inclusions[i], xDistribution, yDistribution, rDistribution);  
				bbox = inclusions[i]->getBoundingBox() ;
			}
			
		}
		
		if(tries < triesMax)
		{
			if(i%100 == 0)
				std::cout << "\rplaced " << i << " particles (tries " << tries << "/" << triesMax << ")" << std::flush ;
			if(scale > 1.)
			{
				Point s(1./scale, 1./scale) ;
				inclusions[i]->transform( SCALE , s) ;
			}
			ret.push_back(inclusions[i]);
		}
		

	}
	
	std::cout << "\n" << ret.size() << " inclusions placed after " << tries << " tries" << std::endl ;
	
	return ret ;
}


std::vector<Feature *> Amie::placement(const Geometry * box, std::vector<Feature *> inclusions, int *nombreGranulatsPlaces, int nombreGranulatsDejaPlaces, int triesMax, std::vector<Geometry *> exclusionZones, bool verbose)
{
	int tries = 0 ;
	
	RandomNumber gen ;

	double volume = 0 ;
	std::vector<Feature *> ret ;
	
	if(box->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		if(verbose)
			std::cout << "placing..." << std::endl ;
		Point offset = box->getCenter() ;
		std::vector<Point> boundingBox = box->getBoundingBox() ;
		double longueurX = std::abs(boundingBox[2].getX()-boundingBox[0].getX());
		double longueurY = std::abs(boundingBox[0].getY()-boundingBox[2].getY());
		if(verbose)
			std::cout << longueurX << ", " << longueurY << std::endl ;
		Grid grid(longueurX, longueurY, 10, box->getCenter()) ;
		longueurX*=1.2 ;
		longueurY*=1.2 ;
		for(size_t i = 0 ; i < nombreGranulatsDejaPlaces ; i++)
		{
			grid.add(inclusions[i]) ;
			ret.push_back(inclusions[i]) ;
		}
		
		for(size_t i=nombreGranulatsDejaPlaces ; i < inclusions.size() && tries < triesMax ; i++) 
		{
			tries++ ;
			double ix = longueurX - 2.1*inclusions[i]->getRadius() ;
			double iy = longueurY - 2.1*inclusions[i]->getRadius() ;
			Point newCentre(inclusions[i]->getCenter()) ;
			if(i >= nombreGranulatsDejaPlaces)
			{
				newCentre.getX() = gen.uniform(ix) - ix/2. + offset.getX() ;
				newCentre.getY() = gen.uniform(iy) - iy/2. + offset.getY() ;
			}
			inclusions[i]->setCenter(newCentre) ;
			std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
			
			while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || intersections(inclusions[i], exclusionZones) )/*|| inclusions[0]->in(inclusions[i]->getCenter()) || (inclusions[0]->in(bbox[0]) || inclusions[0]->in(bbox[1]) || inclusions[0]->in(bbox[2]) || inclusions[0]->in(bbox[3])))*/
			{
				Point newCentre(gen.uniform(ix) - ix/2. + offset.getX() , gen.uniform(iy) - iy/2. + offset.getY()) ;
				inclusions[i]->setCenter(newCentre) ;
				bbox = inclusions[i]->getBoundingBox() ;
			}

			while(!grid.add(inclusions[i]) && tries < triesMax)
			{
				tries++ ;
				Point newCentre(gen.uniform(ix) - ix/2. + offset.getX() , gen.uniform(iy) - iy/2. + offset.getY()) ;
				inclusions[i]->setCenter(newCentre) ;
				bbox = inclusions[i]->getBoundingBox() ;
				while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) || intersections(inclusions[i], exclusionZones))/* || inclusions[0]->in(inclusions[i]->getCenter()) || (inclusions[0]->in(bbox[0]) || inclusions[0]->in(bbox[1]) || inclusions[0]->in(bbox[2]) || inclusions[0]->in(bbox[3])))*/
				{
					Point newCentre(gen.uniform(ix) - ix/2. + offset.getX() , gen.uniform(iy) - iy/2. + offset.getY()) ;
					inclusions[i]->setCenter(newCentre) ;
					bbox = inclusions[i]->getBoundingBox() ;
				}
			}

// 			if(tries == triesMax)
// 				std::cout << " triesmax" << std::endl ;
			
			if(tries< triesMax)
			{
				if(verbose)
					if(i%100 == 0)
						std::cout << "\rplaced " << i << " particles" << std::flush ;
				ret.push_back(inclusions[i]) ;
				volume += inclusions[i]->area() ;
				tries = 0 ;
			}
			else
				break ;
		}
		
		if(verbose)
			std::cout << "\n placed aggregate volume = " << volume << std::endl ;
		
		return ret ;
	}
	else
	{
		Point offset = box->getCenter() ;
		std::vector<Point> boundingBox = box->getBoundingBox() ;
		double longueurX = std::abs(boundingBox[0].getX()-boundingBox[7].getX());
		double longueurY = std::abs(boundingBox[0].getY()-boundingBox[7].getY());
		double longueurZ = std::abs(boundingBox[0].getZ()-boundingBox[7].getZ());
		double ndiv = 1 ; //round(((longueurX+longueurY+longueurZ)/3.)/(inclusions[inclusions.size()/2]->getRadius()*2.)) ;

		Grid3D *grid = new Grid3D(longueurX, longueurY, longueurZ, ndiv, offset) ;
		
		for(size_t i=nombreGranulatsDejaPlaces ; i < inclusions.size() && tries < triesMax ; i++) 
		{
			tries++ ;
			Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(), 
			                chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY(),
			                chiffreAleatoire(longueurZ-2.1*inclusions[i]->getRadius())-(longueurZ-2.1*inclusions[i]->getRadius())/2. + offset.getZ()
			               ) ;
			inclusions[i]->setCenter(newCentre) ;
			while(!box->in(inclusions[i]->getCenter()) && !box->intersects(inclusions[i]) )
			{
				Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(), 
				                chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY(),
				                chiffreAleatoire(longueurZ-2.1*inclusions[i]->getRadius())-(longueurZ-2.1*inclusions[i]->getRadius())/2. + offset.getZ()
				               ) ;
				inclusions[i]->setCenter(newCentre) ;
			}
			
			while(!grid->add(inclusions[i]) && tries < triesMax)
			{
				tries++ ;
				Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(), 
				                chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY(),
				                chiffreAleatoire(longueurZ-2.1*inclusions[i]->getRadius())-(longueurZ-2.1*inclusions[i]->getRadius())/2. + offset.getZ()
				               ) ;
				inclusions[i]->setCenter(newCentre) ;
				while(!box->in(inclusions[i]->getCenter()) && !box->intersects(inclusions[i]) )
				{
					Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(), 
					                chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY(),
					                chiffreAleatoire(longueurZ-2.1*inclusions[i]->getRadius())-(longueurZ-2.1*inclusions[i]->getRadius())/2. + offset.getZ()
					               ) ;
					inclusions[i]->setCenter(newCentre) ;
				}
			}

// 		for(size_t i=0 ; i < inclusions.size() && tries < triesMax ; i++) 
// 		{
// 			tries++ ;
// 			inclusions[i]->setCenter( grid->randomFreeCenter()) ;
// 			while(!box->in(inclusions[i]->getCenter()) || box->intersects(inclusions[i]) )
// 			{
// 				inclusions[i]->setCenter(grid->randomFreeCenter()) ;
// 			}
// 			bool go =true ;
// 			for(size_t j = 0 ; j < i ; j++)
// 			{
// 				go = !(inclusions[i]->intersects(inclusions[j])) ;
// 				if(!go)
// 				  break ;
// 			}
// 			while(!go && tries < triesMax)
// 			{
// 				tries++ ;
// 				inclusions[i]->setCenter(grid->randomFreeCenter()) ;
// 				while(!box->in(inclusions[i]->getCenter()) || box->intersects(inclusions[i]) )
// 				{
// 					inclusions[i]->setCenter(grid->randomFreeCenter()) ;
// 				}
// 				
// 				go =true ;
// 				for(size_t j = 0 ; j < i ; j++)
// 				{
// 					go = !(inclusions[i]->intersects(inclusions[j])) ;
// 					if(!go)
// 					  break ;
// 				}
// 			}
			
			if(tries< triesMax)
			{
				if(verbose)
					if(i%100 == 0)
						std::cout << "\rplaced " << i << " particles" << std::flush ;
				ret.push_back(inclusions[i]) ;
				volume += inclusions[i]->volume() ;
				tries = 0 ;
			}
			else
			{
				if(verbose)
					std::cout << "\rplaced " << i << " particles" << std::flush ;
				break ;			
			}
			}
		
		if(verbose)
			std::cout << "\n placed aggregate volume = " << volume << std::endl ;
		delete grid ;
		return ret ;
		
	}
		
}

std::vector<Amie::EllipsoidalInclusion *> Amie::placement_with_rotation(const Geometry * box, std::vector<EllipsoidalInclusion *> inclusions, int *nombreGranulatsPlaces, int triesMax, bool verbose) 
{
	int tries = 0 ;
	int changeAxis = 0 ;
	int hasBeenChanged = 0 ;

	double volume = 0 ;
	std::vector<EllipsoidalInclusion *> ret ;
	
	if(verbose)
		std::cout << "placing..." << std::endl ;
	Point offset = box->getCenter() ;
	std::vector<Point> boundingBox = box->getBoundingBox() ;
	double longueurX = std::abs(boundingBox[2].getX()-boundingBox[0].getX());
	double longueurY = std::abs(boundingBox[0].getY()-boundingBox[2].getY());
	if(verbose)
		std::cout << longueurX << ", " << longueurY << std::endl ;
	Grid grid(longueurX, longueurY, 10, box->getCenter()) ;
	longueurX*=1.2 ;
	longueurY*=1.2 ;
	for(size_t i=0 ; i < inclusions.size() && tries < triesMax ; i++) 
	{
		tries++ ;
		Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(), chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY()) ;
		inclusions[i]->Ellipse::setCenter(newCentre) ;
		std::vector<Point> bbox = inclusions[i]->getBoundingBox() ;
		while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) )
		{
			Point newCentre(
			                 chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())
			                 -(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(),
			                 chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())
			                 -(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY()
			               ) ;
			inclusions[i]->Ellipse::setCenter(newCentre) ;
			bbox = inclusions[i]->getBoundingBox() ;
		}
		while(!grid.add(inclusions[i]->getPrimitive()) && tries < triesMax)
		{
			tries++ ;
			changeAxis++ ;
			Point newCentre(
			                 chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())
			                 - (longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(),
			                 chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())
			                 - (longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY()
			               ) ;
			inclusions[i]->Ellipse::setCenter(newCentre) ;
/*			if(changeAxis > triesMax/5)
			{
				Point newAxis((double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX) ;
				inclusions[i] = new EllipsoidalInclusion(inclusions[i]->getMajorRadius(),
									inclusions[i]->getMinorRadius(),
									newCentre,newAxis) ;
				changeAxis = 0 ;
				std::cout <<i << " => Axis has been changed" << std::endl ;
				hasBeenChanged++ ;
			}*/
			bbox = inclusions[i]->getBoundingBox() ;
			while(!box->in(inclusions[i]->getCenter()) || !(box->in(bbox[0]) && box->in(bbox[1]) && box->in(bbox[2]) && box->in(bbox[3])) )
			{
				changeAxis++ ;
				Point newCentre(
				                 chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())
				                 - (longueurX-2.1*inclusions[i]->getRadius())/2. + offset.getX(),
				                 chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())
				                 - (longueurY-2.1*inclusions[i]->getRadius())/2. + offset.getY()
				               ) ;
				inclusions[i]->Ellipse::setCenter(newCentre) ;
/*				if(changeAxis > triesMax/5)
				{
					Point newAxis((double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX) ;
					inclusions[i] = new EllipsoidalInclusion(inclusions[i]->getMajorRadius(),
										inclusions[i]->getMinorRadius(),
										newCentre,newAxis) ;
					changeAxis = 0 ;
					std::cout << i << " => Axis has been changed" << std::endl ;
					hasBeenChanged++ ;
				}*/
				bbox = inclusions[i]->getBoundingBox() ;
			}
		}
		
		if(tries< triesMax)
		{
			if(verbose)
				if(i%100 == 0)
					std::cout << "\rplaced " << i << " particles" << std::flush ;
			ret.push_back(inclusions[i]) ;
			volume += inclusions[i]->area() ;
			tries = 0 ;
			changeAxis = 0 ;
		}
		else
			break ;
	}
	
	if(verbose)
		std::cout << "\n placed aggregate volume = " << volume << std::endl ;
	if(verbose)
		std::cout << hasBeenChanged << " aggregates have been rotated from initial major axis" << std::endl ;
	
	return ret ;
}




