
// Author:  Jérôme Krebs <jerome.krebs@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//

#include "placement.h"


using namespace Mu ;


double Mu::chiffreAleatoire(double longueur) // fonction qui retourne une valeur aléatoire
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

std::vector<Inclusion *> Mu::placement(double longueurX, double longueurY, std::vector<Inclusion *> inclusions, int *nombreGranulatsPlaces, int triesMax)
{
	int tries = 0 ;

	double volume = 0 ;
	std::vector<Inclusion *> ret ;
	
	Grid grid(longueurX, longueurY, 20) ;
	
	for(size_t i=0 ; i < inclusions.size() && tries < triesMax ; i++) 
	{
		tries++ ;

		inclusions[i]->getCenter().x = chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2.;
		inclusions[i]->getCenter().y = chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2.;
		while(!grid.add(inclusions[i]) && tries < triesMax)
		{
			tries++ ;
			inclusions[i]->getCenter().x = chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2.;
			inclusions[i]->getCenter().y = chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2.;
		}
		
		if(tries< triesMax)
		{
			if(i%100 == 0)
				std::cout << "\rplaced " << i << " particles" << std::flush ;
			ret.push_back(inclusions[i]) ;
			volume += inclusions[i]->area() ;
			tries = 0 ;
		}
		else
			break ;
	}
	
	std::cout << "\n placed aggregate volume = " << volume << std::endl ;
	
	return ret ;
		
}

std::vector<Feature *> Mu::placement(const Geometry * box, std::vector<Feature *> inclusions, int *nombreGranulatsPlaces, int triesMax)
{
	int tries = 0 ;

	double volume = 0 ;
	std::vector<Feature *> ret ;
	
	if(box->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		std::cout << "placing..." << std::endl ;
		Point offset = box->getCenter() ;
		std::vector<Point> boundingBox = box->getBoundingBox() ;
		double longueurX = std::abs(boundingBox[2].x-boundingBox[0].x);
		double longueurY = std::abs(boundingBox[0].y-boundingBox[2].y);
		std::cout << longueurX << ", " << longueurY << std::endl ;
		Grid grid(longueurX, longueurY, 1) ;
		longueurX*=1.2 ;
		longueurY*=1.2 ;
		for(size_t i=0 ; i < inclusions.size() && tries < triesMax ; i++) 
		{
			tries++ ;
			Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.x, chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.y) ;
			inclusions[i]->setCenter(newCentre) ;
			while(!box->in(inclusions[i]->getCenter()) || !box->intersection(inclusions[i]).empty() )
			{
				Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.x, chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.y) ;
				inclusions[i]->setCenter(newCentre) ;
			}

			inclusions[i]->setCenter(inclusions[i]->getCenter()- offset) ;
			
			while(!grid.add(inclusions[i]) && tries < triesMax)
			{
				tries++ ;
				Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.x, chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.y) ;
				inclusions[i]->setCenter(newCentre) ;
				while(!box->in(inclusions[i]->getCenter()) || !box->intersection(inclusions[i]).empty() )
				{
					Point newCentre(chiffreAleatoire(longueurX-2.1*inclusions[i]->getRadius())-(longueurX-2.1*inclusions[i]->getRadius())/2. + offset.x, chiffreAleatoire(longueurY-2.1*inclusions[i]->getRadius())-(longueurY-2.1*inclusions[i]->getRadius())/2. + offset.y) ;
					inclusions[i]->setCenter(newCentre) ;
				}
			}
			
			if(tries< triesMax)
			{
				if(i%100 == 0)
					std::cout << "\rplaced " << i << " particles" << std::flush ;
				inclusions[i]->setCenter(inclusions[i]->getCenter() + offset) ;
				ret.push_back(inclusions[i]) ;
				volume += inclusions[i]->area() ;
				tries = 0 ;
			}
			else
				break ;
		}
		
		std::cout << "\n placed aggregate volume = " << volume << std::endl ;
		
		return ret ;
	}
	else
	{
		Point offset = box->getCenter() ;
		std::vector<Point> boundingBox = box->getBoundingBox() ;
		double longueurX = std::abs(boundingBox[0].x-boundingBox[7].x);
		double longueurY = std::abs(boundingBox[0].y-boundingBox[7].y);
		double longueurZ = std::abs(boundingBox[0].z-boundingBox[7].z);
		double ndiv = 1 ;
		Grid3D *grid = new Grid3D(longueurX, longueurY, longueurZ, ndiv) ;
		longueurX*=1.2 ;
		longueurY*=1.2 ;
		longueurZ*=1.2 ;
		for(size_t i=0 ; i < inclusions.size() && tries < triesMax ; i++) 
		{
			tries++ ;
			inclusions[i]->setCenter( grid->randomFreeCenter()) ;
			while(!box->in(inclusions[i]->getCenter()) || box->intersects(inclusions[i]) )
			{
				inclusions[i]->setCenter(grid->randomFreeCenter()) ;
			}
			
			inclusions[i]->setCenter(inclusions[i]->getCenter()- offset) ;
			
			while(!grid->add(inclusions[i]) && tries < triesMax)
			{
				tries++ ;
				inclusions[i]->setCenter(grid->randomFreeCenter()) ;
				while(!box->in(inclusions[i]->getCenter()) || box->intersects(inclusions[i]) )
				{
					inclusions[i]->setCenter(grid->randomFreeCenter()) ;
				}
			}
			
			if(tries< triesMax)
			{
				if(i%100 == 0)
					std::cout << "\rplaced " << i << " particles" << std::flush ;
				inclusions[i]->setCenter(inclusions[i]->getCenter() + offset) ;
				ret.push_back(inclusions[i]) ;
				volume += inclusions[i]->volume() ;
				tries = 0 ;
			}
			else
				break ;
				
		}
		
		std::cout << "\n placed aggregate volume = " << volume << std::endl ;
		delete grid ;
		return ret ;
	}
		
}

