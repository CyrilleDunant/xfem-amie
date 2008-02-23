
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
	
	Grid grid(longueurX, longueurY, 100) ;
	
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

