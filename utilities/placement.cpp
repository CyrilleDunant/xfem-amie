
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
	
	Grid grid(longueurX, longueurY, 40) ;
	
	for(size_t i=0 ; i < inclusions.size() && tries < triesMax ; i++) 
	{
		tries++ ;

		inclusions[i]->getCenter().x = chiffreAleatoire(longueurX-2.*inclusions[i]->getRadius())-(longueurX-2.*inclusions[i]->getRadius())/2.;
		inclusions[i]->getCenter().y = chiffreAleatoire(longueurY-2.*inclusions[i]->getRadius())-(longueurY-2.*inclusions[i]->getRadius())/2.;
		while(!grid.add(inclusions[i]) && tries < triesMax)
		{
			tries++ ;
			inclusions[i]->getCenter().x = chiffreAleatoire(longueurX-2.*inclusions[i]->getRadius())-(longueurX-2.*inclusions[i]->getRadius())/2.;
			inclusions[i]->getCenter().y = chiffreAleatoire(longueurY-2.*inclusions[i]->getRadius())-(longueurY-2.*inclusions[i]->getRadius())/2.;
		}
		
		if(tries< triesMax)
		{
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

Pixel::Pixel()
{
	filled = false ;
}

Pixel::Pixel(double x, double y, double s) : tl(x-s*.5, y+s*.5), tr(x+s*.5, y+s*.5), bl(x-s*.5, y-s*.5), br(x+s*.5, y-s*.5), filled(false) { } ;

const std::vector<Mu::Inclusion *> & Pixel::getInclusions() const
{
	return inclusions ;
}

std::vector<Mu::Inclusion *> & Pixel::getInclusions()
{
	return inclusions ;
}

bool Pixel::in(const Point & p) const
{
	double size = tr.x-bl.x;
// 		std::cout << (p.x >= tl.x-size)  << (p.x <= br.x+size) << (p.y >= br.y-size) << (p.y <= tl.y+size) << std::endl ;
	return (p.x >= tl.x-size)  && (p.x <= br.x+size) && (p.y >= br.y-size) && (p.y <= tl.y+size);
}

bool Pixel::coOccur(const Inclusion * inc)
{
	
	return inc->in(tl) || inc->in(tr) || inc->in(br) || inc->in(bl) || in(inc->getCenter()) ;
}

void Pixel::remove(Inclusion * inc)
{
	inclusions.erase(std::find(inclusions.begin(), inclusions.end(), inc)) ;
	filled = false ;
}

bool Pixel::add(Inclusion * inc)
{
	if(filled)
		return false;

	if(!inclusions.empty())
	{
		for(size_t i = 0 ; i < inclusions.size() ; i++)
		{
			if(
			    static_cast< Circle *>(inclusions[i])->intersects(static_cast< Circle *>(inc))
			    || inc->in(inclusions[i]->getCenter())
			    || inclusions[i]->in(inc->getCenter())
			  )
				return false;
		}
		inclusions.push_back(inc) ;
		
		return true;
	}
	else
	{
		if(inc->in(tl) && inc->in(tr) && inc->in(br) && inc->in(bl))
			filled = true ;
		
		inclusions.push_back(inc) ;
		
		return true ;
	}
}

Grid::Grid(double sizeX, double sizeY, int div ) : pixels((size_t)round(div*div*std::max(sizeY/sizeX, sizeX/sizeY))), x(sizeX), y(sizeY) 
{
	if(x>y)
	{
		lengthX = div*(x/y) ;
		lengthY = div ;
	}
	else
	{
		lengthX = div ;
		lengthY = div*(y/x) ;
	}
	
	int iterator = 0 ;
	
	double psize = x/lengthX;
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			pixels[iterator] = Pixel(x*(double)(i)/(double)lengthX-x*.5+psize*.5, y*(double)(j)/(double)lengthY-y*.5+psize*.5, psize) ;
			iterator++ ;
		}
	}
}

bool Grid::add(Inclusion * inc)
{
	bool ret = true ;
	std::vector<Pixel *> cleanup ;
	for(size_t i = 0 ; i < pixels.size() ; i++)
	{
		if(pixels[i].coOccur(inc))
		{
			if(pixels[i].add(inc))
			{
				cleanup.push_back(&pixels[i]) ;
			}
			else
			{
				for(size_t j = 0 ; j < cleanup.size() ; j++)
				{
					cleanup[j]->remove(inc) ;
				}
				return false ;
			}
			
		}
		
	}
	
	return ret ;
}
