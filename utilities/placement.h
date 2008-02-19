
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __PLACEMENT_H__
#define __PLACEMENT_H__

#include <iostream>
#include <vector>

#include "../features/inclusion.h"
namespace Mu
{

class Pixel
{
protected:
	std::vector<Mu::Inclusion *> inclusions ;
	Point tl ;
	Point tr ;
	Point bl ;
	Point br ;
	bool filled ;
public:
	Pixel();
	
	Pixel(double x, double y, double s) ;

	const std::vector<Mu::Inclusion *> & getInclusions() const;
	
	std::vector<Mu::Inclusion *> & getInclusions();
	
	bool in(const Point & p) const;

	bool coOccur(const Inclusion * inc);
	
	void remove(Inclusion * inc);
	
	bool add(Inclusion * inc);

} ;

class Grid
{
protected:
	std::valarray<Pixel> pixels;
	double x ;
	double y ;
	size_t lengthX ;
	size_t lengthY ;
public:
		
	Grid(double sizeX, double sizeY, int div );
	
	bool add(Inclusion * inc);
} ;

	std::vector<Mu::Inclusion *> placement(double longueurX, double longueurY, std::vector<Mu::Inclusion *> inclusions, int *nombreGranulatsPlaces, int triesMax) ;
// 	double masseInitiale;
// 	double densite;


	double chiffreAleatoire(double );
	bool bord(double , double , double , double , double );
} ;

#endif
