
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
// namespace Mu
// {
// 
// 
// std::vector<Inclusion *> placement(double longueurX, double longueurY, std::vector<Inclusion *> inclusions, int *nombreGranulatsPlaces, int triesMax) ;
// // 	double masseInitiale;
// // 	double densite;
// 
// 
// 	double chiffreAleatoire(double );
// 	bool bord(double , double , double , double , double );
// 
// // 	std::vector <Inclusion *> place(double ,double, std::vector<Inclusion *>);//Inclusion *
// // 
// // 	
// 
// 
// } ;

std::vector<Mu::Inclusion *> placement(double longueurX, double longueurY, std::vector<Mu::Inclusion *> inclusions, int *nombreGranulatsPlaces, int triesMax) ;
// 	double masseInitiale;
// 	double densite;


	double chiffreAleatoire(double );
	bool bord(double , double , double , double , double );

#endif
