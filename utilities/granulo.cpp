//
// C++ Implementation: granulo
//
// Description:
//
//
// Author:  <>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <cmath>
#include <vector>
#include "granulo.h"
#include <iostream>  // I/O 
#include <fstream>   // file I/O
#include <iomanip>   // format manipulation
#include "placement.h"
#include "../geometry/geometry_base.h"

using namespace Mu ;



Granulo::Granulo(double a, double b,double masseInitiale, double densite)//,double masseInitiale, double densite
{
    c = a ; //constructeur
    n = b ;
    this->masseInitiale=masseInitiale;
    this->densite=densite;


}

Granulo::Granulo()//,double masseInitiale, double densite
{
    c = 0.1 ; //constructeur
    n = 1. ;
    this->masseInitiale=1000;
    this->densite=2.1;


}

std::vector <Inclusion *> Granulo::operator()(double rayon_granulat, double pourcentMasseMin, int inclusionNumber , double itzSize)
{
    std::vector<Inclusion *> rayon;


    double pourcentMasse=1;
    int i=0;

    double masseReste = masseInitiale;
    double pourcentInitial=1.-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.

// 	rayon.push_back(rayon_granulat);
    rayon.push_back( new Inclusion(rayon_granulat, 0, 0)) ;
    std::cout<< "n°"<<"masse granulat" << "  " <<"rayon granulat"<<"  " <<"volume granulat" << "  "<<"pourcent masse" <<"   "<<"masse restante"<<std::endl;
    double volumeGranulatsPasPlaces =0;
    while((pourcentMasse>pourcentMasseMin) && (i < inclusionNumber))
    {

        double masse_granulat =pow(( densite*rayon_granulat*rayon_granulat*M_PI), .666666666666666); // masse plus plus gros granulat
//        std::cout<<"m granulat  " << masse_granulat<<std::endl;
        if (masse_granulat > masseInitiale)
        {
            std::cerr<<"La masse du granultat pese plus que la masse totale!"<<std::endl;
            return rayon ;
        }
        double volume = rayon_granulat*rayon_granulat*M_PI; // volume du plus gros granulat
        //std::cout<<"V granulat  " <<volume<<std::endl;

        if (rayon_granulat<=0.51)
        {
            volumeGranulatsPasPlaces +=volume;

        }
        masseReste = masseReste-masse_granulat; // soustraction de la masse du granulat à la masse de granultats restante
// 		std::cout<< i<<"   "<<masse_granulat<<"  " <<rayon_granulat<<"  " <<volume << "  " <<pourcentMasse<< "   " << masseReste << std::endl;
        pourcentMasse = masseReste/masseInitiale; // calcul du poucentage massique qu'il reste
        //std::cout<<"pourcentmasse  "<<pourcentMasse <<volume<<std::endl;



        rayon_granulat = (pow(-(log(1.-(pourcentMasse*pourcentInitial))/c),1./n))/2.; // diametre du granulat suivant

	if(i%100 == 0)
		std::cout << rayon_granulat << std::endl ;

//        if (rayon_granulat<0.07)
  //          return rayon;
// 		rayon.push_back(rayon_granulat);
        rayon.push_back(new Inclusion(rayon_granulat, 0, 0)) ;
// 		std::cout << "rayon granulat reste  "<<rayon[i]->getRadius() << std::endl<<std::endl;

        i++;

    }
    std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;

    return  rayon;

}

std::vector<EllipsoidalInclusion *> Granulo::operator()(bool ell, Point * biais, double rayon_granulat, double pourcentMasseMin, double rfactor, int inclusionNumber , double itzSize)
{
    std::vector<EllipsoidalInclusion *> rayon;


    double pourcentMasse=1;
    int i=0;

    double sfactor = rfactor ;
    double masseReste = masseInitiale;
    double pourcentInitial=1.-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.

// 	rayon.push_back(rayon_granulat);
    rayon.push_back( new EllipsoidalInclusion(rayon_granulat/sfactor, rayon_granulat,0., 0.,biais->x,biais->y)) ;
	std::cout << rayon_granulat << std::endl ;
//    std::cout<< "n°"<<"masse granulat" << "  " <<"rayon granulat"<<"  " <<"volume granulat" << "  "<<"pourcent masse" <<"   "<<"masse restante"<<std::endl;
    double volumeGranulatsPasPlaces =0;
    double alea = 0.1 ;
    Point b_(biais->x, biais->y) ;
    while (pourcentMasse>pourcentMasseMin && i < inclusionNumber)
    {
	if(ell)
	{
		alea = (double)rand()/(double)RAND_MAX ;
		sfactor = (rfactor + (1 - rfactor) * alea) ; // * alea ;
		alea = (double)rand()/(double)RAND_MAX ;
		b_.setX(biais->x + 0.3 * alea - 0.15) ;
		alea = (double)rand()/(double)RAND_MAX ;
		b_.setY(biais->y + 0.3 * alea - 0.15) ;		
		//std::cout << b_.x << " " << b_.y << std::endl ;
		//std::cout << sfactor << std::endl ;
	}

        double masse_granulat =pow(( densite*rayon[i]->area()), .666666666666666); // masse plus plus gros granulat
//        std::cout<<"m granulat  " << masse_granulat<<std::endl;
        if (masse_granulat > masseInitiale)
        {
            std::cerr<<"La masse du granultat pese plus que la masse totale!"<<std::endl;
            return rayon ;
        }
        double volume = rayon[i]->area() ; // volume du plus gros granulat
        //std::cout<<"V granulat  " <<volume<<std::endl;

        if (rayon_granulat<=0.51)
        {
            volumeGranulatsPasPlaces +=volume;

        }
        masseReste = masseReste-masse_granulat; // soustraction de la masse du granulat à la masse de granultats restante
// 		std::cout<< i<<"   "<<masse_granulat<<"  " <<rayon_granulat<<"  " <<volume << "  " <<pourcentMasse<< "   " << masseReste << std::endl;
        pourcentMasse = masseReste/masseInitiale; // calcul du poucentage massique qu'il reste
        //std::cout<<"pourcentmasse  "<<pourcentMasse <<volume<<std::endl;



        rayon_granulat = (pow(-(log(1.-(pourcentMasse*pourcentInitial))/c),1./n))/2.; // diametre du granulat suivant

	if(i%100 == 0)
		std::cout << rayon_granulat << std::endl ;

//        if (rayon_granulat<0.07)
  //          return rayon;
// 		rayon.push_back(rayon_granulat);
        rayon.push_back(new EllipsoidalInclusion(rayon_granulat/sfactor,rayon_granulat, 0., 0.,b_.x,b_.y)) ;
// 		std::cout << "rayon granulat reste  "<<rayon[i]->getRadius() << std::endl<<std::endl;

        i++;

    }
    std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
    return  rayon;

}


GranuloBolome::GranuloBolome( double masseInitiale_, double densite_, TypeGranulo t)
{ //constructeur
    this->masseInitiale=masseInitiale_;
    this->densite=densite_;
    this->type = t ;
}


std::vector <Inclusion *> GranuloBolome::operator()(double rayonGranulatMax, double pourcentMasseMin, int inclusionNumber, double itzSize)
{

//     std::cout << "generating 3D distribution" << std::endl ;
// 
// 		std::vector <Inclusion3D *> threeDInclusions = (*this)(true,rayonGranulatMax, pourcentMasseMin, pow(inclusionNumber, 3./2.), 0) ;
// 		std::vector<Feature *> feats ;
// 		feats.insert(feats.end(), threeDInclusions.begin(), threeDInclusions.end()) ;
// 		Hexahedron hex(.04, .04, .04, 0,0,0) ;
// 		int tot ;
// 		feats= placement(&hex, feats, &tot, 6400);

    std::vector<Inclusion *> rayon;

//     for (size_t i =0 ; i < feats.size() ;i++)
//     {
// 			if(std::abs(feats[i]->getCenter().z) < feats[i]->getRadius())
// 			{
// 				rayon.push_back(new Inclusion(sqrt(feats[i]->getRadius()*feats[i]->getRadius()-feats[i]->getCenter().z*feats[i]->getCenter().z), feats[i]->getCenter().x, feats[i]->getCenter().y)) ;
// 			}
// 			delete feats[i] ;
//     }
// 
//     return rayon ;

    if (!inclusionNumber)
        return rayon ;
    double rayon_granulat=rayonGranulatMax;

    double pourcentMasse=100;
    int i=0;

    double masseReste = masseInitiale;
// 	double pourcentInitial=1-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.
    double v = 0 ;

// 	rayon.push_back(rayon_granulat);
    rayon.push_back( new Inclusion(rayon_granulat+itzSize, 0, 0)) ;
// 	std::cout<< "n°"<<"masse granulat" << "  " <<"rayon granulat"<<"  " <<"volume granulat" << "  "<<"pourcent masse" <<"   "<<"masse restante"<<std::endl;
// 	size_t h=0;
    double volumeGranulatsPasPlaces = 0;
    while (pourcentMasse>pourcentMasseMin)
    {

        double masse_granulat = densite*4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat; // masse plus plus gros granulat

        if (masse_granulat > masseInitiale)
        {
            std::cerr << "La masse du granulat " << masse_granulat << " pese plus que la masse totale " <<  masseInitiale << "!"  << std::endl;
            return rayon ;
        }
        double volume = 4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat ;//4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI/(3.); // volume du plus gros granulat
        v +=  M_PI*rayon_granulat*rayon_granulat;
        if (rayon_granulat<=0.000000051)
        {
            volumeGranulatsPasPlaces +=volume;
        }

        masseReste = masseReste-masse_granulat;
        pourcentMasse = masseReste/masseInitiale*100.;

        switch (type)
        {
        case BOLOME_A:
            rayon_granulat = 0.5*(rayonGranulatMax+4.*pourcentMasse*rayonGranulatMax
                                  - sqrt((rayonGranulatMax +8.*pourcentMasse*rayonGranulatMax)
                                         * (rayonGranulatMax +8.*pourcentMasse*rayonGranulatMax)
                                         - 16.*pourcentMasse*pourcentMasse*rayonGranulatMax*rayonGranulatMax));
            rayon_granulat = (rayon_granulat*.5*sqrt(rayon_granulat*rayon_granulat*.75) + rayon_granulat*rayon_granulat*atan(rayon_granulat*.5/sqrt(rayon_granulat*rayon_granulat*.75 ))*.5 )/rayon_granulat/100.;
// 			std::cout << rayon_granulat << std::endl ;
            break ;
        case BOLOME_B:
            // diametre du granulat suivant
            rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
            break ;
        case BOLOME_C:
            if (rayon_granulat<0.2)
                rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;

            else
                rayon_granulat = 1.05*pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
            break ;
        case BOLOME_D:
            double m;//pente
            double b;//ordonnée en x=0

            if (pourcentMasse> 0. && pourcentMasse<10)
            {
                Point A(.000150,0.);
                Point B(.000315,10.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = std::max(7.5e-05, (pourcentMasse-b)/m*.5);
            }
            if (pourcentMasse<20 && pourcentMasse>= 10)
            {
                Point A(.000315,10.);
                Point B(.00063,20.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }
            if (pourcentMasse<45 && pourcentMasse> 20)
            {
                Point A(.00063,20.);
                Point B(.00125,45.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }
            if (pourcentMasse<70 && pourcentMasse> 45)
            {
                Point A(.00125,45.);
                Point B(.0025,70.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }
            if (pourcentMasse<100. && pourcentMasse> 70)
            {
                Point A(.0025,70.);
                Point B(.005,100.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }

            rayon_granulat = (rayon_granulat*.5*sqrt(rayon_granulat*rayon_granulat*.75) + rayon_granulat*rayon_granulat*atan(rayon_granulat*.5/sqrt(rayon_granulat*rayon_granulat*.75 ))*.5 )/rayon_granulat;
            break;
        default:
            break ;
        }

        if ((int)rayon.size() >= inclusionNumber || pourcentMasse < pourcentMasseMin)
        {
            std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
            std::cout<<"volumeAgg "<<v<<std::endl;
            std::cout << rayon.size() <<" particles" << std::endl ;
            return rayon;

        }

        rayon.push_back(new Inclusion(rayon_granulat+itzSize, 0, 0)) ;
        i++;


    }
    std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
    std::cout<<"volumeAgg "<<v<<std::endl;
    std::cout << rayon.size() << " particles" << std::endl ;
    return  rayon;

}

std::vector <Inclusion3D *> GranuloBolome::operator()(bool,double rayonGranulatMax, double pourcentMasseMin, int inclusionNumber, double itzSize)
{
    std::vector<Inclusion3D *> rayon;
    double rayon_granulat=rayonGranulatMax;

    double pourcentMasse=1.;
    int i=0;
    if (!inclusionNumber)
        return rayon ;

		
    double masseReste = masseInitiale;
// 	double pourcentInitial=1-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.
    double v = 0 ;

    rayon.push_back( new Inclusion3D(rayon_granulat+itzSize, 0, 0, 0)) ;

    double volumeGranulatsPasPlaces =0;
    while (pourcentMasse*100.>pourcentMasseMin)
    {

        double masse_granulat = densite*4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat; // masse plus plus gros granulat
        if (masse_granulat > masseInitiale)
        {
            std::cerr<<"La masse du granulat pese plus que la masse totale!"<<std::endl;
            return rayon ;
        }
        double volume = 4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat ;//4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI/(3.); // volume du plus gros granulat
        v +=  volume;
        if (rayon_granulat<=0.000000051)
        {
            volumeGranulatsPasPlaces +=volume;
        }

        masseReste -= masse_granulat;
        pourcentMasse = masseReste/masseInitiale;

        switch (type)
        {
        case BOLOME_A:
        {
            double A_B = 8. ;
            double val = (pourcentMasse*100. - A_B) / (100. - A_B) ;
            rayon_granulat = rayonGranulatMax*val*val;
            break ;
        }
        case BOLOME_B:
            // diametre du granulat suivant
            rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
            break ;
        case BOLOME_C:
            if (rayon_granulat<0.2)
                rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;

            else
                rayon_granulat = 1.05*pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
            break ;
        case BOLOME_D:
            double m;//pente
            double b;//ordonnée en x=0


            if (pourcentMasse> 0. && pourcentMasse<=.10)
            {
                Point A(.000150,0.);
                Point B(.000315,.10);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = std::max(7.5e-05, (pourcentMasse-b)/m*.5);
            }
            if (pourcentMasse<=.20 && pourcentMasse> .10)
            {
                Point A(.000315,.10);
                Point B(.00063,.20);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }
            if (pourcentMasse<=.45 && pourcentMasse> .20)
            {
                Point A(.00063,.20);
                Point B(.00125,.45);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }
            if (pourcentMasse<=.70 && pourcentMasse> .45)
            {
                Point A(.00125,.45);
                Point B(.0025,.70);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }
            if (pourcentMasse<=1. && pourcentMasse> .70)
            {
                Point A(.0025,.70);
                Point B(.005,1.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                rayon_granulat = (pourcentMasse-b)/m*.5;
            }

            break;
        default:
            break ;
        }

        if ((int)rayon.size() >= inclusionNumber || pourcentMasse < pourcentMasseMin/100.)
        {
            std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
            std::cout<<"volumeAgg "<<v<<std::endl;
            std::cout << rayon.size() <<" particles" << std::endl ;
            return rayon;

        }

        rayon.push_back(new Inclusion3D(rayon_granulat+itzSize, 0, 0, 0)) ;
        i++;


    }
    std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
    std::cout<<"volumeAgg "<<v<<std::endl;
    std::cout << rayon.size() << " particles" << std::endl ;
    return  rayon;

}

GranuloFromFile::GranuloFromFile(std::string fname, std::vector<double> coef, double mm, double kg)
{
	filename = fname ;
	double tot = 0 ;
	for(int i=0 ; i < coef.size() ; i++)
		tot += coef[i] ;
	if(tot == 0)
		tot = 1 ;
	std::vector<double> coef_good ;
	for(int i=0 ; i < coef.size() ; i++)
		coef_good.push_back(coef[i]/tot) ;

	char tsz [256] ;
	char tm [256] ;
	double temp_mass = 0 ;

	std::ifstream granuloreader  ;
	granuloreader.open(filename.c_str(), std::ios::in) ;	
	while(!granuloreader.eof())
	{
		temp_mass = 0 ;
		granuloreader >> tsz ;
		size.push_back(atof(tsz)*mm) ;
		for(int i=0 ; i < coef_good.size() ; i++)
		{
			granuloreader >> tm ;
			temp_mass += coef_good[i]*atof(tm)*kg ;
		}
		mass.push_back(temp_mass) ;
	}

	granuloreader.close() ;
	size.pop_back() ;
	mass.pop_back() ;

	while(mass[1]<0.0001)	
	{
		size.erase(size.begin()) ;
		mass.erase(mass.begin()) ;
	}
	int i = mass.size() ;
	while(mass[i-1]<0.0001)	
	{
		size.pop_back() ;
		mass.pop_back() ;
		i = mass.size() ;
	}
	mass.erase(mass.begin()) ;
	std::cout << size.size() << " ; " << mass.size() << std::endl ;
}

std::vector<Inclusion *> GranuloFromFile::getCircleInclusion(double density, int n_agg, double area)
{
	double mass_goal = density * pow(area,1.5) ;
	double mass_tot = 0 ;
	for(int i = 0 ; i < mass.size() ; i++)
		mass_tot += mass[i] ;
	for(int i = 0 ; i < mass.size() ; i++)
		mass[i] = mass[i] * mass_goal / mass_tot ;	
	double first_mass = (4/3)*M_PI*size[0]*size[0]*size[0] ;
	for(int i = 1 ; i < mass.size() ; i++)
		mass[i] = mass[i] * first_mass / mass[0] ;
	mass[0] = first_mass ;
	double radius = 0 ;
	double placedMass = 0 ;
	double placedMassIteration = 0 ;
//	std::vector<std::vector<Inclusion *> > inc_full ;
	std::vector<double> i_full ;
	std::vector<Inclusion *> inc ;
	int i_agg = 0 ;
	for(int i = 0 ; i < mass.size() ; i++)
	{
		int i_iter = 0 ;
//		std::vector<Inclusion *> inc_iter ;
		placedMassIteration = 0 ;
		while(placedMassIteration < mass[i])
		{
			radius = size[i+1] + ((double)rand()/(double)RAND_MAX) * (size[i] - size[i+1]) ;
//			inc_iter.push_back(new Inclusion(radius,0,0)) ;
			placedMassIteration += (4/3)*M_PI*radius*radius*radius*density ;
			i_iter++ ;
			i_agg++ ;
		}
		placedMass += placedMassIteration ;
//		inc_full.push_back(inc_iter) ;
		i_full.push_back(i_iter) ;
		std::cout << i << " ; " << i_iter << " ; " << placedMassIteration << "/" << mass[i] << std::endl ;
	}
	double n_norm = 0 ;
	for(int i = 0 ; i < i_full.size() ; i++)
		n_norm += pow(i_full[i], 0.333333) ;
	double factor = (double) n_agg / (double) n_norm ;
	for(int i = 0 ; i < mass.size() ; i++)
	{
		int n_iter = pow(i_full[i], 0.333333) * factor ;
		if(n_iter < 1)
			n_iter = 1 ;
		for(int j = 0 ; j < n_iter ; j++)
		{
			radius = size[i+1] + ((double)rand()/(double)RAND_MAX) * (size[i] - size[i+1]) ;
			inc.push_back(new Inclusion(radius,0,0)) ;
		}		
		std::cout << n_iter << std::endl ;
	}
	while(inc.size() < n_agg)
	{
		std::cout << "here" << std::endl ;
		radius = size[size.size()] + ((double)rand()/(double)RAND_MAX) * (size[size.size() - 1] - size[size.size()]) ;
		inc.push_back(new Inclusion(radius,0,0)) ;
	}
	double area_tot = 0 ;
	for(int i = 0 ; i < inc.size() ; i++)
		area_tot += M_PI*(inc[i]->getRadius())*(inc[i]->getRadius()) ;
	std::cout << "area for " << n_agg << " inclusions => " << area_tot << " (goal = " << area << " )" << std::endl ;
	while(area_tot > area * 1.01)
	{
		area_tot = 0 ;
		int n_rand = (inc.size() - 1) * (double)rand() / (double) RAND_MAX ;
		inc.erase(inc.begin()+n_rand) ;		
		for(int i = 0 ; i < inc.size() ; i++)
			area_tot += M_PI*(inc[i]->getRadius())*(inc[i]->getRadius()) ;
	}
	std::cout << inc.size() << " inclusions placed over " << n_agg << std::endl ;
	return inc ;
}

void GranuloFromFile::resize(double newSize)
{
	double logNewSize = log10(newSize) ;
	double logMaxSize = log10(size[0]) ;
	double logThisSize = 0 ;
	size[0] = newSize ;
//	double thisRadius = newSize ;
//	double thisUnitMass = 0 ;
//	int thisUnit = 0 ;
	for(int i = 0 ; i < mass.size() ; i++)
	{
		logThisSize = log10(size[i+1]) ;
//		thisRadius = 0.5 * (size[i] + size[i+1]) ;
//		thisUnitMass = (4/3) * M_PI * thisRadius * thisRadius * thisRadius ;
//		thisUnit = mass[i] / thisUnitMass ;
		size[i+1] = pow(10,logNewSize - logMaxSize + logThisSize) ;
//		thisRadius = 0.5 * (size[i] + size[i+1]) ;
//		thisUnitMass = (4/3) * M_PI * thisRadius * thisRadius * thisRadius ;
//		mass[i] = thisUnit * thisUnitMass ;
	}

}

