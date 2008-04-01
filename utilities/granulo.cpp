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
#include <iostream>

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

std::vector <Inclusion *> Granulo::operator()(double rayon_granulat,double pourcentMasseMin)
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
	while (pourcentMasse>pourcentMasseMin)
	{
		
		double masse_granulat =pow(( densite*rayon_granulat*rayon_granulat*M_PI), .666666666666666); // masse plus plus gros granulat
		std::cout<<"m granulat  " << masse_granulat<<std::endl;
		if (masse_granulat > masseInitiale)
		{
			std::cerr<<"La masse du granultat pese plus que la masse totale!"<<std::endl;
			return rayon ;
		}
		double volume = rayon_granulat*rayon_granulat*M_PI; // volume du plus gros granulat
		//std::cout<<"V granulat  " <<volume<<std::endl;
		
		if(rayon_granulat<=0.51)
		{
		volumeGranulatsPasPlaces +=volume;
	
		}
		masseReste = masseReste-masse_granulat; // soustraction de la masse du granulat à la masse de granultats restante
// 		std::cout<< i<<"   "<<masse_granulat<<"  " <<rayon_granulat<<"  " <<volume << "  " <<pourcentMasse<< "   " << masseReste << std::endl;
		pourcentMasse = masseReste/masseInitiale; // calcul du poucentage massique qu'il reste
		//std::cout<<"pourcentmasse  "<<pourcentMasse <<volume<<std::endl;
	 	

		
		rayon_granulat = (pow(-(log(1.-(pourcentMasse*pourcentInitial))/c),1./n))/2.; // diametre du granulat suivant
		
		if(rayon_granulat<0.07)
			return rayon;
// 		rayon.push_back(rayon_granulat);
		rayon.push_back(new Inclusion(rayon_granulat, 0, 0)) ;
// 		std::cout << "rayon granulat reste  "<<rayon[i]->getRadius() << std::endl<<std::endl;
	
		i++;
	
	}
	std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
	return  rayon;

}

GranuloBolome::GranuloBolome( double masseInitiale, double densite, TypeGranulo t)
{ //constructeur
	this->masseInitiale=masseInitiale;
	this->densite=densite;
	this->type = t ;
}


std::vector <Inclusion *> GranuloBolome::operator()(double rayonGranulatMax, double pourcentMasseMin, int inclusionNumber, double itzSize)
{	

	std::vector<Inclusion *> rayon;
	if(!inclusionNumber)
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
	double volumeGranulatsPasPlaces =0;
	while (pourcentMasse>pourcentMasseMin)
	{
		
		double masse_granulat = densite*4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat; // masse plus plus gros granulat
		
		if (masse_granulat > masseInitiale)
		{
			std::cerr<<"La masse du granulat pese plus que la masse totale!"<<std::endl;
			return rayon ;
		}
		double volume = 4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat ;//4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI/(3.); // volume du plus gros granulat
		v +=  M_PI*rayon_granulat*rayon_granulat;
		if(rayon_granulat<=0.000000051)
		{
			volumeGranulatsPasPlaces +=volume;
		}

		masseReste = masseReste-masse_granulat; 
		pourcentMasse = masseReste/masseInitiale*100.; 
		
		switch(type)
		{
		case BOLOME_A:
			rayon_granulat = 0.5*(rayonGranulatMax+4.*pourcentMasse*rayonGranulatMax
          - sqrt((rayonGranulatMax +8.*pourcentMasse*rayonGranulatMax) * (rayonGranulatMax +8.*pourcentMasse*rayonGranulatMax) 
          - 16.*pourcentMasse*pourcentMasse*rayonGranulatMax*rayonGranulatMax));
			break ;
		case BOLOME_B:
			// diametre du granulat suivant
			rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
			break ;
		case BOLOME_C:
			if(rayon_granulat<0.2)
				rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
				
			else 
				rayon_granulat = 1.05*pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
			break ;
		case BOLOME_D:
			double m;//pente
			double b;//ordonnée en x=0
			
// 			if(pourcentMasse<0.045 && pourcentMasse> 0.03)
// 			{
// 				Point A(0.1,3.);
// 				Point B(0.16,4.5);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
// 			if(pourcentMasse<0.07 && pourcentMasse> 0.045)
// 			{
// 				Point A(0.16,4.5);
// 				Point B(0.25,7.);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
// 			if(pourcentMasse<0.31 && pourcentMasse> 0.07)
// 			{
// 				Point A(0.25,7.0);
// 				Point B(4.,31.);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
// 			if(pourcentMasse<0.35 && pourcentMasse> 0.31)
// 			{
// 				Point A(4.,31.);
// 				Point B(5.,35.);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
			if(pourcentMasse> 0. && pourcentMasse<10)
			{
				Point A(.000150,0.);
				Point B(.000315,10.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = std::max(7.5e-05, (pourcentMasse-b)/m*.5);
			}
			if(pourcentMasse<20 && pourcentMasse>= 10)
			{
				Point A(.000315,10.);
				Point B(.00063,20.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (pourcentMasse-b)/m*.5;
			}
			if(pourcentMasse<45 && pourcentMasse> 20)
			{
				Point A(.00063,20.);
				Point B(.00125,45.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (pourcentMasse-b)/m*.5;
			}
			if(pourcentMasse<70 && pourcentMasse> 45)
			{
				Point A(.00125,45.);
				Point B(.0025,70.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (pourcentMasse-b)/m*.5;
			}
			if(pourcentMasse<100. && pourcentMasse> 70)
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
		
		if((int)rayon.size() >= inclusionNumber || pourcentMasse < pourcentMasseMin)
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
	
	double pourcentMasse=100;
	int i=0;

	double masseReste = masseInitiale;
// 	double pourcentInitial=1-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.
	double v = 0 ;

	rayon.push_back( new Inclusion3D(rayon_granulat+itzSize, 0, 0, 0)) ;

	double volumeGranulatsPasPlaces =0;
	while (pourcentMasse>pourcentMasseMin)
	{
		
		double masse_granulat = densite*4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat; // masse plus plus gros granulat
		
		if (masse_granulat > masseInitiale)
		{
			std::cerr<<"La masse du granulat pese plus que la masse totale!"<<std::endl;
			return rayon ;
		}
		double volume = 4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat ;//4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI/(3.); // volume du plus gros granulat
		v +=  volume;
		if(rayon_granulat<=0.000000051)
		{
			volumeGranulatsPasPlaces +=volume;
		}

		masseReste = masseReste-masse_granulat; 
		pourcentMasse = masseReste/masseInitiale*100.; 
		
		switch(type)
		{
		case BOLOME_A:
			rayon_granulat = 0.5*(rayonGranulatMax+4.*pourcentMasse*rayonGranulatMax
          - sqrt((rayonGranulatMax +8.*pourcentMasse*rayonGranulatMax) * (rayonGranulatMax +8.*pourcentMasse*rayonGranulatMax) 
          - 16.*pourcentMasse*pourcentMasse*rayonGranulatMax*rayonGranulatMax));
			break ;
		case BOLOME_B:
			// diametre du granulat suivant
			rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
			break ;
		case BOLOME_C:
			if(rayon_granulat<0.2)
				rayon_granulat = pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
				
			else 
				rayon_granulat = 1.05*pourcentMasse*pourcentMasse/20000.*rayonGranulatMax;
			break ;
		case BOLOME_D:
			double m;//pente
			double b;//ordonnée en x=0
			
// 			if(pourcentMasse<0.045 && pourcentMasse> 0.03)
// 			{
// 				Point A(0.1,3.);
// 				Point B(0.16,4.5);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
// 			if(pourcentMasse<0.07 && pourcentMasse> 0.045)
// 			{
// 				Point A(0.16,4.5);
// 				Point B(0.25,7.);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
// 			if(pourcentMasse<0.31 && pourcentMasse> 0.07)
// 			{
// 				Point A(0.25,7.0);
// 				Point B(4.,31.);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
// 			if(pourcentMasse<0.35 && pourcentMasse> 0.31)
// 			{
// 				Point A(4.,31.);
// 				Point B(5.,35.);
// 				m= (B.y-A.y)/(B.x-A.x);
// 				b=A.y-(m*A.x);
// 				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
// 			}
			if(pourcentMasse> 0. && pourcentMasse<10)
			{
				Point A(.000150,0.);
				Point B(.000315,10.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = std::max(7.5e-05, (pourcentMasse-b)/m*.5);
			}
			if(pourcentMasse<20 && pourcentMasse>= 10)
			{
				Point A(.000315,10.);
				Point B(.00063,20.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (pourcentMasse-b)/m*.5;
			}
			if(pourcentMasse<45 && pourcentMasse> 20)
			{
				Point A(.00063,20.);
				Point B(.00125,45.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (pourcentMasse-b)/m*.5;
			}
			if(pourcentMasse<70 && pourcentMasse> 45)
			{
				Point A(.00125,45.);
				Point B(.0025,70.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (pourcentMasse-b)/m*.5;
			}
			if(pourcentMasse<100. && pourcentMasse> 70)
			{
				Point A(.0025,70.);
				Point B(.005,100.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (pourcentMasse-b)/m*.5;
			}

			break;
		default:
			break ;
		}
		
		if((int)rayon.size() >= inclusionNumber || pourcentMasse < pourcentMasseMin)
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
