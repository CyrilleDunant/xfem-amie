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
		
		double masse_granulat =( densite*4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI)/(3.); // masse plus plus gros granulat
		std::cout<<"m granulat  " << masse_granulat<<std::endl;
		if (masse_granulat > masseInitiale)
		{
			std::cerr<<"La masse du granultat pese plus que la masse totale!"<<std::endl;
			return rayon ;
		}
		double volume = 4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI/(3.); // volume du plus gros granulat
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


std::vector <Inclusion *> GranuloBolome::operator()(double rayonGranulatMax, double pourcentMasseMin)
{	
	std::vector<Inclusion *> rayon;
	double rayon_granulat=rayonGranulatMax;
	
	double pourcentMasse=1;
	int i=0;

	double masseReste = masseInitiale;
// 	double pourcentInitial=1-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.
	

// 	rayon.push_back(rayon_granulat);
	rayon.push_back( new Inclusion(rayon_granulat, 0, 0)) ;
// 	std::cout<< "n°"<<"masse granulat" << "  " <<"rayon granulat"<<"  " <<"volume granulat" << "  "<<"pourcent masse" <<"   "<<"masse restante"<<std::endl;
// 	size_t h=0;
	double volumeGranulatsPasPlaces =0;
	while (pourcentMasse>pourcentMasseMin)
	{
		
		double masse_granulat =( densite*4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI)/(3.); // masse plus plus gros granulat
// 		std::cout<<"m granulat  " << masse_granulat<<std::endl;
		if (masse_granulat > masseInitiale)
		{
			std::cerr<<"La masse du granultat pese plus que la masse totale!"<<std::endl;
			return rayon ;
		}
		double volume = 4.*rayon_granulat*rayon_granulat*rayon_granulat*M_PI/(3.); // volume du plus gros granulat
		//std::cout<<"V granulat  " <<volume<<std::endl;

		if(rayon_granulat<=0.0051)
		{
		volumeGranulatsPasPlaces +=volume;
	
	}


		masseReste = masseReste-masse_granulat; // soustraction de la masse du granulat à la masse de granultats restante
// 		std::cout<< i<<"   "<<masse_granulat<<"  " <<rayon_granulat<<"  " <<volume << "  " <<pourcentMasse<< "   " << masseReste << std::endl;
		pourcentMasse = masseReste/masseInitiale; // calcul du poucentage massique qu'il reste
// 		std::cout<<"pourcentmasse  "<<pourcentMasse <<volume<<std::endl;
	 	

		switch(type)
		{
		case BOLOME_A:
			rayon_granulat = 0.5*(2.*rayonGranulatMax+(2.*100.*pourcentMasse*2.*rayonGranulatMax/50.)-pow(((2.*rayonGranulatMax +(2.*100.*pourcentMasse*2.*rayonGranulatMax/50.))*(2.*rayonGranulatMax +(2.*100.*pourcentMasse*2.*rayonGranulatMax/50.))) -(4.*10000.*pourcentMasse*pourcentMasse*4.*rayonGranulatMax*rayonGranulatMax/2500.),0.5))/2.;
			break ;
		case BOLOME_B:
			rayon_granulat = (pourcentMasse*pourcentMasse*10000.*2.*rayonGranulatMax)/20000.;   // diametre du granulat suivant
// 			std::cout<<rayon_granulat<<"    " <<pourcentMasse<<std::endl;
			break ;
		case BOLOME_C:
			if(rayon_granulat<0.2)
				rayon_granulat = (pourcentMasse*pourcentMasse*10000.*2.*rayonGranulatMax)/20000.;
				
			else 
				rayon_granulat = 1.05*(pourcentMasse*pourcentMasse*10000.*2.*rayonGranulatMax)/20000.;
			break ;
		case BOLOME_D:
			double m;//pente
			double b;//ordonnée en x=0
			
			if(pourcentMasse<0.045 && pourcentMasse> 0.03)
			{
				Point A(0.1,3.);
				Point B(0.16,4.5);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.07 && pourcentMasse> 0.045)
			{
				Point A(0.16,4.5);
				Point B(0.25,7.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.31 && pourcentMasse> 0.07)
			{
				Point A(0.25,7.0);
				Point B(4.,31.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.35 && pourcentMasse> 0.31)
			{
				Point A(4.,31.);
				Point B(5.,35.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.42 && pourcentMasse> 0.35)
			{
				Point A(5.,35.);
				Point B(6.3,42.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.51 && pourcentMasse> 0.42)
			{
				Point A(6.3,42.);
				Point B(8.,51.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.57 && pourcentMasse> 0.51)
			{
				Point A(8.,51.);
				Point B(10.,57.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.82 && pourcentMasse> 0.57)
			{
				Point A(10.,57.);
				Point B(16.,82.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<0.97 && pourcentMasse> 0.82)
			{
				Point A(16.,82.);
				Point B(25.,97.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			if(pourcentMasse<1. && pourcentMasse> 0.97)
			{
				Point A(25.,97.);
				Point B(32.,100.);
				m= (B.y-A.y)/(B.x-A.x);
				b=A.y-(m*A.x);
				rayon_granulat = (100.*pourcentMasse-b)/m/2.;
			}
			break;
		default:
			break ;
		}
		
		if(rayon_granulat<0.00008)
		{
			std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
			return rayon;
			
		}
// 	
		rayon.push_back(new Inclusion(rayon_granulat, 0, 0)) ;
		
// 		double chiffreAleatoire = (double) rand()/RAND_MAX*532;
// 		if(chiffreAleatoire<18*rayon_granulat)//0.1:: 50x50 :32.2; 100x100: 19.5;  200x200: 18
// 		{
// 		rayon.push_back(new Inclusion(rayon_granulat, 0, 0)) ;	
// 		}

// 		std::cout << "rayon granulat reste  "<<rayon[i]->getRadius() << std::endl;
		
		i++;

	
	}
	std::cout<<"volumeGranulatsPasPlaces "<<volumeGranulatsPasPlaces<<std::endl;
	return  rayon;

}