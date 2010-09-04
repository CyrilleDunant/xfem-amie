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
#include <cstring>
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

std::vector <Inclusion *> Granulo::operator()(double rayon_granulat, double pourcentMasseMin, double inclusionNumber , double itzSize, bool verbose)
{
    std::vector<Inclusion *> rayon;


    double pourcentMasse=1;
    int i=0;

    double remainderMass = masseInitiale;
    double pourcentInitial=1.-exp(-c*pow(rayon_granulat,n));

    rayon.push_back( new Inclusion(rayon_granulat, 0, 0)) ;
    double volumeGranulatsPasPlaces =0;
    while((pourcentMasse>pourcentMasseMin) && (i < inclusionNumber))
    {

        double aggregateMass =pow(( densite*rayon_granulat*rayon_granulat*M_PI), .666666666666666); 
        if (aggregateMass > masseInitiale)
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
        remainderMass = remainderMass-aggregateMass; 
        pourcentMasse = remainderMass/masseInitiale; 



        rayon_granulat = (pow(-(log(1.-(pourcentMasse*pourcentInitial))/c),1./n))/2.; 

	if(i%100 == 0)
		std::cout << rayon_granulat << std::endl ;


        rayon.push_back(new Inclusion(rayon_granulat, 0, 0)) ;

        i++;

    }
    if(verbose)
		std::cout<<"Unused Aggregate Volume "<<volumeGranulatsPasPlaces<<std::endl;

    return  rayon;

}

std::vector<EllipsoidalInclusion *> Granulo::operator()(bool ell, Point * biais, double rayon_granulat, double pourcentMasseMin, double rfactor, int inclusionNumber , double itzSize, bool verbose)
{
    std::vector<EllipsoidalInclusion *> rayon;


    double pourcentMasse=1;
    int i=0;

    double sfactor = rfactor ;
    double remainderMass = masseInitiale;
    double pourcentInitial=1.-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.

// 	rayon.push_back(rayon_granulat);
    rayon.push_back( new EllipsoidalInclusion(Point(0.,0.),Point(biais->x,biais->y)*rayon_granulat,sfactor)) ;
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

        double aggregateMass =pow(( densite*rayon[i]->area()), .666666666666666); // masse plus plus gros granulat
//        std::cout<<"m granulat  " << aggregateMass<<std::endl;
        if (aggregateMass > masseInitiale)
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
        remainderMass = remainderMass-aggregateMass; // soustraction de la masse du granulat à la masse de granultats restante
// 		std::cout<< i<<"   "<<aggregateMass<<"  " <<rayon_granulat<<"  " <<volume << "  " <<pourcentMasse<< "   " << remainderMass << std::endl;
        pourcentMasse = remainderMass/masseInitiale; // calcul du poucentage massique qu'il reste
        //std::cout<<"pourcentmasse  "<<pourcentMasse <<volume<<std::endl;



        rayon_granulat = (pow(-(log(1.-(pourcentMasse*pourcentInitial))/c),1./n))/2.; // diametre du granulat suivant

	if(i%100 == 0)
		std::cout << rayon_granulat << std::endl ;

//        if (rayon_granulat<0.07)
  //          return rayon;
// 		rayon.push_back(rayon_granulat);
        rayon.push_back(new EllipsoidalInclusion(Point(0.,0.),Point(b_.x,b_.y),sfactor)) ;
// 		std::cout << "rayon granulat reste  "<<rayon[i]->getRadius() << std::endl<<std::endl;

        i++;

    }
    if(verbose)
		std::cout<<"Unused Aggregate Volume "<<volumeGranulatsPasPlaces<<std::endl;
    return  rayon;

}


GranuloBolome::GranuloBolome( double masseInitiale_, double densite_, TypeGranulo t)
{ //constructeur
    this->masseInitiale=masseInitiale_;
    this->densite=densite_;
    this->type = t ;
}


std::vector <Inclusion *> GranuloBolome::operator()(double rmax, double pourcentMasseMin, double inclusionNumber, double itzSize, bool verbose)
{

    std::vector<double> radius;
	double minDiameter = 0.00000120 ;
    if (!inclusionNumber)
        return std::vector <Inclusion *>() ;
	
    double diameter=rmax*2.;
	double dmax = rmax*2. ;

    double pourcentMasse=100;
    double remainderMass = masseInitiale;
	
    radius.push_back( diameter*.5+itzSize) ;

	while (pourcentMasse > pourcentMasseMin)
    {
		if(radius.size() >= inclusionNumber)
			goto out ;
// 		std::cout << pourcentMasse << "  " <<std::flush ;
// 		std::cout << diameter << std::endl ;
        double aggregateMass = /*densite*M_PI*diameter*diameter*.25/(ratio*ratio) ;*/densite*(1./6.)*M_PI*diameter*diameter*diameter; 

        if (aggregateMass > masseInitiale && diameter > minDiameter)
        {
            std::cerr << "Aggregate mass " << aggregateMass << " is larger than total mass " <<  masseInitiale << "!"  << std::endl;
            return std::vector <Inclusion *>() ;
        }

        remainderMass = remainderMass-aggregateMass;
		if(remainderMass < 0)
			goto out ;
        pourcentMasse = remainderMass/masseInitiale*100.;
        switch (type)
        {
        case BOLOME_A:
		{

			double b = 1.+pourcentMasse/25. ;

			diameter = dmax*(b - sqrt(b*b-pourcentMasse*pourcentMasse/(25.*25.)))*.5;
			
			if(diameter <= 100.*POINT_TOLERANCE)
				 goto out ;
            break ;
		}
        case BOLOME_B:
		{
            // diametre du granulat suivant
            diameter = pourcentMasse*pourcentMasse/20000.*dmax;
            break ;
		}
        case BOLOME_C:
		{
            if (diameter<0.2)
                diameter = pourcentMasse*pourcentMasse/20000.*dmax;

            else
                diameter = 1.05*pourcentMasse*pourcentMasse/20000.*dmax;
            break ;
		}
        case BOLOME_D:
		{
            double m;//pente
            double b;//ordonnée en x=0

            if (pourcentMasse> 0. && pourcentMasse<10)
            {
                Point A(.000150,0.);
                Point B(.000315,10.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                diameter = std::max(7.5e-05, (pourcentMasse-b)/m);
            }
            if (pourcentMasse<20 && pourcentMasse>= 10)
            {
                Point A(.000315,10.);
                Point B(.00063,20.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                diameter = (pourcentMasse-b)/m*.5;
            }
            if (pourcentMasse<45 && pourcentMasse> 20)
            {
                Point A(.00063,20.);
                Point B(.00125,45.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                diameter = (pourcentMasse-b)/m;
            }
            if (pourcentMasse<70 && pourcentMasse> 45)
            {
                Point A(.00125,45.);
                Point B(.0025,70.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                diameter = (pourcentMasse-b)/m;
            }
            if (pourcentMasse<100. && pourcentMasse> 70)
            {
                Point A(.0025,70.);
                Point B(.005,100.);
                m= (B.y-A.y)/(B.x-A.x);
                b=A.y-(m*A.x);
                diameter = (pourcentMasse-b)/m;
            }

            break;
		}
        default:
            break ;
        }
        
        radius.push_back(diameter*.5) ;
		if(radius.size() >= inclusionNumber)
			break ;
    }
    
 out:   
    std::vector<Inclusion *> particles ;
	if(radius.empty())
		return  particles;
		
	std::sort(radius.begin(), radius.end()) ;
	std::reverse(radius.begin(), radius.end());
	
	double v = 0 ;
	for(size_t j = 0 ; j < radius.size() ; j++)
	{
		particles.push_back(new Inclusion(radius[j], 0, 0));
		v += particles.back()->area() ;
	}

	if(verbose)
	{
		std::cout<<"Total Aggregate Volume "<< v << std::endl;
		std::cout << "n = " << particles.size() << ", largest r = " << particles.front()->getRadius() 
				<< ", smallest r =" << particles.back()->getRadius() ;
	}

    return  particles;

}

std::vector <Inclusion3D *> GranuloBolome::operator()(bool,double rayonGranulatMax, double pourcentMasseMin, double inclusionNumber, double itzSize, bool verbose)
{
    std::vector<Inclusion3D *> rayon;
    double rayon_granulat=rayonGranulatMax;

    double pourcentMasse=1.;
    int i=0;
    if (!inclusionNumber)
        return rayon ;

		
    double remainderMass = masseInitiale;
// 	double pourcentInitial=1-exp(-c*pow(rayon_granulat,n));//calcul la constante pour que le diamètre initiale corresponde à un pourcentage de 1.
    double v = 0 ;

    rayon.push_back( new Inclusion3D(rayon_granulat+itzSize, 0, 0, 0)) ;

    double volumeGranulatsPasPlaces =0;
    while (pourcentMasse*100.>pourcentMasseMin)
    {

        double aggregateMass = densite*4./3.*M_PI*rayon_granulat*rayon_granulat*rayon_granulat; // masse plus plus gros granulat
        if (aggregateMass > masseInitiale)
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

        remainderMass -= aggregateMass;
        pourcentMasse = remainderMass/masseInitiale;

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
			if(verbose)
			{
				std::cout<<"Unused Aggregate volume "<<volumeGranulatsPasPlaces<<std::endl;
				std::cout<<"Total Aggregate Volume "<<v<<std::endl;
				std::cout << rayon.size() <<" particles" << std::endl ;
			}
            return rayon;

        }

        rayon.push_back(new Inclusion3D(rayon_granulat+itzSize, 0, 0, 0)) ;
        i++;


    }
    if(verbose)
	{
		std::cout<<"Unused Aggregate volume "<<volumeGranulatsPasPlaces<<std::endl;
		std::cout<<"Total Aggregate Volume "<<v<<std::endl;
		std::cout << rayon.size() << " particles" << std::endl ;
	}
    return  rayon;

}

GranuloFromFile::GranuloFromFile(std::string fname, std::vector<std::string> columns)
{
    std::cout << "importing file: " << fname << std::endl ;
    this->filename = fname ;
    for(int i = 0 ; i < columns.size() ; i++)
        this->fields.push_back(columns[i]) ;

    std::fstream filereader ;
    filereader.open(this->filename.c_str(),std::ios::in) ;
    char buff [256] ;
    while(!filereader.eof())
    {
        filereader >> buff ;
        this->values.push_back(atof(buff)) ;
    }
    filereader.close() ;
    std::cout << "done..." << std::endl ;
}

bool GranuloFromFile::verifyField(std::vector<std::string> columns)
{
    for(int i = 0 ; i < columns.size() ; i++)
    {
        if(this->getFieldNumber(columns[i]) == -1)
            return false ;
    }
    return true ;
}

int GranuloFromFile::getFieldNumber(std::string column)
{
        for(int j = 0 ; j < this->fields.size() ; j++)
        {
            if(strcmp(column.c_str(),fields[j].c_str()) == 0)
            {
                return j ;
            }
        }
    std::cerr << "error: field not found" << std::endl ;
    return -1 ;
}

std::vector<double> GranuloFromFile::getFieldValues(std::string column)
{
    std::vector<double> val ;
    int f = this->getFieldNumber(column) ;
    if(f > -1)
        val = this->getFieldValues(f) ;
    return val ;
}

std::vector<double> GranuloFromFile::getFieldValues(int f)
{
    std::vector<double> val ;
    int nv = this->values.size() ;
    int nf = this->fields.size() ;
    int i = 0 ;
    while(i * nf < nv)
    {
        val.push_back(values[i * nf + f]) ;
        i++ ;
    }
    return val ;
}

std::vector<Feature *> GranuloFromFile::getFeatures(TypeInclusion type, int ninc)
{
    std::vector<Feature *> inc ;
    std::vector<std::string> columns ;
    switch(type)
    {
        case 0:
            // inclusions
            columns.push_back("radius") ;
            columns.push_back("center_x") ;
            columns.push_back("center_y") ;
            break ;
        case 1:
            // inclusions 3D
            columns.push_back("radius") ;
            columns.push_back("center_x") ;
            columns.push_back("center_y") ;
            columns.push_back("center_z") ;
            break ;
        case 2:
            // ellipses
            columns.push_back("radius_a") ;
            columns.push_back("radius_b") ;
            columns.push_back("center_x") ;
            columns.push_back("center_y") ;
            columns.push_back("axis_x") ;
            columns.push_back("axis_y") ;
            break ;
    }
    std::vector< std::vector<double> > fieldvalues ;
    std::cout << "extracting columns" ;
    for(int i = 0 ; i < columns.size() ; i++)
    {
        std::cout << " ... " << columns[i] ;
        std::vector<double> val = this->getFieldValues(columns[i]) ;
        fieldvalues.push_back(val) ;
    }
    std::cout << std::endl ;
    switch(type)
    {
        case CIRCLE_INCLUSION:
            // inclusions
            std::cout << "creating inclusions..." << std::endl ;
            for(int i = 0 ; i < fieldvalues[0].size() && i < ninc ; i++)
                inc.push_back(new Inclusion(fieldvalues[0][i], fieldvalues[1][i], fieldvalues[2][i])) ;
            break ;
        case SPHERE_INCLUSION:
            // inclusions 3D
            std::cout << "creating 3D inclusions..." << std::endl ;
            for(int i = 0 ; i < fieldvalues[0].size() && i < ninc ; i++)
                inc.push_back(new Inclusion3D(fieldvalues[0][i], fieldvalues[1][i], fieldvalues[2][i], fieldvalues[3][i])) ;
            break ;
        case ELLIPSE_INCLUSION:
            // ellipses
            std::cout << "creating ellipses..." << std::endl ;
            for(int i = 0 ; i < fieldvalues[0].size() && i < ninc ; i++)
            {
                Point center(fieldvalues[2][i], fieldvalues[3][i]) ;
                Point a(fieldvalues[4][i], fieldvalues[5][i]) ;
                a *= fieldvalues[0][i] ;
                double b = fieldvalues[1][i]/fieldvalues[0][i] ;
                inc.push_back(new EllipsoidalInclusion(center,a,b)) ;
            }
            break ;
		default:
			break;
    }
    inc.pop_back() ;
    return inc ;
}

std::vector<Inclusion3D *> GranuloFromFile::getInclusion3D(int ninc)
{
    std::vector<Inclusion3D *> inc ;
    std::vector<std::string> columns ;
            columns.push_back("radius") ;
            columns.push_back("center_x") ;
            columns.push_back("center_y") ;
            columns.push_back("center_z") ;
    std::vector< std::vector<double> > fieldvalues ;
    std::cout << "extracting columns" ;
    for(int i = 0 ; i < columns.size() ; i++)
    {
        std::cout << " ... " << columns[i] ;
        std::vector<double> val = this->getFieldValues(columns[i]) ;
        fieldvalues.push_back(val) ;
    }
            std::cout << "creating 3D inclusions..." << std::endl ;
            for(int i = 0 ; i < fieldvalues[0].size() && i < ninc ; i++)
                inc.push_back(new Inclusion3D(fieldvalues[0][i], fieldvalues[1][i], fieldvalues[2][i], fieldvalues[3][i])) ;
    inc.pop_back() ;
    return inc ;
}
