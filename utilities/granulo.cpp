//
// C++ Implementation: granulo
//
// Description:
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2007-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <cmath>
#include <vector>
#include <cstring>
#include "granulo.h"
#include "random.h"
#include <iostream>  // I/O 
#include <fstream>   // file I/O
#include "placement.h"
#include "../geometry/geometry_base.h"
#include "../features/sample.h"

using namespace Mu ;

ParticleSizeDistribution::ParticleSizeDistribution() 
{
  
}

ParticleSizeDistribution * ParticleSizeDistribution::getPSD(PSDType type)
{
	switch(type)
	{
	  case BOLOME_A:
	    return new PSDBolomeA() ;
	  case BOLOME_B:
	    return new PSDBolomeB() ;
	  case BOLOME_C:
	    return new PSDBolomeC() ;
	  case BOLOME_D:
	    return new PSDBolomeD() ;
	}
	return new ParticleSizeDistribution() ;
}
	
std::vector<Inclusion *> ParticleSizeDistribution::get2DInclusions(double rmax, double mass, PSDType type, PSDEndCriteria crit) 
{
	ParticleSizeDistribution * psd = ParticleSizeDistribution::getPSD(type) ;  
	std::vector<double> radii ;
	double diameter = rmax*2. ;
	double remainingMass = mass ;
	double remainingFraction = 1. ;
	
	while(!crit.meets(diameter*0.5, remainingFraction, radii.size()))
	{
	      radii.push_back(diameter*0.5) ;
	      remainingMass -= diameter*diameter*M_PI*1./4. ;
	      if(remainingMass < 0)
	      {
		    radii.pop_back() ;
		    break ;
	      }
	      remainingFraction = remainingMass / mass ;
	      diameter = psd->getNext2DDiameter(diameter, remainingFraction, rmax*2.) ;
	}
	crit.print(diameter*0.5, remainingFraction, radii.size()) ;
	
	std::sort(radii.begin(), radii.end()) ;
	std::reverse(radii.begin(), radii.end());
	
	std::cout << "rmin = " << radii[radii.size()-1] << "\t" << "rmax = " << radii[0] << std::endl ;
	std::cout << radii.size() << " particles generated, covering a surface of " << mass-remainingMass << std::endl ;
	
	std::vector<Inclusion *> incs ;
	for(size_t i = 0 ; i < radii.size() ; i++)
	      incs.push_back(new Inclusion(radii[i],0.,0.)) ;
	delete psd ;
	return incs ;
}

std::vector<Inclusion3D *> ParticleSizeDistribution::get3DInclusions(double rmax, double mass, PSDType type, PSDEndCriteria crit) 
{
	ParticleSizeDistribution * psd = ParticleSizeDistribution::getPSD(type) ;  
	std::vector<double> radii ;
 	double diameter = rmax*2. ;
	double remainingMass = mass ;
	double remainingFraction = 1. ;
	
	while(!crit.meets(diameter*0.5, remainingFraction, radii.size()))
	{
	      radii.push_back(diameter*0.5) ;
	      remainingMass -= diameter*diameter*diameter*M_PI*1./6. ;
	      if(remainingMass < 0)
	      {
		    radii.pop_back() ;
		    break ;
	      }
	      remainingFraction = remainingMass / mass ;
	      diameter = psd->getNext3DDiameter(diameter, remainingFraction, rmax*2.) ;
	} 
	crit.print(diameter*0.5, remainingFraction, radii.size()) ;

	std::sort(radii.begin(), radii.end()) ;
	std::reverse(radii.begin(), radii.end());
	
	std::cout << "rmin = " << radii[radii.size()-1] << "\t" << "rmax = " << radii[0] << std::endl ;
	std::cout << radii.size() << " particles generated, filling a volume of " << mass-remainingMass << std::endl ;
	
	std::vector<Inclusion3D *> incs ;
	for(size_t i = 0 ; i < radii.size() ; i++)
	      incs.push_back(new Inclusion3D(radii[i],0.,0.,0.)) ;
	delete psd ;
	return incs ;  
}

std::vector<Inclusion *> ParticleSizeDistribution::get2DMortar(double rmax, double width, size_t n, PSDType type) 
{
	return ParticleSizeDistribution::get2DInclusions(rmax, width*width*0.65, type, PSDEndCriteria(-1, 0.01, n)) ;
}

std::vector<Inclusion *> ParticleSizeDistribution::get2DConcrete(double rmax, double width, size_t n, PSDType type) 
{
	return ParticleSizeDistribution::get2DInclusions(rmax, width*width*0.8, type, PSDEndCriteria(-1, 0.01, n)) ;
}

std::vector<Inclusion *> ParticleSizeDistribution::get2DConcrete(FeatureTree * F, Form * behaviour, double rmax, size_t n, PSDType type, size_t tries, size_t seed) 
{
	Feature * box = F->getFeature(0) ;
	std::vector<Inclusion *> inc = ParticleSizeDistribution::get2DConcrete(rmax, sqrt(box->area()), n, type) ;
	std::vector<Feature *> feats ;
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		feats.push_back(inc[i]) ;
	}
	inc.clear() ;
	int nAgg = 1 ;
	srand(seed) ;
	feats = placement( dynamic_cast<Rectangle *>(box), feats, &nAgg, 0, tries ) ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		inc.push_back(dynamic_cast<Inclusion *>(feats[i])) ;
	}
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		inc[i]->setBehaviour(behaviour) ;
		F->addFeature(box, inc[i]) ;
	}
	return inc ;
}

std::vector<Inclusion *> ParticleSizeDistribution::get2DMortar(FeatureTree * F, Form * behaviour, double rmax, size_t n, PSDType type, size_t tries, size_t seed) 
{
	Feature * box = F->getFeature(0) ;
	std::vector<Inclusion *> inc = ParticleSizeDistribution::get2DMortar(rmax, sqrt(box->area()), n, type) ;
	std::vector<Feature *> feats ;
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		feats.push_back(inc[i]) ;
	}
	inc.clear() ;
	int nAgg = 1 ;
	srand(seed) ;
	feats = placement( dynamic_cast<Rectangle *>(box), feats, &nAgg, 0, tries ) ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		inc.push_back(dynamic_cast<Inclusion *>(feats[i])) ;
	}
	for(size_t i = 0 ; i < inc.size() ; i++)
	{
		inc[i]->setBehaviour(behaviour) ;
		F->addFeature(box, inc[i]) ;
	}
	return inc ;
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > ParticleSizeDistribution::get2DExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> incs, StiffnessWithImposedDeformation * behaviour, double radius, size_t n, size_t max) 
{
  	Feature * box = F->getFeature(0) ;
	Sample * sample = dynamic_cast<Sample *>(box) ;
	RandomNumber gen ;
  	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	double aggregateArea = 0 ;
	
	std::vector<ExpansiveZone *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		double w = sample->width()*0.5-radius*60 ;
		double h = sample->height()*0.5-radius*60 ;
		Point pos(gen.uniform(-w,w),gen.uniform(-h,h)) ;
		pos += sample->getCenter() ;
		bool alone  = true ;
		for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
		{
			if (squareDist(pos, zonesToPlace[j]->Circle::getCenter()) < (radius*60.+radius*60.)*(radius*60.+radius*60.))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.x, pos.y, behaviour)) ;
	}

	std::map<Inclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		for(int j = 0 ; j < incs.size() ; j++)
		{
			Circle circle(incs[j]->getRadius()*0.95 - radius*100, incs[j]->getCenter()) ;
			if(circle.in(zonesToPlace[i]->getCenter()) && incs[j]->getBoundingPoints().size() > 10)
			{
				zonesPerIncs[incs[j]]++ ; ;
				F->addFeature(incs[j],zonesToPlace[i]) ;
				ret.push_back(std::make_pair(zonesToPlace[i],incs[j])) ;
				placed = true ;
				break ;
			}
		}
		if(!placed)
			delete zonesToPlace[i] ;
		
		if(ret.size() == max)
		  break ;
	}
	
// 	exit(0) ;
	int count = 0 ;
	for(auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i)
	{
		aggregateArea+= i->first->area() ;
		count+= i->second ;
	}
	
	std::cout << ret.size() << " zones placed on reactive aggregate area of " << aggregateArea << std::endl ;
// 	std::cout << "initial Reacted Area = " << M_PI *radius *radius *ret.size() << " in " << ret.size() << " zones" << std::endl ;
// 	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;	
	
}



double ParticleSizeDistribution::getNext2DDiameter(double diameter, double, double) 
{
	return diameter ;
}

double ParticleSizeDistribution::getNext3DDiameter(double diameter, double, double) 
{
	return diameter ;
}

double PSDBolomeA::getNext2DDiameter(double diameter, double fraction, double dmax) 
{
 	double b = -(4.*fraction+1.) ;
 	double delta = 8.*fraction + 1. ;
 	return std::max(15.e-5,(- b - std::sqrt(delta))/2.*dmax/**2./M_PI*/) ;
//  	double b = 1.+fraction/.25 ;
// 	if((b - std::sqrt(b*b-fraction*fraction/(0.25*0.25)))*0.5 > 1)
// 	    std::cout << "aggregate larger than dmax!" << std::endl ;
//  	return std::max( 15e-5,dmax*(b - std::sqrt(b*b-fraction*fraction/(0.25*0.25)))*0.5) ;
}

double PSDBolomeA::getNext3DDiameter(double diameter, double fraction, double dmax) 
{
	double b = 1.+fraction/.25 ;
	return dmax*(b - std::sqrt(b*b-fraction*fraction/(0.25*0.25))) ;
}

double PSDBolomeB::getNext2DDiameter(double diameter, double fraction, double dmax) 
{
	return fraction*fraction*dmax ; 
}

double PSDBolomeB::getNext3DDiameter(double diameter, double fraction, double dmax) 
{
	return fraction*fraction*dmax ;
}

double PSDBolomeC::getNext2DDiameter(double diameter, double fraction, double dmax) 
{
	double next = fraction*fraction*dmax ;
	if(next > 0.004)
	      next *= 1.05 ;
	return next ; 
}

double PSDBolomeC::getNext3DDiameter(double diameter, double fraction, double dmax) 
{
	double next = fraction*fraction*dmax ;
	if(next > 0.004)
	      next *= 1.05 ;
	return next ; 
}

double PSDBolomeD::getNext2DDiameter(double diameter, double fraction, double dmax) 
{
	double m;//pente
	double b;//ordonnée en x=0
	if (fraction> 0. && fraction<0.1)
	{
	    Point A(.000150,0.);
	    Point B(.000315,0.1);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return std::max(15.e-05, (fraction-b)/m/**2./M_PI*/) ;
	}
	if (fraction<0.2 && fraction>= 0.1)
	{
	    Point A(.000315,0.1);
	    Point B(.00063,0.2);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m/**2./M_PI*/ ;
	}
	if (fraction<0.45 && fraction>= 0.2)
	{
	    Point A(.00063,0.2);
	    Point B(.00125,0.45);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m/**2./M_PI*/ ;
	}
	if (fraction<0.7 && fraction>= 0.45)
	{
	    Point A(.00125,0.45);
	    Point B(.0025,0.70);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m/**2./M_PI*/;
	}
	if (fraction<1. && fraction>= 0.7)
	{
	    Point A(.0025,.7);
	    Point B(.005,1.);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m/**2./M_PI*/ ;
	}
	return 15e-5 ;
}

double PSDBolomeD::getNext3DDiameter(double diameter, double fraction, double dmax) 
{
double m;//pente
	double b;//ordonnée en x=0
	if (fraction> 0. && fraction<0.1)
	{
	    Point A(.000150,0.);
	    Point B(.000315,0.1);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return std::max(7.5e-05, (fraction-b)/m);
	}
	if (fraction<0.2 && fraction>= 0.1)
	{
	    Point A(.000315,0.1);
	    Point B(.00063,0.2);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m*.5;
	}
	if (fraction<0.45 && fraction> 0.2)
	{
	    Point A(.00063,0.2);
	    Point B(.00125,0.45);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m;
	}
	if (fraction<0.7 && fraction> 0.45)
	{
	    Point A(.00125,0.45);
	    Point B(.0025,0.70);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m;
	}
	if (fraction<1. && fraction> 0.7)
	{
	    Point A(.0025,.7);
	    Point B(.005,1.);
	    m= (B.y-A.y)/(B.x-A.x);
	    b=A.y-(m*A.x);
	    return (fraction-b)/m;
	}
	return 7.5e-5 ;

}





GranuloFromFile::GranuloFromFile(std::string fname, std::vector<std::string> columns)
{
    std::cout << "importing file: " << fname << std::endl ;
    this->filename = fname ;
    for(int i = 0 ; i < columns.size() ; i++)
        this->fields.push_back(columns[i]) ;

    std::fstream filereader ;
    filereader.open(this->filename.c_str(),std::ios::in); 
    
    while(!filereader.eof())
    {
	double buff ;
        filereader >> buff ;
        this->values.push_back(buff) ;
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

std::vector<Inclusion3D *> GranuloFromFile::getInclusion3D(int ninc, double scale)
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
                inc.push_back(new Inclusion3D(fieldvalues[0][i]*scale, fieldvalues[1][i]*scale, fieldvalues[2][i]*scale, fieldvalues[3][i]*scale)) ;
    std::cout << "done" << std::endl ;
    inc.pop_back() ;
    return inc ;
}
