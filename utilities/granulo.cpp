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
#include "../features/microstructuregenerator.h"

using namespace Amie ;

PSDGenerator::PSDGenerator() 
{
  
}

// ParticleSizeDistribution * ParticleSizeDistribution::getPSD(PSDType type)
// {
// 	switch(type)
// 	{
// 	  case BOLOME_A:
// 	    return new PSDBolomeA() ;
// 	  case BOLOME_B:
// 	    return new PSDBolomeB() ;
// 	  case BOLOME_C:
// 	    return new PSDBolomeC() ;
// 	  case BOLOME_D:
// 	    return new PSDBolomeD() ;
// 	  case PSD_UNIFORM:
// 			return new ParticleSizeDistribution() ;
// 	}
// 	return new ParticleSizeDistribution() ;
// }
// 	
std::vector<Inclusion *> PSDGenerator::get2DInclusions(double rmax, double mass, ParticleSizeDistribution * type, PSDEndCriteria crit) 
{
	std::vector<double> radii ;
	double diameter = rmax*2. ;
	double remainingMass = mass ;
	double remainingFraction = 1. ;
	int tries = 0 ;
	while(!crit.meets(diameter*0.5, remainingFraction, tries++))
	{
		radii.push_back(diameter*0.5) ;
		remainingMass -= diameter*diameter*M_PI*.25 ;
		if(remainingMass < 0)
		{
			radii.pop_back() ;
			break ;
		}
		remainingFraction = remainingMass / mass ;
		diameter = type->getNext2DDiameter(diameter, remainingFraction, rmax*2.) ;
	}
// 	crit.print(diameter*0.5, remainingFraction, radii.size()) ;
	
	std::sort(radii.begin(), radii.end()) ;
	std::reverse(radii.begin(), radii.end());
	
// 	if(!radii.empty())
// 	{
// 		std::cout << "rmin = " << radii[radii.size()-1] << "\t" << "rmax = " << radii[0] << std::endl ;
// 		std::cout << radii.size() << " particles generated, covering a surface of " << mass-remainingMass << std::endl ;
// 	}
	std::vector<Inclusion *> incs ;
	for(size_t i = 0 ; i < radii.size()-1 ; i++)
	{
		incs.push_back(new Inclusion(radii[i],0.,0.)) ;
	}
	return incs ;
}

std::vector<Inclusion3D *> PSDGenerator::get3DInclusions(double rmax, double mass, ParticleSizeDistribution * type, PSDEndCriteria crit) 
{
	ParticleSizeDistribution * psd = type ;  
	std::vector<double> radii ;
 	double diameter = rmax*2. ;
	double remainingMass = mass ;
	double remainingFraction = 1. ;
	int tries = 0 ;
	
	while(!crit.meets(diameter*0.5, remainingFraction, tries++))
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
// 	crit.print(diameter*0.5, remainingFraction, radii.size()) ;

	std::sort(radii.begin(), radii.end()) ;
	std::reverse(radii.begin(), radii.end());
	
// 	std::cout << "rmin = " << radii[radii.size()-1] << "\t" << "rmax = " << radii[0] << std::endl ;
// 	std::cout << radii.size() << " particles generated, filling a volume of " << mass-remainingMass << std::endl ;
	
	std::vector<Inclusion3D *> incs ;
	for(size_t i = 0 ; i < radii.size() ; i++)
	      incs.push_back(new Inclusion3D(radii[i],0.,0.,0.)) ;
	delete psd ;
	return incs ;  
}

std::vector<Inclusion *> PSDGenerator::get2DMortar(double rmax, double width, size_t n, ParticleSizeDistribution * type) 
{
	if(!type)
		type = new PSDBolomeD() ;
	return get2DInclusions(rmax, width*width*0.65, type, PSDEndCriteria(-1, 0.01, n)) ;
}

std::vector<Inclusion *> PSDGenerator::get2DConcrete(double rmax, double width, size_t n, ParticleSizeDistribution * type, double percent) 
{
	if(!type)
		type = new PSDBolomeA() ;
	return get2DInclusions(rmax, width*width*percent, type, PSDEndCriteria(-1, 0.01, n)) ;
}

std::vector<Feature *> PSDGenerator::get2DConcrete(FeatureTree * F, Form * behaviour, size_t n, double rmax, double itz, ParticleSizeDistribution * type, GeometryType geo, double aspectRatio, double orientation, size_t tries,double fraction, Geometry * placement, std::vector<Geometry *> exclusionZones, size_t seed) 
{
	if(!type)
		type = new PSDBolomeA() ;
	Feature * box = F->getFeature(0) ;
	double ar = sqrt(box->area()) ;
	if(placement != nullptr)
		ar = sqrt(placement->area()) ;

	std::vector<Inclusion *> inc = get2DConcrete(rmax, ar, n, type,fraction) ;
	std::vector<Feature *> feats ;
	if(geo == ELLIPSE || geo == TRIANGLE)
	{
		InclusionConverter converter(geo) ;
		converter.setAspectRatio(aspectRatio) ;
//		converter.setOrientation(orientation) ;
		feats = converter.convert(inc) ;
	}
	else 
	{
		if(geo != CIRCLE)
			std::cout << "unknown geometry inclusion - fall back on circles" << std::endl ;
		for(size_t i = 0 ; i < inc.size() ; i++)
		{
			feats.push_back(inc[i]) ;
		}
	}
	inc.clear() ;
	srand(seed) ;
	if(placement)
		feats = placement2D( placement, feats, itz, 0, tries, orientation, exclusionZones ) ;
	else
		feats = placement2D( dynamic_cast<Rectangle *>(box), feats, itz, 0, tries, orientation, exclusionZones ) ;
	double area = 0 ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		area += feats[i]->area() ;
		if(behaviour)
			feats[i]->setBehaviour(behaviour) ;
// 		feats[i]->isVirtualFeature = true ;
//		feats[i]->isUpdated = false ;
		F->addFeature(box, feats[i]) ;
	}
	return feats ;
}

std::vector<Feature *> PSDGenerator::get2DEmbeddedInclusions(FeatureTree * F, Form * behaviour, std::vector<Feature *> base, size_t n, double rmax, double itz, ParticleSizeDistribution * type, GeometryType geo, double aspectRatio, double orientation, size_t tries,double fraction, Geometry * placement, std::vector<Geometry *> exclusionZones, size_t seed) 
{
	if(!type)
		type = new PSDBolomeA() ;
	Feature * box = F->getFeature(0) ;
	double ar = sqrt(box->area()) ;
	if(placement != nullptr)
		ar = sqrt(placement->area()) ;

	std::vector<Inclusion *> inc = get2DConcrete(rmax, ar, n, type,fraction) ;
	std::vector<Feature *> feats ;
	if(geo == ELLIPSE || geo == TRIANGLE)
	{
		InclusionConverter converter(geo) ;
		converter.setAspectRatio(aspectRatio) ;
		feats = converter.convert(inc) ;
	}
	else 
	{
		if(geo != CIRCLE)
			std::cout << "unknown geometry inclusion - fall back on circles" << std::endl ;
		for(size_t i = 0 ; i < inc.size() ; i++)
		{
			feats.push_back(inc[i]) ;
		}
	}
	inc.clear() ;
	srand(seed) ;
	std::vector<Geometry *> geom ;
	for(size_t i = 0 ; i < base.size() ; i++)
		geom.push_back( dynamic_cast<Geometry *>(base[i])) ;
	if(placement)
		feats = placement2DInInclusions( placement, geom, feats, itz, 0, tries, orientation, exclusionZones ) ;
	else
		feats = placement2DInInclusions( dynamic_cast<Rectangle *>(box), geom, feats, itz, 0, tries, orientation, exclusionZones ) ;
	double area = 0 ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		area += feats[i]->area() ;
		if(behaviour)
			feats[i]->setBehaviour(behaviour) ;
		for(size_t j = 0 ; j < base.size() ; j++)
		{
			if(base[j]->in(feats[i]->getCenter()))
			{
				F->addFeature(base[j], feats[i]) ;
				break ;
			}
		}
	}
	return feats ;
}
std::vector<Inclusion *> PSDGenerator::get2DMortar(FeatureTree * F, Form * behaviour, double rmax, size_t n, ParticleSizeDistribution * type, size_t tries, size_t seed) 
{
	if(!type)
		type = new PSDBolomeD() ;
	Feature * box = F->getFeature(0) ;
	std::vector<Inclusion *> inc = get2DMortar(rmax, sqrt(box->area()), n, type) ;
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

std::vector<std::pair<ExpansiveZone *, Inclusion *> > PSDGenerator::get2DExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> incs, StiffnessWithImposedDeformation * behaviour, double radius, size_t n, size_t max, int maxPerAgg) 
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
			zonesToPlace.push_back(new ExpansiveZone(nullptr, radius, pos.getX(), pos.getY(), behaviour)) ;
	}

	std::map<Inclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		if(ret.size() < max)
		{
			  for(int j = 0 ; j < incs.size() ; j++)
			  {
				  Circle circle(incs[j]->getRadius()*0.95 , incs[j]->getCenter()) ;
				  if(circle.in(zonesToPlace[i]->getCenter()))
				  {
					if((maxPerAgg < 0 || zonesPerIncs[incs[j]] < maxPerAgg) )					
					{
					  zonesPerIncs[incs[j]]++ ; 
					  F->addFeature(incs[j],zonesToPlace[i]) ;
//					  F->addPoint(&(zonesToPlace[i]->getCenter())) ;
					  ret.push_back(std::make_pair(zonesToPlace[i],incs[j])) ;
					  placed = true ;
					  break ;
					}
				  }
			  }
		}
		if(!placed)
			delete zonesToPlace[i] ;
		
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

std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion *> > PSDGenerator::get2DGrowingExpansiveZonesInAggregates(FeatureTree * F, std::vector<Inclusion *> incs, ViscoelasticityAndImposedDeformation * behaviour, Function radius, double rmax, size_t n, size_t max, int maxPerAgg) 
{
	Feature * box = F->getFeature(0) ;
	Sample * sample = dynamic_cast<Sample *>(box) ;
	RandomNumber gen ;
  	std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion *> > ret ;
	double aggregateArea = 0 ;
	
	std::vector<TimeDependentHomogenisingInclusion *> zonesToPlace ;
	
	for(size_t i = 0 ; i < n ; i++)
	{
		double w = sample->width()*0.5 - 2*rmax ;
		double h = sample->height()*0.5 - 2*rmax ;
		Point pos(gen.uniform(-w,w),gen.uniform(-h,h)) ;
		pos += sample->getCenter() ;
		bool alone  = true ;
		for(size_t j = 0 ; j< zonesToPlace.size() ; j++)
		{
			if (squareDist(pos, zonesToPlace[j]->getCenter()) < (rmax*rmax*9))
			{
				alone = false ;
				break ;
			}
		}
		if (alone)
			zonesToPlace.push_back(new TimeDependentHomogenisingInclusion(nullptr, radius, pos.getX(), pos.getY(), behaviour)) ;
	}

	std::map<Inclusion *, int> zonesPerIncs ; 
	for(size_t i = 0 ; i < zonesToPlace.size() ; i++)
	{
		bool placed = false ;
		if(ret.size() < max)
		{
			  for(int j = 0 ; j < incs.size() ; j++)
			  {
				  if(incs[j]->getRadius() > rmax*4.)
				  {
					  Circle circle(incs[j]->getRadius() - rmax*3. , incs[j]->getCenter()) ;
					  if(circle.in(zonesToPlace[i]->getCenter()))
					  {
						if((maxPerAgg < 0 || zonesPerIncs[incs[j]] < maxPerAgg) )
						{
						  zonesPerIncs[incs[j]]++ ; 
						  F->addFeature(incs[j],zonesToPlace[i]) ;
						  ret.push_back(std::make_pair(zonesToPlace[i],incs[j])) ;
						  placed = true ;
						  break ;
						}
					  }
				  }
			  }
		}
		if(!placed)
			delete zonesToPlace[i] ;
		
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
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return std::max(15.e-05, (fraction-b)/m/**2./M_PI*/) ;
	}
	if (fraction<0.2 && fraction>= 0.1)
	{
	    Point A(.000315,0.1);
	    Point B(.00063,0.2);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return (fraction-b)/m/**2./M_PI*/ ;
	}
	if (fraction<0.45 && fraction>= 0.2)
	{
	    Point A(.00063,0.2);
	    Point B(.00125,0.45);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return (fraction-b)/m/**2./M_PI*/ ;
	}
	if (fraction<0.7 && fraction>= 0.45)
	{
	    Point A(.00125,0.45);
	    Point B(.0025,0.70);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return (fraction-b)/m/**2./M_PI*/;
	}
	if (fraction<1. && fraction>= 0.7)
	{
	    Point A(.0025,.7);
	    Point B(.005,1.);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
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
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return std::max(7.5e-05, (fraction-b)/m);
	}
	if (fraction<0.2 && fraction>= 0.1)
	{
	    Point A(.000315,0.1);
	    Point B(.00063,0.2);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return (fraction-b)/m*.5;
	}
	if (fraction<0.45 && fraction> 0.2)
	{
	    Point A(.00063,0.2);
	    Point B(.00125,0.45);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return (fraction-b)/m;
	}
	if (fraction<0.7 && fraction> 0.45)
	{
	    Point A(.00125,0.45);
	    Point B(.0025,0.70);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return (fraction-b)/m;
	}
	if (fraction<1. && fraction> 0.7)
	{
	    Point A(.0025,.7);
	    Point B(.005,1.);
	    m= (B.getY()-A.getY())/(B.getX()-A.getX());
	    b=A.getY()-(m*A.getX());
	    return (fraction-b)/m;
	}
	return 7.5e-5 ;

}


GranuloFromFile::GranuloFromFile(const std::string & fname, std::vector<std::string> columns)
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



GranuloFromCumulativePSD::GranuloFromCumulativePSD(const std::string & filename, PSDSpecificationType t, double factor, double cutOffUp, double cutOffDown) 
{ 
	std::fstream file(filename) ;

	if(!file.is_open())
	{
		std::cout << "file " << filename << " doesn't exist!" << std::endl ;
		exit(0) ;
	}
	bool foundZeroFrac = false ;
	double zeroFracRad = 0 ;
	double maxfrac = 0 ;
	double minfrac = 100 ;
	do {
		double frac ;
		double rad ;
		file >> frac >> rad ;
		
		if((cutOffUp > 0 && rad*factor <= cutOffUp || cutOffUp < 0)  && 
			 (cutOffDown > 0 && rad*factor >= cutOffDown || cutOffDown < 0) 
		)
		{
// 			if(frac < 1e-13 )
// 			{
// 				foundZeroFrac = true ;
// 				zeroFracRad = std::max(rad*factor, zeroFracRad) ;
// 				minfrac = std::min(minfrac, frac) ;
// 			}
// 			else
// 			{
				fraction.push_back(frac);
				radius.push_back(rad*factor);
				maxfrac = std::max(maxfrac, frac) ;
				minfrac = std::min(minfrac, frac) ;
// 			}
		}
	}while(!file.eof()) ;
	
	fraction.pop_back();
	radius.pop_back();
	
	//collate possible list of zero-radius or fraction
// 	if(radius.front() > radius.back())
// 	{
// 		fraction.push_back(0);
// 		radius.push_back(zeroFracRad);
// 	}
// 	else
// 	{
// 		fraction.insert(fraction.begin(), 0.);
// 		radius.insert(radius.begin(), zeroFracRad) ;
// 	}
	
	if(cutOffUp > 0  || cutOffDown > 0)
	{
		if(t == CUMULATIVE_PERCENT)
			t = CUMULATIVE_ABSOLUTE ;
		if(t == CUMULATIVE_PERCENT_REVERSE)
			t = CUMULATIVE_ABSOLUTE_REVERSE ;
	}
	
	for(size_t i = 0 ; i < fraction.size() ; i++)
	{
		fraction[i] -= minfrac ;
	}
	
	switch(t)
	{
		case CUMULATIVE_PERCENT:
		{
			for(size_t i = 0 ; i < fraction.size() ; i++)
				fraction[i] /= 100. ;
			break ;
		}
		case CUMULATIVE_FRACTION:
		{
			break ;
		}
		case CUMULATIVE_ABSOLUTE:
		{
			double maxv = std::max(fraction.front(), fraction.back()) ;
			for(size_t i = 0 ; i < fraction.size() ; i++)
				fraction[i] /= maxv ;
			break ;
		}
		case CUMULATIVE_PERCENT_REVERSE:
		{
			for(size_t i = 0 ; i < fraction.size() ; i++)
				fraction[i] = (100. - fraction[i])/100. ;
			break ;
		}
		case CUMULATIVE_FRACTION_REVERSE:
		{
			for(size_t i = 0 ; i < fraction.size() ; i++)
				fraction[i] = 1. - fraction[i] ;
			break ;
		}
		case CUMULATIVE_ABSOLUTE_REVERSE:
		{
			double maxv = std::max(fraction.front(), fraction.back()) ;
			for(size_t i = 0 ; i < fraction.size() ; i++)
			{
				fraction[i] /= maxv ;
				fraction[i] = 1. - fraction[i] ;
			}
			break ;
		}
	}
	if(fraction.back() < fraction.front())
	{
		std::reverse(fraction.begin(), fraction.end());
		std::reverse(radius.begin(), radius.end());
	}
// 	for( size_t i = 0  ; i < radius.size() ; i++)
// 		std::cout << fraction[i] << "  " << radius[i] << std::endl ;
// 	
// 	exit(0) ;
	
}

double GranuloFromCumulativePSD::getNext2DDiameter(double diameter, double frac, double dmax)
{

	if(diameter < 2.*radius.front()+POINT_TOLERANCE_2D)
	{
		return 2.*radius.front() ;
	}
		
	for( size_t i = 0  ; i < radius.size()-1 ; i++)
	{
		if(fraction[i] > frac )
		{
			double df = (exp(frac)-exp(fraction[i-1]))/(exp(fraction[i])-exp(fraction[i-1])) ;
			double v = std::max(2.*(radius[i-1]*(1.-df) + radius[i]*df), 2.*radius.front())  ;
			return v ;
		}
	}
	double df = (exp(frac)-exp(fraction[radius.size()-2]))/(exp(fraction[radius.size()-1])-exp(fraction[radius.size()-2])) ;
	double v = std::max(2.*(radius[radius.size()-2]*(1.-df) + radius[radius.size()-1]*df), 2.*radius.front())  ;
	return v ;

	
}
double GranuloFromCumulativePSD::getNext3DDiameter(double diameter, double frac, double dmax)
{
	if(diameter < 2.*radius.front()+POINT_TOLERANCE_2D)
	{
		return 2.*radius.front() ;
	}
		
	for( size_t i = 0  ; i < radius.size()-1 ; i++)
	{
		if(fraction[i] > frac )
		{
			double df = (frac-fraction[i-1])/(fraction[i]-fraction[i-1]) ;
			double v = std::max(2.*(radius[i-1]*(1.-df) + radius[i]*df), 2.*radius.front())  ;
			return v ;
		}
	}
	double df = (frac-fraction[radius.size()-2])/(fraction[radius.size()-1]-fraction[radius.size()-2]) ;
	double v = std::max(2.*(radius[radius.size()-2]*(1.-df) + radius[radius.size()-1]*df), 2.*radius.front())  ;
	return v ;
}

