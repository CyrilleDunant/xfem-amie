// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2014
//
// Copyright: See COPYING file that comes with this distribution
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/stiffness.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/stiffness_with_imposed_deformation.h"
#include "../physics/diffusion.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../geometry/geometry_with_exclusion.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/itoa.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/random.h"
#include "../utilities/writer/triangle_writer.h"
#include "../utilities/writer/voxel_writer.h"


#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 

using namespace Mu ;


int main(int argc, char *argv[])
{
	
	std::vector<Inclusion *> incs =PSDGenerator::get2DInclusions(20, 100000, new GranuloFromCumulativePSD("testPSD.txt", CUMULATIVE_PERCENT), PSDEndCriteria(.01, 0., 10000) ) ;
	for(size_t i = 0 ; i < incs.size() ; i++)
		std::cout << incs[i]->getRadius() << std::endl ;
	
	exit(0) ;
	
	double npores = 1500 ;
	double poreArea = 0. ;
	std::vector<double> radii ; 
	//p=0.165, u=-5.15, s=0.3
	//p=0.205, u=-5.00, s=0.4
	//p=0.240, u=-4.80, s=0.7
	
	//P0=45 - p=0.200, u=-5.35, s=0.4
	//P0=80 - p=0.175, u=-5.20, s=0.4
	std::lognormal_distribution<double> distribution(atof(argv[1]),atof(argv[2])) ;
	std::default_random_engine generator;
	
	double poreFraction = atof(argv[3]); 
	double bulk = 13.9 ;
	double shear = 8.75 ;
	double E = .33333*(9.*bulk*shear)/(3.*bulk+shear) ;
	double nu = (3.*bulk-2.*shear)/(2.*(3.*bulk+shear)) ;
	
	for(size_t i = 0 ; i < npores ; i++)
	{
		double test = 0 ;
		do{
			test = distribution(generator) ;
		}while(test*1000 < .05 || test*1000 > 20) ;
		radii.push_back(test*1000.);
		poreArea += M_PI*radii.back()*radii.back() ;
	}
	
	std::sort(radii.begin(), radii.end()) ;
	std::reverse(radii.begin(), radii.end());
	
	std::vector<Feature *> pores ;
	for(size_t i = 0 ; i < npores ; i++)
		pores.push_back(new Pore(radii[i], 0., 0.));

	
	double sampleSide = sqrt(poreArea/poreFraction) ;
	
	Sample s(sampleSide, sampleSide, 0., 0.) ;
	s.setBehaviour(new ElasticOnlyPasteBehaviour(E, nu, SPACE_TWO_DIMENSIONAL) );
	
	pores = placement2D(s.getPrimitive(), pores, 2, 0, 100000) ;
	
	FeatureTree ft(&s) ;
	
	for(size_t i = 0 ; i < pores.size() ; i++)
		ft.addFeature(&s, pores[i]);
	
	ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, LEFT, 0.));
	ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM, 0.));
	
//	 Pour le système peu expansif P0=45.45 [MPa]
//	 Celui expansif P0=80.4 [MPa]
	
	std::vector<GeometryDefinedBoundaryCondition *> boundaryConditions ;
	for(size_t i = 0 ; i < pores.size() ; i++)
	{
		boundaryConditions.push_back(new GeometryDefinedBoundaryCondition(SET_NORMAL_STRESS, static_cast<Geometry *>(pores[i]), 0)) ;
		ft.addBoundaryCondition(boundaryConditions.back()) ;
	}

	
	ft.setOrder(LINEAR) ;
	ft.setSamplingNumber(192);
	std::fstream outfile(argv[4], std::ios::out | std::ios::app) ;
	outfile << "# p = " << atof(argv[3])<< ", u = " << atof(argv[1]) << ", s = " << atof(argv[2])  << std::endl ;
	outfile << "# P0  R_crit  Vol  dx  dy \n" << std::endl ;
	outfile.close() ;
	for(double p = 10 ; p < 110 ; p+=10)
	{
		std::fstream outfile(argv[4], std::ios::out | std::ios::app) ;
		double criticalRadius = 50 ;
		double poresUnderPressure = 0 ;
		for(size_t i = 0 ; i < pores.size() ; i++)
		{
			double pressure = p-200./(pores[i]->getRadius()-2) ;
			
			if(pressure > 0 && pores[i]->getRadius() > 2)
			{
				pressure = std::min(200./(pores[i]->getRadius()-2), p) ;
				boundaryConditions[i]->setData(pressure*1e-3) ;
				criticalRadius = std::min(pores[i]->getRadius(), criticalRadius) ;
				poresUnderPressure += pores[i]->area() ;
			}
		}
		
		ft.step() ;
		std::vector<double> apparentStrain = ft.getMacroscopicStrain(s.getPrimitive()) ;
		outfile<< p << "   " << criticalRadius << "   "<< poresUnderPressure/poreArea << "   " << apparentStrain[0] << "   " << apparentStrain[1] << std::endl ;
		outfile.close() ;
		
		MultiTriangleWriter writerm( "displacements_pores", "displacements_layer", nullptr ) ;
		writerm.reset( &ft ) ;
	// 	writerm.getField( PRINCIPAL_STRAIN_FIELD ) ;
	// 	writerm.getField( PRINCIPAL_REAL_STRESS_FIELD ) ;
		writerm.append() ;
		writerm.writeSvg(50, true) ;
		
	}
	
	exit(0) ;
		
  return 0 ;
}
