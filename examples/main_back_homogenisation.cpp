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

using namespace Amie ;


int main(int argc, char *argv[])
{
  
  std::fstream psdfile(argv[1]) ;

  double maxRadius = 38.15785 ;
  double poreFraction = 0.5;
  double minRadius = 0.1 ;
  double delta = atof(argv[2]) ;
  double factor = 1 ;
  int set = atoi(argv[3]) ;

  std::vector<Inclusion3D *> incs = PSDGenerator::get3DInclusions(maxRadius, 500.*500.*500./poreFraction, new GranuloFromCumulativePSD(argv[1], CUMULATIVE_PERCENT, factor, -1, -1), PSDEndCriteria(0.5, 0.1, 10000000) ) ;
  
  //E_anh = 3 ; E_inner = 1.5 ; E_matrix = 1	
  //E_anh = 1 ; E_inner = 1.5 ; E_matrix = 3	
  //E_anh = 1.5 ; E_inner = 3 ; E_matrix = 1

  
  double E_anyhdrous = 3 ;
  double E_inner = 1. ;
  double E_matrix = delta ;
  double nu = 0.2 ;
  if(set == 2)
  {
    E_anyhdrous = .01 ;
    E_inner = 1.5 ;
    E_matrix = 3. ; 
  }
  if(set == 3)
  {
    E_anyhdrous = 1.5 ;
    E_inner = 3 ;
    E_matrix = .01 ; 
  }
  
//   std::cout << "n = "<< incs.size() << ", maxr = " << maxRadius << ", porefrac = " << poreFraction << std::endl ;
  
  std::vector<Feature *> inner ;
  double innerVolume = 0 ;
  double innerVolumeReal = 0 ;
  double anhydrousVolume = 0 ;
  for(size_t i = 0 ; i < incs.size() ; i++)
  {
    innerVolumeReal += incs[i]->volume() ;
    inner.push_back(incs[i]);
    incs[i]->setBehaviour(new ElasticOnlyPasteBehaviour(E_inner, nu, SPACE_THREE_DIMENSIONAL)) ;
    innerVolume += incs[i]->volume() ;
  }
  
  double sampleSide = 500. ;
//   std::cout << "sampleSide = " << sampleSide << std::endl ;
  
  Sample3D s(sampleSide, sampleSide, sampleSide, 0., 0., 0.) ;
  s.setBehaviour(new ElasticOnlyPasteBehaviour(E_matrix, nu, SPACE_THREE_DIMENSIONAL) );
  
  int dummy ;
  inner = placement3D(s.getPrimitive(), inner, .05, 0, 500) ;
  
  FeatureTree ft(&s) ;
  
  for(size_t i = 0 ; i < inner.size() ; i++)
  {
    if(!i)
        ft.addFeature(&s, inner[i]);
    else
        ft.addFeature(inner[i-1], inner[i]) ;
//     ft.setSamplingFactor(inner[i], 2.) ;
//     double currentRadius = inner[i]->getRadius() ;
//     double anyhdrouRadius = currentRadius - delta ;
//     if(anyhdrouRadius > 1)
//     {
//         Inclusion3D * anh = new Inclusion3D(anyhdrouRadius, inner[i]->getCenter()) ;
//         anh->setBehaviour(new ElasticOnlyPasteBehaviour(E_anyhdrous, nu, SPACE_THREE_DIMENSIONAL)) ;
//         ft.addFeature(inner[i], anh);
//         
//         anhydrousVolume += anh->volume() ;
//         innerVolumeReal -= anh->volume() ;
//     }
  }
  
  ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, LEFT, 0.));
  ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM, 0.));
  ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, FRONT, -0.1*sampleSide));
//   ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, LEFT, 0.));
//   ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM, 0.));
  ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ZETA, BACK, 0));
  
//   ft.setOrder(LINEAR) ;
//   ft.setSamplingRestriction(SAMPLE_RESTRICT_8) ;
  ft.setSamplingNumber(1024*6);

    ft.step() ;
    Vector epsilon = ft.getAverageField(STRAIN_FIELD) ;
    Vector sigma = ft.getAverageField(REAL_STRESS_FIELD) ;
    
//     std::cout << "phi_inner = " << innerVolume/s.volume() << ", phi_anh = " << anhydrousVolume/s.volume() << ", phi_matrix = " << (s.volume()-anhydrousVolume-innerVolume)/s.volume() << std::endl ;
//     std::cout << "stress = " << sigma[0] << ", " << sigma[1] << ", " << sigma[2] << std::endl ;
//     std::cout << "strain = " << epsilon[0] << ", " << epsilon[1] << ", " << epsilon[2] << std::endl ;
    std::cout << "   "<< delta << "   " << sigma[2]/epsilon[2] << std::endl ;
    
  return 0 ;
}
