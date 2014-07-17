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
  
  std::fstream psdfile(argv[1]) ;
  double thresholdRadius = 30. ;
  if(argc > 3)
    thresholdRadius = std::max(atof(argv[3]), 30.) ;
  double maxRadius ;
  double poreFraction ;
  double waterLayerDepth = 1.5 ;
  double factor = 1000. ;
  double minRadius = waterLayerDepth*.25 ;
	double shapeFactor = 200 ;
  
  psdfile >> poreFraction >>  maxRadius;
  while(!psdfile.eof())
  {
    double dummy ;
    psdfile >> poreFraction >> dummy ;
  }
  poreFraction /= 100. ;
  maxRadius = std::min(maxRadius*factor, thresholdRadius) ;
  
  std::vector<Inclusion *> incs = PSDGenerator::get2DInclusions(maxRadius, maxRadius*maxRadius*M_PI*60., new GranuloFromCumulativePSD(argv[1], CUMULATIVE_ABSOLUTE_REVERSE, 1000., maxRadius, minRadius), PSDEndCriteria(.1, 0., 50000) ) ;
  double bulk = 13.9 ;
  double shear = 8.75 ;
  double E = .33333*(9.*bulk*shear)/(3.*bulk+shear) ;
  double nu = (3.*bulk-2.*shear)/(2.*(3.*bulk+shear)) ;
  
  
  std::cout << "maxr = " << maxRadius << ", porefrac = " << poreFraction << std::endl ;
  
  std::vector<Feature *> pores ;
  double poreArea = 0 ;
  for(size_t i = 0 ; i < incs.size() ; i++)
  {
    pores.push_back(incs[i]);
    incs[i]->setBehaviour(new VoidForm()) ;
    poreArea += incs[i]->area() ;
  }

  
  double sampleSide = sqrt(poreArea/poreFraction) ;
  
  Sample s(sampleSide, sampleSide, 0., 0.) ;
  s.setBehaviour(new ElasticOnlyPasteBehaviour(E, nu, SPACE_TWO_DIMENSIONAL) );
  
  pores = placement2D(s.getPrimitive(), pores, waterLayerDepth*4., 0, 2000000) ;
  
  FeatureTree ft(&s) ;
  
  for(size_t i = 0 ; i < pores.size() ; i++)
    ft.addFeature(&s, pores[i]);
  
  ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, LEFT, 0.));
  ft.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM, 0.));
  
//   Pour le système peu expansif P0=45.45 [MPa]
//   Celui expansif P0=80.4 [MPa]
  
  std::vector<GeometryDefinedBoundaryCondition *> boundaryConditions ;
  for(size_t i = 0 ; i < pores.size() ; i++)
  {
    boundaryConditions.push_back(new GeometryDefinedBoundaryCondition(SET_NORMAL_STRESS, static_cast<Geometry *>(pores[i]), 0)) ;
    ft.addBoundaryCondition(boundaryConditions.back()) ;
  }

  
  ft.setOrder(LINEAR) ;
  ft.setSamplingNumber(192);
  std::fstream outfile(argv[2], std::ios::out | std::ios::app) ;
  outfile << "# psd = " << argv[1] << ", outfile = " << argv[2]  << std::endl ;
  outfile << "# P0  R_crit  Vol  dx  dy \n" << std::endl ;
  outfile.close() ;
  for(double p = 10 ; p < 110 ; p+=10)
  {
    std::fstream outfile(argv[2], std::ios::out | std::ios::app) ;
    double criticalRadius = thresholdRadius ;
    double poresUnderPressure = 0 ;
    for(size_t i = 0 ; i < pores.size() ; i++)
    {
      double pressure = p-shapeFactor/(pores[i]->getRadius()-waterLayerDepth) ;
      
      if(pressure > 0 && pores[i]->getRadius() > waterLayerDepth)
      {
        pressure = std::min(shapeFactor/(pores[i]->getRadius()-waterLayerDepth), p) ;
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
  //   writerm.getField( PRINCIPAL_STRAIN_FIELD ) ;
    writerm.getField( PRINCIPAL_REAL_STRESS_FIELD ) ;
    writerm.append() ;
    writerm.writeSvg(50, true) ;
    
  }
    
  return 0 ;
}
