// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
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

  // scale factor (to be adjusted for the meshing)
  double scale = 1. ;
  if(argc > 1)
    scale = atof(argv[1]) ;

  // creates a 3D box of width, height and depth = 0.04, and centered on the point 0,0,0
  // (length are in meters)
  Sample boxd( nullptr, 0.02, 0.02, 0.,0.) ;
	Sample boxm( nullptr, 0.02, 0.02, 0.,0.) ;
	boxd.setBehaviour( new HydratingDiffusionCementPaste() ) ;
	

  // creates main object
  FeatureTree Fd(&boxd) ;
	Fd.addFeature(&boxd, new Pore(0.005,0.,0.));
	Fd.setInitialValue(.99);
	FeatureTree Fm(&boxm) ;
	Fm.addFeature(&boxm, new Pore(0.005,0.,0.));
	boxm.setBehaviour( new HydratingMechanicalCementPaste(&Fd) ) ;
  Fd.setOrder(LINEAR_TIME_LINEAR) ;
	
	Fm.setOrder(QUADRATIC);


  // sampling criteria
  Fd.setSamplingNumber(64); //512*16 ) ;
	Fm.setSamplingNumber(64); //512*16 ) ;


  // add boundary conditions
//   Fd.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 0));
  Fd.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, TOP, 0.7, 0));
	Fd.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, RIGHT, 0.7, 0));
// 	Fd.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_FLUX_XI, RIGHT, 7e-14, 0));

	Fm.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, LEFT, 0));
//   Fm.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI, RIGHT, 0));
	Fm.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, BOTTOM, 0));

  // assemble and solve problem
	Fd.setDeltaTime(1./24.) ;
	Fm.setDeltaTime(1./24.) ;
	srand(0) ;
	srandom(0) ;
	Fd.step();
	srand(0) ;
	srandom(0) ;
	Fm.step();
	MultiTriangleWriter writer( "saturation_field", "saturation_field_layer", nullptr ) ;
	writer.reset( &Fd ) ;
	writer.getField( TWFT_SCALAR ) ;
	writer.getField( TWFT_DOH ) ;
	writer.append() ;
	writer.writeSvg(0., true) ;

	MultiTriangleWriter writerm( "displacements", "displacements_layer", nullptr ) ;
	writerm.reset( &Fm ) ;
	writerm.getField( STRAIN_FIELD ) ;
	writerm.getField( REAL_STRESS_FIELD ) ;
	writerm.getField( TWFT_STIFFNESS ) ;
	writerm.append() ;
	writerm.writeSvg(200., true) ;
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		Fd.step();
		Fm.step();
	}

	writer.reset( &Fd ) ;
	writer.getField( TWFT_SCALAR ) ;
	writer.getField( TWFT_DOH ) ;
	writer.append() ;
	writer.writeSvg(0., true) ;

	writerm.reset( &Fm ) ;
	writerm.getField( STRAIN_FIELD ) ;
	writerm.getField( REAL_STRESS_FIELD ) ;
	writerm.getField( TWFT_STIFFNESS ) ;
	writerm.append() ;
	writerm.writeSvg(200., true) ;
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		Fd.step();
		Fm.step();
	}

	writer.reset( &Fd ) ;
	writer.getField( TWFT_SCALAR ) ;
	writer.getField( TWFT_DOH ) ;
	writer.append() ;
	writer.writeSvg(0., true) ;

	writerm.reset( &Fm ) ;
	writerm.getField( STRAIN_FIELD ) ;
	writerm.getField( REAL_STRESS_FIELD ) ;
	writerm.getField( TWFT_STIFFNESS ) ;
	writerm.append() ;
	writerm.writeSvg(200., true) ;
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		Fd.step();
		Fm.step();
	}
	
	writer.reset( &Fd ) ;
	writer.getField( TWFT_SCALAR ) ;
	writer.getField( TWFT_DOH ) ;
	writer.append() ;
	writer.writeSvg(0., true) ;

	writerm.reset( &Fm ) ;
	writerm.getField( STRAIN_FIELD ) ;
	writerm.getField( REAL_STRESS_FIELD ) ;
	writerm.getField( TWFT_STIFFNESS ) ;
	writerm.append() ;
	writerm.writeSvg(200., true) ;

	for(size_t i = 0 ; i < 3 ; i++)
	{
		Fd.step();
		Fm.step();
	}
	
	writer.reset( &Fd ) ;
	writer.getField( TWFT_SCALAR ) ;
	writer.getField( TWFT_DOH ) ;
	writer.append() ;
	writer.writeSvg(0., true) ;

	writerm.reset( &Fm ) ;
	writerm.getField( STRAIN_FIELD ) ;
	writerm.getField( REAL_STRESS_FIELD ) ;
	writerm.getField( TWFT_STIFFNESS ) ;
	writerm.append() ;
	writerm.writeSvg(200., true) ;
	
	for(size_t i = 0 ; i < 3 ; i++)
	{
		Fd.step();
		Fm.step();
	}
	
	writer.reset( &Fd ) ;
	writer.getField( TWFT_SCALAR ) ;
	writer.getField( TWFT_DOH ) ;
	writer.append() ;
	writer.writeSvg(0., true) ;

	writerm.reset( &Fm ) ;
	writerm.getField( STRAIN_FIELD ) ;
	writerm.getField( REAL_STRESS_FIELD ) ;
	writerm.getField( TWFT_STIFFNESS ) ;
	writerm.append() ;
	writerm.writeSvg(200., true) ;
  return 0 ;
}
