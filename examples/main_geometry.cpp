// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "main.h"
#include "../features/features.h"
#include "../physics/stiffness.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../features/sample3d.h"
#include "../features/inclusion3d.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/writer/voxel_writer.h"


#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h> 

using namespace Amie ;


int main(int argc, char *argv[])
{	
	
	VirtualMachine vm ;
	
	Function f("2 2 + x * x *") ;
	
	std::cout << vm.deval(f,XI ,.1) << std::endl ;
	
	
	
	exit(0) ;
	
	
  // scale factor (to be adjusted for the meshing)
  
	double nu = 0.2 ;
	double E = 1 ;
	
	Matrix m0(6,6) ;
	m0[0][0] = 1. - nu ; m0[0][1] = nu ; m0[0][2] = nu ;
	m0[1][0] = nu ; m0[1][1] = 1. - nu ; m0[1][2] = nu ;
	m0[2][0] = nu ; m0[2][1] = nu ; m0[2][2] = 1. - nu ;
	m0[3][3] = 0.5 - nu ;
	m0[4][4] = 0.5 - nu ;
	m0[5][5] = 0.5 - nu ;
	m0 *= E/((1.+nu)*(1.-2.*nu)) ;
	
  // creates a 3D box of width, height and depth = 0.04, and centered on the point 0,0,0
  // (length are in meters)
  Sample3D box( nullptr, 0.04, 0.04, 0.04, 0.,0.,0.) ;
  
  // creates 3D inclusions with rmax = 0.002, covering a volume of 0.000252, using BOLOME type B particle sizs distribution
  // psd ends when smallest radius reaches 0.00025 or when 1000 inclusions have been generated
  // covering volume 0.000042
  //std::vector<Inclusion3D *> incs = ParticleSizeDistribution::get3DInclusions(0.002, 0.000042, BOLOME_B, PSDEndCriteria(0.00025, -1, 3200)) ;
  std::vector<Inclusion3D *> incs = PSDGenerator::get3DInclusions(0.002, 0.000023104, new PSDBolomeB(), PSDEndCriteria(0.00025, -1, 3200)) ;
  // attributes mechanical behaviour to the box and the aggregates
  
	box.setBehaviour( new ElasticOnlyPasteBehaviour( 12e9, 0.3, SPACE_THREE_DIMENSIONAL ) ) ;
  for(size_t i = 0 ; i < incs.size() ; i++)
    incs[i]->setBehaviour( new ElasticOnlyAggregateBehaviour( 60e9, 0.3, SPACE_THREE_DIMENSIONAL ) ) ;

  // creates main object
  FeatureTree F(&box) ;
  F.setOrder(LINEAR) ;

  // place the inclusions in the box
  std::vector<Feature *> feats ;
  for(size_t i = 0 ; i < incs.size() ; i++)
    feats.push_back(incs[i]) ;
  feats = placement3D(box.getPrimitive(), feats, 0.) ;
  incs.clear() ;
	
	
  for(size_t i = 0 ; i < feats.size() ; i++)
    {
      incs.push_back( dynamic_cast<Inclusion3D *>( feats[i] ) ) ;
      F.addFeature( &box, incs[i] ) ;
    }
	
  // sampling criteria
  F.setSamplingNumber(32000); //512*16 ) ;

  std::ofstream inclusions;
  inclusions.open("inclusions.csv");
  auto inc_it = incs.begin();
  auto inc_end = incs.end();
  inclusions << incs.size() << std::endl;
  inclusions  <<  "#radius,x,y,z" << std::endl;
  for (; inc_it != inc_end; ++inc_it) {
    Point & center = (*inc_it)->getCenter();
    inclusions  << (*inc_it)->getRadius() << "," << center.getX() << "," << center.getY() << "," << center.getZ()  << std::endl;
  }
  inclusions.close();
	
  // generate elements and get list of tetrahedrons
  
  unsigned int node_counter = 1;
  typedef std::map<Point *, unsigned int> nodes_map;
  nodes_map node_index_mapping;
  std::vector<Point *> nodes;

  std::ofstream mesh_file; mesh_file.open("mesh.msh");
  mesh_file << "$MeshFormat" << std::endl
	    << "2.2 0 " << sizeof(double) << std::endl
	    << "$EndMeshFormat" << std::endl;

  // tets[i]->getBoundingPoints().size();
  //   tets[i]->getBoundingPoints(j).getId();
  // tets[i]->getBoundingPoints().size()

  for (auto it = F.get3DMesh()->begin() ;it != F.get3DMesh()->begin(); it++) {
    for (unsigned int n = 0; n < it->size(); ++n) {
      nodes_map::iterator pid = node_index_mapping.find(&it->getBoundingPoint(n));
      if(pid == node_index_mapping.end()) {
	nodes.push_back(&it->getBoundingPoint(n));
	node_index_mapping[&it->getBoundingPoint(n)] = node_counter;
	node_counter++;
      }
    }
  }

  std::vector<Point *>::iterator point_it  = nodes.begin();
  std::vector<Point *>::iterator point_end = nodes.end();

  mesh_file << "$Nodes" << std::endl
	    << nodes.size() << std::endl;
  for (unsigned int n = 1; point_it != point_end; ++point_it, ++n) {
    Point & p = **point_it;
    mesh_file << n << " " << p.getX() << " " << p.getY() << " " << p.getZ() << std::endl;
  }
  mesh_file << "$EndNodes" << std::endl;

  mesh_file << "$Elements" << std::endl
	    << F.get3DMesh()->begin().size() << std::endl;

  auto it = F.get3DMesh()->begin();
  for (unsigned int t = 1; it != F.get3DMesh()->end(); ++it, ++t) {
    unsigned int msh_type = 0;
    if(it->size() == 4) msh_type = 4;
    else if(it->size() == 10) msh_type = 11;
    assert(msh_type != 0);

    mesh_file << t << " " << msh_type << " 2 0 0";
    for (unsigned int n = 0; n < it->size(); ++n) {
      nodes_map::iterator pid = node_index_mapping.find(&it->getBoundingPoint(n));
      mesh_file << " " << pid->second;
    }
    mesh_file << std::endl;
  }
  

  mesh_file << "$EndElements" << std::endl;

  mesh_file.close();

  // get nodes
  //std::vector<Point> nodes = F.getNodes() ;

  // add boundary conditions
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM));
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA, TOP,0.0001));
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT));
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, FRONT));

  // assemble and solve problem
  F.step();
	
	
	VoxelWriter vw("sphere_stiffness_0.05M_nomin", 200) ;
	vw.getField(&F, VWFT_STIFFNESS) ;
	vw.write();
  // calculate averaged stress and strain
	Vector x = F.
	
	
	
	getDisplacements() ;
  Vector strain = F.getAverageField(STRAIN_FIELD);
  Vector stress = F.getAverageField(REAL_STRESS_FIELD);

  // write results in textfile
  
  std::ofstream strain_file; strain_file.open("strain.getT()xt");
  strain_file << "STRAIN_FIELD: " << std::endl;

  for (size_t i=0; i<strain.size(); ++i){
    strain_file << strain[i];
    strain_file << std::endl;
  } 

  strain_file.close();

  std::ofstream stress_file; stress_file.open("stress.getT()xt");
  stress_file << "REAL_STRESS_FIELD: " << std::endl;

  for (size_t i=0; i<stress.size(); ++i){
    stress_file << stress[i];
    stress_file << std::endl;
  } 

  stress_file.close();


  return 0 ;
}
