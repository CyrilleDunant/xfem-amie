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

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <time.h> 
#define DEBUG 

#define ID_QUIT 1
#define ID_ZOOM 5
#define ID_UNZOOM 6
#define ID_NEXT10 7
#define ID_NEXT100 3
#define ID_NEXT1000 4
#define ID_NEXT 2
#define ID_NEXT_TIME 0
#define ID_REFINE 8
#define ID_AMPLIFY 9
#define ID_DEAMPLIFY 10

#define ID_DISP 11
#define ID_STRAIN_XX 12
#define ID_STRAIN_XY 13
#define ID_STRAIN_YY 14
#define ID_STRESS_XX 15
#define ID_STRESS_XY 16
#define ID_STRESS_YY 17
#define ID_STIFNESS 18
#define ID_ELEM 19
#define ID_VON_MISES 20
#define ID_ANGLE 22
#define ID_ENRICHMENT 21

#define DISPLAY_LIST_DISPLACEMENT 1
#define DISPLAY_LIST_ELEMENTS 2
#define DISPLAY_LIST_STRAIN_XX 3
#define DISPLAY_LIST_STRAIN_YY 4
#define DISPLAY_LIST_STRAIN_XY 5
#define DISPLAY_LIST_STRESS_XX 6
#define DISPLAY_LIST_STRESS_YY 7
#define DISPLAY_LIST_STRESS_XY 8
#define DISPLAY_LIST_CRACK 9
#define DISPLAY_LIST_STIFFNESS 10
#define DISPLAY_LIST_VON_MISES 11
#define DISPLAY_LIST_ANGLE 23
#define DISPLAY_LIST_ENRICHMENT 12
#define DISPLAY_LIST_STIFFNESS_DARK 24

using namespace Mu ;

FeatureTree * featureTree ;
std::vector<DelaunayTriangle *> triangles ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double x_max = 0 ;
double y_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;

double timepos = 0.00 ;
double percent = 0.10 ;
double placed_area = 0 ;

double stress = 15e6 ;

Sample sample(nullptr, 0.004, 0.004, 0.0, 0.0) ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

std::pair<std::vector<Inclusion * >, std::vector<Pore * > > i_et_p ;

std::vector<std::pair<ExpansiveZone *, EllipsoidalInclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;

Vector b(0) ;
Vector x(0) ;
Vector sigma(0) ; 
Vector sigma11(0) ; 
Vector sigma22(0) ; 
Vector sigma12(0) ; 
Vector epsilon(0) ; 
Vector epsilon11(0) ; 
Vector epsilon22(0) ; 
Vector epsilon12(0) ; 
Vector vonMises(0) ; 
Vector angle(0) ; 


int main(int argc, char *argv[])
{	
	// scale factor (to be adjusted for the meshing)
	double scale = 1. ;
	if(argc > 2)
		scale = atof(argv[2]) ;
  
	// creates a 3D box of width, height and depth = 0.04, and centered on the point 0,0,0
	// (length are in meters)
	Sample3D box( nullptr, 0.04*scale, 0.04*scale, 0.04*scale, 0.,0.,0.) ;
  
	// creates 3D inclusions with rmax = 0.002, covering a volume of 0.000252, using BOLOME type B particle sizs distribution
	// psd ends when smallest radius reaches 0.00025 or when 1000 inclusions have been generated
	std::vector<Inclusion3D *> incs = ParticleSizeDistribution::get3DInclusions(0.002*scale, 0.000252*scale, BOLOME_B, PSDEndCriteria(0.00025*scale, -1, 1000)) ;
  
	// attributes mechanical behaviour to the box and the aggregates
	box.setBehaviour( new ElasticOnlyPasteBehaviour( 12e9, 0.3, SPACE_THREE_DIMENSIONAL ) ) ;
	for(size_t i = 0 ; i < incs.size() ; i++)
		incs[i]->setBehaviour( new ElasticOnlyAggregateBehaviour( 60e9, 0.3, SPACE_THREE_DIMENSIONAL ) ) ;

	// creates main object
	FeatureTree F(&box) ;

	// place the inclusions in the box
	std::vector<Feature *> feats ;
	for(size_t i = 0 ; i < incs.size() ; i++)
		feats.push_back(incs[i]) ;
	int n = 0 ;
	feats = placement(dynamic_cast<Hexahedron *>(&box), feats, &n) ;
	incs.clear() ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	{
		  incs.push_back( dynamic_cast<Inclusion3D *>( feats[i] ) ) ;
		  F.addFeature( &box, incs[i] ) ;
	}
	
	// sampling criteria (to be adjusted for the number of points you need)
	F.setSamplingNumber( atof(argv[1]) ) ;
	
	// generate elements and get list of tetrahedrons
	std::vector<DelaunayTetrahedron *> tets = F.getElements3D() ;
	// get nodes
	std::vector<Point> nodes = F.getNodes() ;
	
	F.step() ;
	
	return 0 ;
}
