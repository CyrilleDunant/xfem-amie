// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/kelvinvoight.h"
#include "../physics/maxwell.h"
#include "../physics/stiffness.h"
#include "../physics/parallel_behaviour.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../physics/homogenization/homogenization_base.h"
#include "../features/sample.h"
#include "../features/sample3d.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/features.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/itoa.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <GL/glut.h>
#include <sys/time.h>
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

Sample box(nullptr, 0.08, 0.08,0.,0.) ;


int main(int argc, char *argv[])
{
	double timestep = atof(argv[1]) ;
	int sampling = (int) atof(argv[2]) ;
	bool fullkv = (bool) atof(argv[3]) ;
	int axis = (int) atof(argv[4]) ;
	int ninc = (int) atof(argv[5]) ;
	int micro = (int) atof(argv[6]) ;
	double itz = atof(argv[7]) ;
	double aspect = 1. ;
	if(argc > 8)
		aspect = atof(argv[8]) ;
	double orientation = M_PI ;
	if(argc > 9)
		orientation = atof(argv[9]) ;
	
	if(micro > 0 && argc < 10)
		exit(0) ;
	
	FeatureTree F(&box) ;
	F.setSamplingNumber(sampling) ;

	Matrix e = (new ElasticOnlyPasteBehaviour())->param ;

	Matrix e1 = Material::cauchyGreen( 12e9, 0., true, SPACE_TWO_DIMENSIONAL) ;

	Matrix kv1e = e1*0.29 ;
	Matrix kv1t = kv1e * 10 ;

	Matrix kv2e = e1*0.25 ;
	Matrix kv2t = kv2e * 1000 ;

	Matrix mxt = e1 * 200. ;

	std::vector<std::pair<Matrix, Matrix> > branches ;
	branches.push_back(std::make_pair(kv1e, kv1t)) ;
	branches.push_back(std::make_pair(kv2e, kv2t)) ;

	Viscoelasticity * pasteKV = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, e, branches ) ;
	Viscoelasticity * pasteBurgers = new Viscoelasticity( BURGER, kv1e, kv1t, e, mxt)  ;

	if(fullkv)
		box.setBehaviour(pasteKV);
	else
		box.setBehaviour(pasteBurgers);

	GeometryType inclusions = CIRCLE ;
	if(micro == 1)
		inclusions = ELLIPSE ;
	if(micro == 2)
		inclusions = TRIANGLE ;

	ParticleSizeDistribution::get2DConcrete(&F, new ViscoElasticOnlyAggregateBehaviour(), ninc, 0.008, itz, BOLOME_A, inclusions, aspect, orientation, ninc*10000) ;
	
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(timestep) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,1)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,2)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,3)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0,4)) ;
 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0,5)) ;
	F.step() ;

	F.getAssembly()->setEpsilon(1e-14) ;
	if(axis == 1)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -10e6, 1)) ;
	else if(axis == 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_XI, RIGHT_AFTER, -10e6, 1)) ;
	F.step() ;
 	
	std::fstream out ;
	std::string name = "visco_" ;
	if(fullkv)
		name.append("kv_") ;
	if(axis == 0)
		name.append("xi_") ;
	if(axis == 1)
		name.append("eta_") ;
	if(micro == 0)
		name.append("circle_") ;
	if(micro == 1)
	{
		name.append("ellipse(") ;
		name.append(argv[6]) ;
		name.append(",") ;
		name.append(argv[7]) ;
		name.append(")_") ;
	}
	if(micro == 2)
	{
		name.append("triangle(") ;
		name.append(argv[8]) ;
		name.append(",") ;
		name.append(argv[9]) ;
		name.append(")_") ;
	}
	name.append(argv[1]) ;
	name.append("_") ;
	name.append(argv[2]) ;
	name.append("_") ;
	name.append(argv[4]) ;
	
	out.open(name.c_str(), std::ios::out) ;

	double time = timestep ;	
	Vector stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
	Vector strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
	Vector rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
	Vector disp = F.getDisplacements() ;
	out << std::setprecision(16) << time << "\t" << disp.max() << "\t" << disp.min() << "\t" <<  stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
		"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
		"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;

	if(true)
	{
		std::string nametrg = name ;
		nametrg.append("_trg_0") ;
		TriangleWriter writer(nametrg, &F, 1) ;
		writer.getField(STRAIN_FIELD) ;
		writer.getField(REAL_STRESS_FIELD) ;
		writer.getField(TWFT_STIFFNESS) ;
		writer.writeSvg(100., false) ;
	}
	
	while(time < 401)
	{
		F.step() ;

		if(time > 10)
		{
			timestep++ ;
			F.setDeltaTime(timestep) ;
		}

		time += timestep ;
		stress = F.getAverageField(REAL_STRESS_FIELD,-1,1) ;
		strain = F.getAverageField(STRAIN_FIELD,-1,1) ;
		rate = F.getAverageField(STRAIN_RATE_FIELD,-1,1) ;
		disp = F.getDisplacements() ;
		out << std::setprecision(16) << time << "\t" << disp.max() << "\t" << disp.min() << "\t" << "\t" << stress[0] << "\t" << stress[1] << "\t" << stress[2] << 
			"\t" << strain[0] << "\t" << strain[1] << "\t" << strain[2] << 
			"\t" << rate[0] << "\t" << rate[1] << "\t" << rate[2] << std::endl ;
		
		if(time == 25 || time == 55 || time == 101 || time == 200 || time == 388)
		{
			std::string nametrg = name ;
			nametrg.append("_trg_") ;
			nametrg.append(std::to_string((int) time)) ;
			TriangleWriter writer(nametrg, &F, 1) ;
			writer.getField(STRAIN_FIELD) ;
			writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.writeSvg(100., false) ;
		}
	}

		
	return 0 ;
}
