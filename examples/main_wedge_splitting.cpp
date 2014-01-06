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
#include "../physics/stiffness_and_fracture.h"
#include "../physics/parallel_behaviour.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../physics/fracturecriteria/creeprupture.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/damagemodels/fiberbasedisotropiclineardamage.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/viscoelasticity.h"
#include "../physics/viscoelasticity_and_fracture.h"
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

#include <omp.h>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Mu ;

// see Cecot 2001, p98, for sample dimensions

double length = 0.2 ; //5.5 ;
double width = 0.06 ;
double depth = 0.03 ;//2.5 ;
double nnotch = 0.085 ;
double nwidth = 0.0025 ;
Sample box(nullptr, length, length,0.,length*0.5) ;
Sample top(nullptr, width, depth, 0.0, length - depth*0.5) ;
Sample notch(nullptr,nwidth, nnotch, 0.0, length - depth - nnotch*0.5) ;

int main(int argc, char *argv[])
{
//	 omp_set_schedule(omp_sched_static, 60) ;
//	omp_set_num_threads(1) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(96) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setMaxIterationsPerStep( 100 ) ;
	double totaltime = atof(argv[1]) ;
	F.setDeltaTime(totaltime/100.) ;
	F.setMinDeltaTime((totaltime/100.)*1e-9) ;
	size_t seed = 0 ;
	ViscoElasticOnlyPasteBehaviour pastenodamage(20e9) ;

	Rectangle * placement= new Rectangle(width,nnotch*1.1, 0., nnotch*0.552) ;

	Sample left(nullptr, (length-width)*.5, length, (length-width)*.25-length*.5, 0) ;
	left.isVirtualFeature  = true ;
	left.setBehaviour(&pastenodamage) ;
	
	Sample middle(nullptr, length, depth+nnotch*0.9, 0., length-(depth+nnotch*0.9)*0.5) ;
	middle.isVirtualFeature = true ;
	middle.setBehaviour(&pastenodamage) ;

	Sample right(nullptr, (length-width)*.5, nnotch*2., length*.5-(length-width)*.25, nnotch) ;
	right.isVirtualFeature  = true ;
	right.setBehaviour(&pastenodamage) ;

	Rectangle refinement( width*1.1,nnotch*1.2, 0., nnotch*0.6) ;
	Rectangle refinement2( width*1,nnotch*1.1, 0., nnotch*0.55) ;
	F.addRefinementZone( &refinement ) ;
//	F.addRefinementZone( &refinement2 ) ;
	
	ViscoDamagePasteBehaviour paste ;
	paste.freeblocks = 0 ;
	if(argv[2] == std::string("strain"))
	{

		paste.up = 0.00035 ;
	}
	if(argv[2] == std::string("stress"))
	{
		paste.ctype = STRESS_CRITERION ;
		paste.up = 0.0003 ;
	}
	if(argv[2] == std::string("mixed"))
	{
		paste.ctype = MIXED_CRITERION ;
		paste.up = 0.00035 ;
		paste.stressFraction = 0.75 ;
	}
	paste.materialRadius = 0.002 ;

	if(argv[2] == std::string("elastic"))
		box.setBehaviour( &pastenodamage ) ;
	else
		box.setBehaviour( &paste ) ;

	top.setBehaviour( new VoidForm() ) ;
	notch.setBehaviour( new VoidForm() ) ;

	ViscoElasticOnlyAggregateBehaviour agg ;
	agg.freeblocks = 0 ;

	std::vector<Geometry *> exclusionZones ;
	exclusionZones.push_back( notch.getPrimitive() ) ;

	std::vector<Feature *> aggregates = ParticleSizeDistribution::get2DConcrete( &F, &agg, 300, 0.008, 0.0003, BOLOME_A, CIRCLE, 1., M_PI, 100000, 0.8, placement, exclusionZones ) ;
	for(size_t i = 30 ; i < aggregates.size() ; i++)
	{
		F.setSamplingFactor(aggregates[i], 1.5) ;
	}

 	F.addFeature(&box, &top) ;
	F.addFeature(&box, &middle) ;
 	F.addFeature(&middle, &notch) ;
	F.setSamplingFactor(&box, .125) ;
//	F.setSamplingFactor(&top, 10.) ;
//	F.setSamplingFactor(&notch, 10.) ;
	F.addFeature(&box, &left) ;
	F.addFeature(&box, &right) ;
	F.setSamplingRestriction( SAMPLE_RESTRICT_4 ) ;

// 	F.addBoundaryCondition(new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, -1e-5, 1e-5, -1e-5, nnotch+1e-5,0, 0 )) ;
 //	F.addBoundaryCondition(new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, -1e-5, 1e-5, -1e-5, nnotch+1e-5, 0, 2 )) ;
 //	F.addBoundaryCondition(new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, -1e-5, 1e-5, -1e-5, nnotch+1e-5, 0, 4 )) ;
// 	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 6 )) ;
 //	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 8 )) ;

	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 1 ) ) ;
 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 3 ) ) ;
 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 5 ) ) ;
	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 0 ) ) ;
 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 2 ) ) ;
 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., 0.), 0, 4 ) ) ;
// 	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., -length*0.5), 0, 7 ) ) ;
 //	F.addBoundaryCondition( new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, Point(0., -length*0.5), 0, 9 ) ) ;

//	BoundingBoxAndRestrictionDefinedBoundaryCondition * disp = new BoundingBoxAndRestrictionDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, -0.1*width, width*0.6, length*0.99, length*1.01, 0., 0 ) ;
	BoundingBoxNearestNodeDefinedBoundaryCondition * dispxr = new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, Point(width*0.5, length), 0., 0) ;
	BoundingBoxNearestNodeDefinedBoundaryCondition * dispyr = new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, Point(width*0.5, length), 0., 1) ;
	BoundingBoxNearestNodeDefinedBoundaryCondition * dispxl = new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, Point(width*(-0.5), length), 0., 0) ;
	BoundingBoxNearestNodeDefinedBoundaryCondition * dispyl = new BoundingBoxNearestNodeDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, Point(width*(-0.5), length), 0., 1) ;
	F.addBoundaryCondition(dispxr) ;
	F.addBoundaryCondition(dispyr) ;
	F.addBoundaryCondition(dispxl) ;
	F.addBoundaryCondition(dispyl) ;

	
	F.step() ;
	F.getAssembly()->setEpsilon(1e-8) ;
 	Vector x = F.getAverageField(STRAIN_FIELD) ;
 	Vector y = F.getAverageField(REAL_STRESS_FIELD) ;
//	std::cout << 0. << "\t" << x[0] << "\t" << y[0] << std::endl ;

	size_t i = 0 ;
	size_t j = 0 ;

	std::fstream out ;
	std::string tata = "wedge_" ;
	tata.append(argv[2]) ;
	tata.append("_") ;
	tata.append(argv[1]) ;
	tata.append(".txt") ;
	out.open(tata.c_str(), std::ios::out) ;

	double speed = 0.0005/totaltime ;
	double totaldisp = 0.;
	bool goOn = true ;
	std::vector<DelaunayTriangle *> trg = F.getElements2D() ;

	while(speed*F.getCurrentTime() <= 0.0005)
	{
		if(goOn)
		{
			i++ ;
			totaldisp = speed * F.getCurrentTime() ;
			dispxr->setData( totaldisp*0.5 ) ;
			dispyr->setData( totaldisp*0.25/(-1.8666) ) ;
			dispxl->setData( totaldisp*(-0.5) ) ;
			dispyl->setData( totaldisp*(0.25)/(-1.8666) ) ;
			F.setDeltaTime(totaltime/100.) ;
			F.setMinDeltaTime((totaltime/100.)*1e-9) ;
			std::cout << totaldisp << std::endl ;
		}

		goOn = F.step() ;

		if(goOn)
		{
			std::string tati = tata ; tati.append("_") ;
			tati.append(itoa(i)) ;
			std::cout << tati << std::endl ;
			TriangleWriter writer(tati, &F, 1) ;
	 		writer.getField(STRAIN_FIELD) ;
	 		writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_DAMAGE) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.write() ;
			j = 0 ;
		}
		if(!goOn)
		{
			j++ ;
			std::string tati = tata ; tati.append("_") ;
			tati.append(itoa(i)) ; tati.append("_inter_") ;
			tati.append(itoa(j)) ;
			std::cout << tati << std::endl ;
			TriangleWriter writer(tati, &F, -1) ;
	 		writer.getField(STRAIN_FIELD) ;
	 		writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_DAMAGE) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.write() ;

		}

 		x = F.getAverageField(STRAIN_FIELD, -1, -1+2*goOn) ;
 		y = F.getAverageField(REAL_STRESS_FIELD, -1, -1+2*goOn) ;
		out << trg[0]->getBoundingPoint(0+3*goOn).t << "\t" << totaldisp << "\t" << x[0] << "\t" << y[0] << "\t" << F.averageDamage << std::endl ;

	}
	
//	F.getAssembly()->print() ;
	
	return 0 ;
}

