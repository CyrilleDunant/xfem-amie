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
#include "../physics/dual_behaviour.h"
#include "../physics/parallel_behaviour.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/orthotropicstiffness.h"
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

Sample box(nullptr, 0.08, 0.08,0.0,0.0) ;

int findIndexOfAggregate( DelaunayTriangle * tri, std::vector<Inclusion *> aggregates, FeatureTree * f)
{
	for(size_t i = 0 ; i < aggregates.size() ; i++)
	{
		std::vector<DelaunayTriangle *> aggTriangles = aggregates[i]->getElements2D(f) ;
		if( std::find( aggTriangles.begin(), aggTriangles.end(), tri) != aggTriangles.end() )
			return i ;
	}
  
	return 0 ;
}



int main(int argc, char *argv[])
{
//	omp_set_num_threads(1) ;

	double timeScale = 450. ;//atof(argv[1]) ;
	double tau = 30. ;
	double timestep = tau ;
	double appliedLoadEta = 0 ;//atof(argv[2])*(-1e6) ;
	double appliedLoadXi = 0 ;//atof(argv[3])*(-1e6) ;
	double degree = atof(argv[1]) ;
	tau /= degree/0.001 ;
		
	int nzones = 400 ;
	int naggregates = 1000 ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(160) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;
	F.setMinDeltaTime(tau*1e-9) ;
	F.setMaxIterationsPerStep(5000) ;
		
//	box.setBehaviour( new ViscoElasticOnlyPasteBehaviour() );
	box.setBehaviour( new ViscoDamagePasteBehaviour() );
	std::vector<Feature *> feats = PSDGenerator::get2DConcrete( &F, new ViscoDamageAggregateBehaviour(), naggregates, 0.008, 0.0001, new PSDBolomeA()) ;
	
	std::vector<Inclusion *> aggregates ;
	for(size_t i = 0 ; i < feats.size() ; i++)
	      aggregates.push_back( dynamic_cast<Inclusion *>(feats[i]) ) ;
	
	double aggregatesArea = 0 ; 
	for(size_t i = 0 ; i < aggregates.size() ; i++)
		aggregatesArea += aggregates[i]->area() ;
	double maxGelArea = aggregatesArea*degree;///nzones ;
	double maxGelRadius = sqrt(maxGelArea/(nzones*M_PI)) ;
	Function radius("t") ;
	radius *= (maxGelRadius / timeScale ) ;
//	radius = f_sqrt(radius) ;
	
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 0)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 2)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 3)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, 4)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, 5)) ;
	
	if(appliedLoadEta < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, appliedLoadEta));

	if(appliedLoadXi < 0)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition( SET_STRESS_XI, RIGHT_AFTER, appliedLoadXi));

	ViscoelasticityAndImposedDeformation * gel = dynamic_cast<ViscoelasticityAndImposedDeformation *>( (new ViscoElasticOnlyGelBehaviour())->getCopy() ) ;
	
	std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion*> > zones = PSDGenerator::get2DGrowingExpansiveZonesInAggregates( &F, aggregates, gel, radius, maxGelRadius, nzones*50, nzones) ;
	
	std::string name = "asr_creep_stress_450_degree_" ;
	name.append(argv[1]) ;

	std::string namebis = name ;
	namebis.append("_inter") ;

	std::ofstream summary ;
	summary.open(name.c_str(), std::ios::out) ;

	std::ofstream sinter ;
	sinter.open(namebis.c_str(), std::ios::out) ;

	std::vector<DelaunayTriangle *> trg ;
	std::vector<Point *> nodes ;
	int i = 0 ;
	while(F.getCurrentTime() < 450)
	{
		i++ ;
//		F.setDeltaTime(tau) ;
//		F.setMinDeltaTime(tau*1e-6) ;
		F.step() ;
		F.getAssembly()->setEpsilon(1e-8) ;
		if(trg.size() == 0)
		{
			trg = F.getElements2D() ;
			nodes = F.getNodes() ;
		}

/*		std::vector<Vector> inter = F.intermediateStates ;
		for(size_t j = 0 ; j < inter.size() ; j++)
		{
			sinter << inter[j][0] << "\t" ;
			double r = VirtualMachine().eval( radius, 0,0,0, inter[j][0] ) ;
			sinter << r << "\t" << nzones*r*r*M_PI/aggregatesArea << "\t" ;
			for(size_t k = 3 ; k < inter[j].size() ; k++)
				sinter << inter[j][k] << "\t" ;
			sinter << std::endl ;
			
		}
		F.intermediateStates.clear() ;*/

		summary << nodes[ nodes.size()-1]->t << "\t" ;
		double r = VirtualMachine().eval( radius, 0,0,0, nodes[ nodes.size()-1]->t ) ;
		summary << r << "\t" << nzones*r*r*M_PI/aggregatesArea << "\t" ;
		
		Vector strain = F.getAverageField(STRAIN_FIELD, -1, 1) ;
		Vector stress = F.getAverageField(REAL_STRESS_FIELD, -1, 1) ;
		summary << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" ;
		summary << stress[0] << "\t" << stress[1] << "\t" << stress[2] << "\t" ;
		summary << F.damageAreaInAggregates(trg) << "\t" << F.damageAreaInPaste(trg) << "\t" << F.averageDamage << "\t" << std::endl ;
		


		std::string nametrg = name ;
		nametrg.append("_trg_") ;
		nametrg.append(itoa(i)) ;
		TriangleWriter w( nametrg, &F, 1) ;
		w.getField( STRAIN_FIELD ) ;
		w.getField( REAL_STRESS_FIELD ) ;
		w.getField( TWFT_STIFFNESS ) ;
		w.getField( TWFT_DAMAGE ) ;
		w.write() ;
		
		tau += timestep ;
	}
	
	
	return 0 ;
}
