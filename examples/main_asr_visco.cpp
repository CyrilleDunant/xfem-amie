// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/features.h"
#include "../utilities/granulo.h"

#include <fstream>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <cmath>
#include <typeinfo>
#include <limits>
#include <sys/time.h>
#define DEBUG


using namespace Amie ;

RectangularFeature box(nullptr, 0.08, 0.08,0.0,0.0) ;

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
	double appliedLoadEta = atof(argv[2])*(-1e6) ;
	double appliedLoadXi = atof(argv[3])*(-1e6) ;
	double degree = atof(argv[1]) ;
	tau /= degree/0.001 ;
		
	int nzones = 600 ;
	int naggregates = 1500 ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(200) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(tau) ;
	F.setMinDeltaTime(tau*1e-9) ;
	F.setMaxIterationsPerStep(5000) ;
		
//	box.setBehaviour( new ViscoElasticOnlyPasteBehaviour() );
	box.setBehaviour( new PasteBehaviour(false, true) );
	std::vector<Feature *> feats = PSDGenerator::get2DConcrete( &F, new AggregateBehaviour(false, true), naggregates, 0.008, 0.0001, new PSDBolomeA()) ;
	
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

	ViscoelasticityAndImposedDeformation * gel = dynamic_cast<ViscoelasticityAndImposedDeformation *>( (new GelBehaviour(true))->getCopy() ) ;
	
	std::vector<std::pair<TimeDependentHomogenisingInclusion *, Inclusion*> > zones = PSDGenerator::get2DGrowingExpansiveZonesInAggregates( &F, aggregates, gel, radius, maxGelRadius, nzones*50, nzones) ;
	
	std::string name = "asr_creep_stress_450_degree_" ;
	name.append(argv[1]) ;

	std::string namebis = name ;
	namebis.append("_inter") ;

	std::ofstream summary ;
	summary.open(name.c_str(), std::ios::out) ;

	std::ofstream sinter ;
	sinter.open(namebis.c_str(), std::ios::out) ;

	std::vector<Point *> nodes ;
	int i = 0 ;
	while(F.getCurrentTime() < 450)
	{
		i++ ;
//		F.setDeltaTime(tau) ;
//		F.setMinDeltaTime(tau*1e-6) ;
		F.step() ;
		F.getAssembly()->setEpsilon(1e-8) ;
		if(nodes.size() == 0)
		{
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

		summary << nodes[ nodes.size()-1]->getT() << "\t" ;
		double r = VirtualMachine().eval( radius, 0,0,0, nodes[ nodes.size()-1]->getT() ) ;
		summary << r << "\t" << nzones*r*r*M_PI/aggregatesArea << "\t" ;
		
		Vector strain = F.getAverageField(TOTAL_STRAIN_FIELD, 1) ;
		Vector stress = F.getAverageField(REAL_STRESS_FIELD, 1) ;
		summary << strain[0] << "\t" << strain[1] << "\t" << strain[2] << "\t" ;
		summary << stress[0] << "\t" << stress[1] << "\t" << stress[2] << "\t" ;
		summary << F.damageAreaInAggregates(F.get2DMesh()->begin(), F.get2DMesh()->end()) << "\t" << F.damageAreaInPaste(F.get2DMesh()->begin(), F.get2DMesh()->end()) << "\t" << F.averageDamage << "\t" << std::endl ;
		


		std::string nametrg = name ;
		nametrg.append("_trg_") ;
		nametrg.append(itoa(i)) ;
		TriangleWriter w( nametrg, &F, 1) ;
		w.getField( TOTAL_STRAIN_FIELD ) ;
		w.getField( REAL_STRESS_FIELD ) ;
		w.getField( TWFT_STIFFNESS ) ;
		w.getField( TWFT_DAMAGE ) ;
		w.write() ;
		
		tau += timestep ;
	}
	
	
	return 0 ;
}

