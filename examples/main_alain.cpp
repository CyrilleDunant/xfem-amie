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
#include "../physics/fracturecriteria/maxstrain.h"
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
#include <sys/time.h>
#define DEBUG
using namespace Amie ;


int main(int argc, char *argv[])
{
	omp_set_num_threads(1) ;

	double duration = atof(argv[1]) ;
	int nsteps = atoi(argv[6]) ;
	bool elastic = (std::string(argv[2]) == std::string("elastic")) ;
	bool simple = (std::string(argv[3]) == std::string("simple")) ;
	double length = 0.02 ; if(simple) { length = 0.01 ; }
	int criterion = atoi(argv[4]) ;
	double deltad = atof(argv[5]) ;
	double timestep = duration/nsteps ;
	double increment = length*0.001/nsteps ;

	std::cout << timestep << std::endl ;

//	double timestep = rate ;//atof(argv[1]) ;
	double absoluteTimeTolerance = 1e-10;//0.1/(86400.*0.00001/0.0001) ;
	double relativeTimeTolerance = absoluteTimeTolerance*timestep ;

	Sample box(nullptr, length, 0.01,0.,0.) ;

	FeatureTree F(&box) ;
	if(simple)
		F.setSamplingNumber(1) ;
	else
		F.setSamplingNumber(6) ;

	ElasticOnlyPasteBehaviour dpaste(16e9, 0.3) ;
	SpaceTimeFiberBasedIsotropicLinearDamage * dampaste = new SpaceTimeFiberBasedIsotropicLinearDamage( deltad, absoluteTimeTolerance, 1.-1e-5) ;//atof(argv[1]) ) ;
	FractureCriterion * critpaste ;
	switch(criterion)
	{
		case 1:
			critpaste =  new SpaceTimeNonLocalLinearSofteningMaximumStrain( 0.00035, 0.00035*16e9, 0.0007) ;
			break ;
		case 2:
			critpaste =  new SpaceTimeNonLocalMaximumStrain( 0.00035, 0.00035*16e9) ;
			break ;
		case 3:
			critpaste =  new SpaceTimeNonLocalMaximumStress( 0.00035, 0.00035*16e9) ;
			break ;
	}
	ViscoelasticityAndFracture * paste ;
	if(elastic)
		paste = new ViscoelasticityAndFracture( PURE_ELASTICITY, dpaste.param, critpaste, dampaste ) ;
	else
		paste = new ViscoelasticityAndFracture( GENERALIZED_KELVIN_VOIGT, dpaste.param, dpaste.param*40./16., dpaste.param*0.005*40/16, critpaste, dampaste ) ;

	paste->getFractureCriterion()->setMaterialCharacteristicRadius(0.05) ;
	box.setBehaviour(paste);

	if(simple)
		F.setSamplingFactor(&box, 1) ;
	else
		F.setSamplingFactor(&box, 2) ;

	F.setOrder(LINEAR_TIME_LINEAR) ;
	F.setDeltaTime(timestep) ;
	if(simple)
		F.setMaxIterationsPerStep(32) ;
	else
		F.setMaxIterationsPerStep(256) ;
	F.setMinDeltaTime(absoluteTimeTolerance) ;

	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,0)) ;
	if(!elastic)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, LEFT, 0,2)) ;
//	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM_LEFT_AFTER, 0,1)) ;
	F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,1)) ;
	if(!elastic)
		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, BOTTOM, 0,3)) ;
	F.step() ;
	std::vector<DelaunayTriangle *> trg = F.getElements2D() ;
	int count = 0 ;
	Rectangle s( 0.002,0.02, 0., 0.) ;

	if(!simple)
	{
		for(size_t i = 0 ; i < trg.size() ; i++)
		{
			if(trg[i]->getBehaviour()->getFractureCriterion())
			{
				if(s.in(trg[i]->getCenter()))
				{
					SpaceTimeNonLocalLinearSofteningMaximumStrain* dam = dynamic_cast<SpaceTimeNonLocalLinearSofteningMaximumStrain *>(trg[i]->getBehaviour()->getFractureCriterion()) ;
					dam->upVal *= 0.99999 ;
					dam->maxstress *= 0.99999 ;
					count++ ;
				}
			}
		}
		std::cout << "paste weakened in " << count << " elements" << std::endl ;
	}


	std::fstream out ;
	std::string name = "tensile_" ;
	name.append(argv[1]) ;
	name.append("_") ;
	name.append(argv[2]) ;
	name.append("_") ;
	name.append(argv[3]) ;
	switch(criterion)
	{
		case 1:
			name.append("_softening") ;
			break ;
		case 2:
			name.append("_strain") ;
			break ;
		case 3:
			name.append("_stress") ;
			break ;
	}
	name.append("_") ;
	name.append(argv[5]) ;
	name.append("_") ;
	name.append(argv[6]) ;
	out.open(name.c_str(), std::ios::out) ;

	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition(SET_ALONG_INDEXED_AXIS, RIGHT_AFTER, 0.0, 0) ;
	F.addBoundaryCondition(disp) ;
	size_t i = 0 ;
	int done = 0 ;
	while(increment*i < length*0.00099)
	{
		i++;
		disp->setData(increment*i) ;
		F.setDeltaTime(timestep) ;
		F.setMinDeltaTime(absoluteTimeTolerance) ;
		while(!F.step()) 
		{ 
			std::string nameinter = name ;
			nameinter.append("_trg_inter") ;
			TriangleWriter writer(nameinter, &F, 1) ;
			writer.getField(STRAIN_FIELD) ;
			writer.getField(REAL_STRESS_FIELD) ;
			writer.getField(TWFT_CRITERION) ;
			writer.getField(TWFT_STIFFNESS) ;
			writer.getField(TWFT_DAMAGE) ;
			writer.write() ;

			double sigma = F.getAverageField(REAL_STRESS_FIELD, -1, 1.)[0] ;
			double epsilon = F.getAverageField(STRAIN_FIELD, -1, 1.)[0] ;
			double Ei = sigma/epsilon ;
			double E0 = 16e9 ;

			if((Ei/E0 < 0.95 && done == 0) || (Ei/E0 < 0.6 && done == 1))
			{
				std::fstream outi ;
				std::string namei = name ;
				namei.append("_") ;
				namei.append(itoa(done)) ;
				outi.open(namei.c_str(), std::ios::out) ;
				for(size_t j = 0 ; j < trg.size() ; j++)
				{
					if(std::abs(trg[j]->getCenter().getY()) < 0.0015)
						outi << trg[j]->getCenter().getX() << "\t" << trg[j]->getBehaviour()->getDamageModel()->getState().max() << std::endl ;
				}
				done++ ;
			}

//			out << increment*i << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, -1.)[0] << "\t" << F.getAverageField(STRAIN_FIELD, -1, -1.)[0] << "\t" << F.getCurrentTime() << "\t 0" << std::endl ;

			std::cout << increment*i << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, -1.)[0] << "\t" << F.getAverageField(STRAIN_FIELD, -1, -1.)[0] << std::endl ;
		}


			if(F.getAverageField(REAL_STRESS_FIELD, -1, 1.)[0] < 0.)
			{

				F.getAssembly()->print() ;
				Vector disp = F.getDisplacements() ;
				for(size_t i = 0 ; i < 32 ; i++)
					std::cout << disp[i] << "\t" ;

				exit(0) ;

				std::string nametrg = name ;
				nametrg.append("_trg_wrong") ;
				TriangleWriter writer(nametrg, &F, 1) ;
				writer.getField(STRAIN_FIELD) ;
				writer.getField(REAL_STRESS_FIELD) ;
				writer.getField(TWFT_STIFFNESS) ;
				writer.getField(TWFT_DAMAGE) ;
				writer.write() ;

			}


		out << increment*i << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1.)[0] << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1.)[0] << "\t" << F.getCurrentTime() << "\t 1" <<  std::endl ;

		std::cout << increment*i << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1.)[0] << "\t" << F.getAverageField(STRAIN_FIELD, -1, 1.)[0] << std::endl ;
	}

	std::fstream outi ;
	std::string namei = name ;
	namei.append("_") ;
	namei.append(itoa(done)) ;
	outi.open(namei.c_str(), std::ios::out) ;
	for(size_t j = 0 ; j < trg.size() ; j++)
	{
		if(std::abs(trg[j]->getCenter().getY()) < 0.0015)
			outi << trg[j]->getCenter().getX() << "\t" << trg[j]->getBehaviour()->getDamageModel()->getState().max() << std::endl ;
	}
	std::string nametrg = name ;
	nametrg.append("_trg_0") ;
	TriangleWriter writer(nametrg, &F, 1) ;
	writer.getField(STRAIN_FIELD) ;
	writer.getField(REAL_STRESS_FIELD) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.getField(TWFT_DAMAGE) ;
	writer.write() ;
		
	return 0 ;
}
