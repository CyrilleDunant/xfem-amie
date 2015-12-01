// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../features/inclusion.h"
#include "../physics/viscoelasticity_and_fracture.h"
#include "../physics/damagemodels/spacetimeisotropiclineardamage.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../utilities/writer/triangle_writer.h"
#include "../utilities/parser.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/stiffness.h"
#include "../features/sample.h"

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

int main(int argc, char *argv[])
{
	timeval time0, time1 ;
	gettimeofday ( &time0, nullptr );

	CommandLineParser parser ;
	parser.addFlag("--check" ) ;
	parser.addFlag("--space-time" ) ;
	parser.addFlag("--fiber" ) ;
	parser.addFlag("--newton" ) ;
	parser.addFlag("--no-openmp" ) ;

	parser.addValue("--yield-strain", 0.0005 ) ;
	parser.addValue("--max-strain", 0.0001 ) ;
	parser.addValue("--young", 10e9 ) ;
	parser.addValue("--radius", 0.01 ) ;
	parser.addValue("--delta-damage", 0.2 ) ;
	parser.addValue("--sampling", 32 ) ;

	parser.parseCommandLine( argc, argv ) ;
	parser.printStatus() ;

	bool check = parser.getFlag("--check") ;
	bool spaceTime = parser.getFlag("--space-time") ;
	bool fiber = parser.getFlag("--fiber") ;
	bool newton = parser.getFlag("--newton") ;
	double yieldstrain = parser.getValue("--yield-strain") ;
	double maxstrain = parser.getValue("--max-strain") ;
	double young = parser.getValue("--young") ;
	double radius = parser.getValue("--radius") ;
	double deltaDamage = parser.getValue("--delta-damage") ;
	int sampling = parser.getValue("--sampling") ;

	if(newton)
	{
		spaceTime = false ;
		fiber = false ;
	}

	#ifdef HAVE_OMP
		omp_set_num_threads(1) ;
	#endif

	std::fstream out ;
	std::string index = "_" ;
	if(check)
		index.append("check_") ;
	if(spaceTime)
	{
		index.append("st") ;
		if(fiber)
			index.append("f") ;
		else
			index.append("i") ;
	}
	if(newton)
		index.append("nr") ;

	std::string filename = "ld"+index  ;
	out.open(filename.c_str(), std::ios::out) ;


        Sample sample(nullptr, check ? 0.1 : 0.2,0.1,0,0) ;

	Matrix stiffness = Stiffness(young, 0.2).param ;

	if(spaceTime)
	{
		SpaceTimeNonLocalLinearSofteningMaximumStrain * fracST = new SpaceTimeNonLocalLinearSofteningMaximumStrain( maxstrain, maxstrain*young, yieldstrain) ;
		fracST->setMaterialCharacteristicRadius( check ? 0.1 : radius ) ;
		SpaceTimeIsotropicLinearDamage * damST = new SpaceTimeIsotropicLinearDamage( 1. ) ;
		SpaceTimeFiberBasedIsotropicLinearDamage * damFST = new SpaceTimeFiberBasedIsotropicLinearDamage( deltaDamage, 0.01, 1. ) ;

		ViscoelasticityAndFracture * pasteST = fiber ? new ViscoelasticityAndFracture( PURE_ELASTICITY, stiffness, fracST, damFST ) : new ViscoelasticityAndFracture( PURE_ELASTICITY, stiffness, fracST, damST ) ;
	
		sample.setBehaviour( pasteST ) ;
	}
	else
	{
		// construct newton-raphson setup here
	}

	FeatureTree F( &sample ) ;
	F.setSamplingNumber( check ? 1 : sampling) ;
	F.setMaxIterationsPerStep(20000) ;
	F.setMinDeltaTime(1e-9) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, spaceTime ? LEFT_AFTER : LEFT) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, spaceTime ? BOTTOM_AFTER : BOTTOM ) ) ;

	if(!check)
	{
		Sample * notch = new Sample( nullptr, 0.002, 0.04, 0., 0.05 ) ;
		notch->setBehaviour( new VoidForm() ) ;
		F.addFeature(&sample, notch ) ;
		F.setSamplingFactor( notch, 1.5 ) ;
	}


	F.step() ;
	F.getAssembly()->setRemoveZeroOnlyLines( false ) ;

	BoundingBoxDefinedBoundaryCondition * disp = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_XI, spaceTime ? RIGHT_AFTER : RIGHT, 0. ) ;
	F.addBoundaryCondition( disp ) ;

	for(size_t i = 0 ; i < 1000 ; i++)
	{
		disp->setData( i*(check ? 0.00000005 : 0.0000001) ) ;
		F.step() ;
		std::cout << F.getCurrentTime() << "\t" << F.getAverageField( STRAIN_FIELD, -1,1 )[0] << "\t" << F.getAverageField( REAL_STRESS_FIELD, -1,1 )[0] << "\t" << F.getAverageField( SCALAR_DAMAGE_FIELD, -1,1 )[0] << std::endl ;
		out << F.getCurrentTime() << "\t" << F.getAverageField( STRAIN_FIELD, -1,1 )[0] << "\t" << F.getAverageField( REAL_STRESS_FIELD, -1,1 )[0] << "\t" << F.getAverageField( SCALAR_DAMAGE_FIELD, -1,1 )[0] << std::endl ;


		std::string trgf = "trg"+index+"_"+itoa(i) ;

		TriangleWriter trg( trgf, &F, 1.) ;
		trg.getField( STRAIN_FIELD ) ;
		trg.getField( REAL_STRESS_FIELD ) ;
		trg.getField( SCALAR_DAMAGE_FIELD ) ;
		trg.getField( TWFT_STIFFNESS ) ;
		trg.write() ;

		if( i > 1 && F.getAverageField( REAL_STRESS_FIELD, -1,1 )[0] < 100)
			break ;


	}

	gettimeofday ( &time1, nullptr );
	double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
	std::cout << "problem solved in " << delta/1000000 << " seconds" << std::endl ;

	return 0 ;
}

