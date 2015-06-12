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
#include "../physics/materials/paste_behaviour.h"
#include "../features/sample.h"

#include <fstream>
#include <omp.h>
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

	bool check = false ;
	bool spaceTime = true ;
	bool fiber = false ;
	bool newton = false ;
	bool noOpenMP = false ;
	double yieldstrain = 0.0005 ;
	double maxstrain = 0.0001 ;
	double young = 10e9 ;
	double radius = 0.01 ;
	int sampling = 32 ;

	if(argc > 1)
	{
		int i = 1 ;
		while(i < argc)
		{
			if(std::string(argv[i]) == std::string("--check"))
			{
				noOpenMP = true ;
				check = true ;
			}
			if(std::string(argv[i]) == std::string("--no-openmp"))
				noOpenMP = true ;
			if(std::string(argv[i]) == std::string("--sampling"))
			{
				i++ ;
				sampling = atoi(std::string(argv[i])) ;
			}
			if(std::string(argv[i]) == std::string("--radius"))
			{
				i++ ;
				radius = atof(std::string(argv[i])) ;
			}
			if(std::string(argv[i]) == std::string("--fiber"))
			{
				fiber = true ;
				spaceTime = true ;
				newton = false ;
			}
			if(std::string(argv[i]) == std::string("--newton"))
			{
				newton = true ;
				spaceTime = false ;
				fiber = false ;
			}
			if(std::string(argv[i]) == std::string("--space-time"))
			{
				spaceTime = true ;
				newton = false ;
			}
			if(std::string(argv[i]) == std::string("--yield-strain"))
			{
				i++ ;
				yieldstrain = atof(std::string(argv[i])) ;
			}
			if(std::string(argv[i]) == std::string("--max-strain"))
			{
				i++ ;
				maxstrain = atof(std::string(argv[i])) ;
			}
			if(std::string(argv[i]) == std::string("--young-modulus"))
			{
				i++ ;
				young = atof(std::string(argv[i])) ;
			}
			i++ ;
		}
	}

	if(noOpenMP)
		omp_set_num_threads(1) ;

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

	Matrix stiffness = ElasticOnlyPasteBehaviour(young).param ;

	if(spaceTime)
	{
		SpaceTimeNonLocalLinearSofteningMaximumStrain * fracST = new SpaceTimeNonLocalLinearSofteningMaximumStrain( maxstrain, maxtrain*young, yielstrain) ;
		fracST->setMaterialCharacteristicRadius( check ? 0.1 : radius ) ;
		SpaceTimeIsotropicLinearDamage * damST = new SpaceTimeIsotropicLinearDamage( 1. ) ;
		SpaceTimeFiberBasedIsotropicLinearDamage * damFST = new SpaceTimeFiberBasedIsotropicLinearDamage( 0.2, 0.01, 1. ) ;

		ViscoelasticityAndFracture * pasteST = fiber ? new ViscoelasticityAndFracture( PURE_ELASTICITY, stiffness, fracST, damFST ) : new ViscoelasticityAndFracture( PURE_ELASTICITY, stiffness, fracST, damST ) ;
	
		sample.setBehaviour( pasteST ) ;
	}
	else
	{
		// construct newton-raphson setup here
	}

	FeatureTree F( &sample ) ;
	F.setSamplingNumber( check ? 1 : sampling) ;
	F.setMaxIterationsPerStep(2000) ;
	F.setMinDeltaTime(1e-9) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, spaceTime ? LEFT_AFTER : LEFT) ) ;
	F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, spaceTime ? BOTTOM_LEFT_AFTER : BOTTOM_LEFT ) ) ;

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

	for(size_t i = 0 ; i < 100 ; i++)
	{
		disp->setData( i*(check ? 0.0000005 : 0.000001) ) ;
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
	}

	gettimeofday ( &time1, nullptr );
	double delta = time1.tv_sec * 1000000 - time0.tv_sec * 1000000 + time1.tv_usec - time0.tv_usec ;
	std::cout << "problem solved in " << delta/1000000 << " seconds" << std::endl ;

	return 0 ;
}

