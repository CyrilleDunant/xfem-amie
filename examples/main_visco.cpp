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
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
#include "../physics/parallel_behaviour.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../features/pore.h"
#include "../utilities/writer/triangle_writer.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/material_laws/material_laws.h"
#include "../physics/material_laws/temperature_material_laws.h"
#include "../physics/orthotropicstiffness.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
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
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"

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
	Sample box(nullptr, 0.08,0.08,0.,0.) ;
	double rate = atof(argv[1]) ;

	FeatureTree F(&box) ;
	F.setSamplingNumber(24) ;
	F.setOrder(LINEAR_TIME_LINEAR) ;
	double time_step = 0.001 ;
	F.setDeltaTime(time_step) ;
	F.setMinDeltaTime(1e-9) ;
	F.setSamplingRestriction( SAMPLE_RESTRICT_4 ) ;

	LogarithmicCreepWithExternalParameters paste("young_modulus = 16e9, poisson_ratio = 0.3, creep_modulus = 30e9, creep_poisson = 0.3, creep_characteristic_time = 1.") ;
	LogarithmicCreepWithExternalParameters aggregates("young_modulus = 60e9, poisson_ratio = 0.2") ;

	Sample obox(nullptr, 0.08, 0.08, 0., 0.) ;
	FeatureTree G(&obox) ;
	G.setSamplingNumber(24) ;
	G.setOrder(LINEAR_TIME_LINEAR) ;
	G.setDeltaTime(time_step) ;
	G.setMinDeltaTime(1e-9) ;
	G.setSamplingRestriction( SAMPLE_RESTRICT_4 ) ;


	Matrix C = paste.C ;
	Matrix E = paste.E ;
	Matrix Cagg = aggregates.C ;

	std::vector<std::pair<Matrix, Matrix> > branches ;
	for(int i = -2 ; i < 6 ; i++)
	{
		double tau = std::pow(10., (double) i) ;
		std::cout << tau << "\t" << std::exp(-1./tau) << std::endl ;
		Matrix Ei = E/(log(10)*std::exp(-1./tau)) ;
		Matrix Zi = Ei*tau ;
		branches.push_back(std::make_pair(Ei, Zi)) ;
	}

	Viscoelasticity * vpaste = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, C, branches) ;
	Viscoelasticity * vaggregates = new Viscoelasticity( PURE_ELASTICITY, Cagg, 8) ;

	for(size_t i = 0 ; i < 8 ; i++)
	{
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, i*2) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, i*2+1) ) ;
	}
	for(size_t i = 0 ; i < 2 ; i++)
	{
		G.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, LEFT_AFTER, 0, i*2) ) ;
		G.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, BOTTOM_AFTER, 0, i*2+1) ) ;
	}

	box.setBehaviour(&vpaste) ;
	obox.setBehaviour(&paste) ;

	std::vector<Feature *> finc = PSDGenerator::get2DConcrete( &F, vaggregates, 50, 0.01, 0.00001, new PSDBolomeA()) ; 
	std::vector<Feature *> ginc = PSDGenerator::get2DConcrete( &G, aggregates, 50, 0.01, 0.00001, new PSDBolomeA()) ; 
	std::vector<Geometry *> fagg ;
	std::vector<Geometry *> gagg ;
	for(size_t i = 0 ; i < finc.size() ; i++)
	{
		fagg.push_back(dynamic_cast<Geometry *>(finc[i])) ;
		gagg.push_back(dynamic_cast<Geometry *>(ginc[i])) ;
	}

	F.step() ;
	G.step() ;

	BoundingBoxDefinedBoundaryCondition * strain = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, 0.004, 1) ;
	BoundingBoxDefinedBoundaryCondition * gstrain = new BoundingBoxDefinedBoundaryCondition( SET_ALONG_INDEXED_AXIS, TOP_AFTER, 0.004, 1) ;

	BoundingBoxDefinedBoundaryCondition * stress = new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6) ;
	BoundingBoxDefinedBoundaryCondition * gstress = new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 1e6) ;

//	F.addBoundaryCondition( stress ) ;
//	G.addBoundaryCondition( gstress ) ;

	F.addBoundaryCondition(strain) ;
	G.addBoundaryCondition(gstrain) ;

	int pasteFIndex = F.get2DMesh()->generateCache( F.getFeature(0)) ;
	int pasteGIndex = G.get2DMesh()->generateCache( G.getFeature(0)) ;
	int aggFIndex = F.get2DMesh()->generateCache( fagg ) ;
	int aggGIndex = G.get2DMesh()->generateCache( gagg ) ;

	std::string outfile = "visco_log_rate_" ;
	outfile.append(argv[1]) ;
	outfile.append(".txt") ;
	std::fstream out ;
	out.open(outfile.c_str(), std::ios::out) ;
	out << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl ;

	double target = 0.0004 ;
	while(F.getCurrentTime()*rate < target)
	{
//		if(time_step < 1.)
			time_step = target/(rate*100) ;
//		else
//			time_step *= 1.001 ;
		strain->setData( (F.getCurrentTime()+time_step*0.5)*rate ) ;
		gstrain->setData( (G.getCurrentTime()+time_step*0.5)*rate ) ;
	

		F.setDeltaTime( time_step ) ;
		F.step() ;
		G.setDeltaTime( time_step ) ;
		G.step() ;
		std::cout << F.getCurrentTime() << "\t" << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1] << "\t" << G.getAverageField( STRAIN_FIELD, -1, 1.)[1] <<  "\t" << G.getAverageField(REAL_STRESS_FIELD, -1, 1)[1] << "\t" << F.get2DMesh()->getField( REAL_STRESS_FIELD, pasteFIndex, -1., 1.)[1] << "\t" << G.get2DMesh()->getField( REAL_STRESS_FIELD, pasteGIndex, -1., 1.)[1]  << std::endl ;
		out << F.getCurrentTime() << "\t" << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << "\t" << F.getAverageField(REAL_STRESS_FIELD, -1, 1)[1] << "\t" << G.getAverageField( STRAIN_FIELD, -1, 1.)[1] <<  "\t" << G.getAverageField(REAL_STRESS_FIELD, -1, 1)[1] << "\t" << F.get2DMesh()->getField( REAL_STRESS_FIELD, pasteFIndex, -1., 1.)[1] << "\t" << G.get2DMesh()->getField( REAL_STRESS_FIELD, pasteGIndex, -1., 1.)[1]  << std::endl ;
	}


	return 0 ;
}

