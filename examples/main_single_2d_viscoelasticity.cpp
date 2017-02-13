// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
//#include "../physics/generalized_spacetime_viscoelasticity.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"
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
#include "../features/polygonSample.h"
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
//    omp_set_num_threads(1) ;
///////////////////////////////////////////////////
// Geometrie du maillage///////////////////////////
///////////////////////////////////////////////////
    Sample samplers(nullptr, 200,200,0.,0.) ;
    FeatureTree F(&samplers) ;
    F.setSamplingNumber(atof(argv[1])) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
//    F.setMinDeltaTime(1e-9) /*Gere l'endommagement'*/;
    F.setSamplingRestriction( 1. ) ;

    std::valarray<Point *> pts(5) ;
    pts[0] = new Point(-80, 50) ;
    pts[1] = new Point(-50.00, 30) ;
    pts[2] = new Point(0.00, 80) ;
    pts[3] = new Point(0.00, 0.000) ;
    pts[4] = new Point(-60, 0.000) ;
//    PolygonalSample inc(&samplers, pts) ;

    std::vector<double> radii1 ;
    radii1.push_back (15.0);
    radii1.push_back (30.0);
//    EnrichmentRing inc1(&samplers, 22, 5, 50.0, -50.0) ;

//    Inclusion inc2(&samplers, 22.0, -50.0, -50) ;

//    Pore inc3(&samplers, 22.0, 50.0, 50) ;
///////////////////////////////////////////////////
// Discretisation temporelle///////////////////////
///////////////////////////////////////////////////
    double time_step = 0.01 ;
    F.setDeltaTime(time_step) ;
//    std::vector<Point> p ;
//    p.push_back( Point( 1e6/1e9, 1e6 ) ) ;
//    p.push_back( Point( 2e6/1e9, 0e6 ) ) ;
//    std::vector<Point> c ;
//    AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion crit( p, c, 1e9/*, 0.1*/ ) ;

///////////////////////////////////////////////////
// Modele et Materiaux      ///////////////////////
///////////////////////////////////////////////////
    double E_rec = 1e10;
    double k_kv = 1.e20;
    double am_kv = 5.*k_kv;
    double k_mx = 2e30;
    double am_mx = 10.0*k_mx ;
    double nu_rec = 0.2 ;
    double delay = 0. ;
    int b = 0 ;
    int a = 0 ;
    Matrix E_kv = Tensor::cauchyGreen( k_kv, nu_rec, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    E_kv.print();
    Matrix C_kv = Tensor::cauchyGreen( am_kv, nu_rec,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    Matrix C_tau = Tensor::cauchyGreen( am_mx, nu_rec,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    Matrix C_eta = Tensor::cauchyGreen( k_mx, nu_rec,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    //TimeUnderLoadLogCreepAccumulator * accu = new TimeUnderLoadLogCreepAccumulator();
    Viscoelasticity * paste = new Viscoelasticity(BURGER, E_kv, C_kv, C_eta, C_tau,1);
    //IterativeMaxwell paste (E, t_char) ;
//    , &crit, new SpaceTimeFiberBasedIsotropicLinearDamage(0.9999, 1e-9, 0.99)) ;
    samplers.setBehaviour( paste );

    //LogarithmicCreepWithExternalParameters aggregates("young_modulus = 60e9, poisson_ratio = 0.2") ;

    F.step() ;

///////////////////////////////////////////////////
// CONDITIONS AUX LIMITES      ////////////////////
// Prendre garde à bien mettre AFTER_ /////////////
///////////////////////////////////////////////////
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_NORMAL_STRESS, TOP_AFTER, 10.0)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT_AFTER)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_AFTER)) ;
    BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, 5e6 ) ;

    F.setMaxIterationsPerStep(256) ;
	std::string test_data = "../../../../amie_data/data_viscoelasticity" ;
    std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::trunc);
    fichier_data << "time ; " << "strain ; " << "stress " << std::endl ;
    fichier_data.close();
/*
    std::string toto = "test_damage";
    std::fstream out ;
    out.open(toto.c_str(), std::ios::out) ; */
    //Chargement
    int i = 0 ;
    while(F.getCurrentTime() < 1)
    {
	i++ ;
	time_step *= 1.1 ;
	F.setDeltaTime(time_step) ;
//	load->setData(0.0001*F.getCurrentTime()) ;
    F.step();
	std::cout << "Pas de temps "<< i << std::endl ;
	std::cout << "\n" << F.getAverageField( TOTAL_STRAIN_FIELD)[1] << "\t" << F.getAverageField( REAL_STRESS_FIELD )[1] << std::endl;
    std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::app);
    if(fichier_data.is_open())
    {
        fichier_data << F.getCurrentTime() << " ; " << F.getAverageField( TOTAL_STRAIN_FIELD )[1] << " ; " << F.getAverageField( REAL_STRESS_FIELD )[1] << std::endl ;
        fichier_data.close();
    }
    else
    {
        std::cerr << "Erreur à l'ouverture !" << std::endl;
    }

	std::string test = "../../../../amie_data/test_single_2d_viscoelasticity" ;
	TriangleWriter writer(test.c_str(), &F, 1.) ;
	writer.getField(TOTAL_STRAIN_FIELD) ;
	writer.getField(PRINCIPAL_REAL_STRESS_FIELD) ;
//	writer.getField(SCALAR_DAMAGE_FIELD) ;
	writer.getField(TWFT_STIFFNESS) ;
	writer.write() ;
    }
    //Dechargement
    F.addBoundaryCondition(load) ;
    time_step = 0.01 ;
    while(F.getCurrentTime() < 150)
    {
        i++ ;

        time_step *= 1.1 ;
        F.setDeltaTime(time_step) ;
    //	load->setData(0.0001*F.getCurrentTime()) ;
        F.step();
        std::cout << "Pas de temps "<< i << std::endl ;
        std::cout << "\n" << F.getAverageField( TOTAL_STRAIN_FIELD )[1] << "\t" << F.getAverageField( REAL_STRESS_FIELD )[1] << std::endl;
        std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::app);
            if(fichier_data.is_open())
            {
                    fichier_data << F.getCurrentTime() << " ; " << F.getAverageField( TOTAL_STRAIN_FIELD )[1] << " ; " << F.getAverageField( REAL_STRESS_FIELD )[1] << std::endl ;
                    fichier_data.close();
            }
            else
            {
                    std::cerr << "Erreur à l'ouverture !" << std::endl;
            }

        std::string test = "../../../../amie_data/test_single_2d_viscoelasticity" ;
        TriangleWriter writer(test.c_str(), &F, 1.) ;
        writer.getField(TOTAL_STRAIN_FIELD) ;
        writer.getField(PRINCIPAL_REAL_STRESS_FIELD) ;
    //	writer.getField(SCALAR_DAMAGE_FIELD) ;
        writer.getField(TWFT_STIFFNESS) ;
        writer.write() ;
    }
    fichier_data.close();
    return 0 ;

}

