// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
//#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../utilities/writer/triangle_writer.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

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
//valeurs a tester pour le fitting
//    omp_set_num_threads(1) ;
//    std::valarray<Point *> pts(5) ;
//    pts[0] = new Point(-80, 50) ;
//    pts[1] = new Point(-50.00, 30) ;
//    pts[2] = new Point(0.00, 80) ;
//    pts[3] = new Point(0.00, 0.000) ;
//    pts[4] = new Point(-60, 0.000) ;
//    PolygonalSample inc(&samplers, pts) ;

//    std::vector<double> radii1 ;
//    radii1.push_back (15.0);
//    radii1.push_back (30.0);
//    EnrichmentRing inc1(&samplers, 22, 5, 50.0, -50.0) ;

//    Inclusion inc2(&samplers, 22.0, -50.0, -50) ;

//    Pore inc3(&samplers, 22.0, 50.0, 50) ;
///////////////////////////////////////////////////
// Discretisation temporelle///////////////////////
///////////////////////////////////////////////////
    int nbr_test = 6 ;
    //donnees sorties
    double poreArea = 0 ;
    double aggArea = 0 ;
    std::vector<double> vec_E_inclusions ;
    for (int k = 0; k < nbr_test; k++){
        double time_step = 0.05 ;
        std::stringstream ss;
        ss << k;
        std::string str = ss.str();
    ///////////////////////////////////////////////////
// Geometrie du maillage///////////////////////////
///////////////////////////////////////////////////
    Sample samplers(nullptr, 0.35,0.35,0.,0.) ;
    FeatureTree F(&samplers) ;
    F.setSamplingNumber(366) ;
    F.setOrder(LINEAR_TIME_LINEAR) ;
//    F.setMinDeltaTime(1e-9) /*Gere l'endommagement'*/;
    F.setSamplingRestriction( SAMPLE_RESTRICT_4 ) ;
        vec_E_inclusions.push_back ((14.0)*1.e9) ;
        F.setDeltaTime(time_step) ;
//    std::vector<Point> p ;
//    p.push_back( Point( 1e6/1e9, 1e6 ) ) ;
//    p.push_back( Point( 2e6/1e9, 0e6 ) ) ;
//    std::vector<Point> c ;
//    AsymmetricSpaceTimeNonLocalMultiLinearSofteningFractureCriterion crit( p, c, 1e9/*, 0.1*/ ) ;

///////////////////////////////////////////////////
// Modele et Materiaux      ///////////////////////
///////////////////////////////////////////////////
// Cement paste
    // Cement paste
    double k_elas = 23.0e9;
    double nu_elas = 0.24 ;
    Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
    std::vector<std::pair<Matrix, Matrix> > branches ;
    double factor_k = 1.28 ;
    std::vector<double> K_chaine_cp = {5.4e10/factor_k, 6.2e10/factor_k,3.4e10/factor_k} ;
	for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
	{
		double tau = std::pow(10., (double) i - 1 );
		Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas, true,  SPACE_TWO_DIMENSIONAL )  ;
		Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas, true,  SPACE_TWO_DIMENSIONAL ) ;
		branches.push_back(std::make_pair(K_i, Am_i)) ;
	}
//Aggregates
    double k_agg = 14.0e9;
    double nu_agg = 0.3;
    double rmax_agg = 0.02;
    double rmin_agg = rmax_agg*0.01;
    double f_agg = 0.4 + (k - 3)*0.015  ;
    // double f_agg = 0.39 ; 
//Voids
    double rmax_pore = 0.004;
    double rmin_pore = rmax_pore*0.01 ;
    double f_pore = 0.1;
    Matrix E_agg = Tensor::cauchyGreen( k_agg, nu_agg, true, SPACE_TWO_DIMENSIONAL) ;
    Viscoelasticity * paste = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
    Viscoelasticity * vaggregates = new Viscoelasticity( PURE_ELASTICITY, E_agg, 3) ;

    samplers.setBehaviour( paste );
    ConstantSizeDistribution * psd_b = new ConstantSizeDistribution() ;
    std::vector<Inclusion *> inc_pores = PSDGenerator::get2DInclusions(rmax_pore,f_pore*0.35*0.35, psd_b, PSDEndCriteria(rmin_pore, 0.01,5001)) ;
    std::vector<Inclusion *> inc_agg = PSDGenerator::get2DInclusions(rmax_agg, f_agg*0.35*0.35, psd_b, PSDEndCriteria(rmin_agg, 0.01,1001 )) ;
    std::vector<Feature *> inc_tot ;
    int i_agg = 0;
    std::cout << "nombre de pores " << inc_pores.size() << std::endl ;
    std::cout << "nombre de granulats " << inc_agg.size() << std::endl ;
    for(size_t i = 0 ; i < (inc_agg.size()); i++)
    {
        inc_tot.push_back(inc_agg[i]);
        inc_tot[i]->setBehaviour(vaggregates);
        aggArea += inc_tot[i]->area() ;
        i_agg += 1 ;
    }
    for(size_t i = 0 ; i < (inc_pores.size()); i++)
    {
        inc_tot.push_back(inc_pores[i]);
        inc_tot[i + i_agg]->setBehaviour(new VoidForm()) ;
        aggArea += inc_tot[i + i_agg]->area() ;
    }
    //Si on une distrib uniforme il faut les merger randomiquement
    //std::random_shuffle(pores.begin(), pores.end());2000000
    //pour la distance entre granulats prendre rmin/20
    inc_tot = placement2D(samplers.getPrimitive(), inc_tot, rmin_pore*0.25, 0, 20000) ;
    for(size_t i = 0 ; i < inc_tot.size() ; i++)
        F.addFeature(&samplers, inc_tot[i]);
    
    F.step();

    //std::vector<Feature *> finc = PSDGenerator::get2DConcrete( &F, vaggregates, 10000, 2e-1, 0.0000000001, psd_b, CIRCLE, 1.0, M_PI, 1.0e5, 0.4) ;
    //LogarithmicCreepWithExternalParameters aggregates("young_modulus = 60e9, poisson_ratio = 0.2") ;

///////////////////////////////////////////////////
// CONDITIONS AUX LIMITES      ////////////////////
// Prendre garde à bien mettre AFTER_ /////////////
///////////////////////////////////////////////////
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, 1.e6)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT_AFTER)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_AFTER)) ;
    BoundingBoxDefinedBoundaryCondition * load = new BoundingBoxDefinedBoundaryCondition( SET_STRESS_ETA, TOP_AFTER, -1.e6 ) ;

    F.setMaxIterationsPerStep(256) ;

    std::string test_data = "../../../../amie_data/data_visco_CASH_CSH1" + str ;
    std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::trunc);
    fichier_data << "frac agg ; " << f_agg << " frac pore ;  " << f_pore  << " nombre agg ;  " << inc_agg.size() << " nombre pores;  " << inc_pores.size() <<  std::endl ;
    fichier_data << "E agg ; " << k_agg << " Creep CSH ; " << K_chaine_cp[0] << ";" << K_chaine_cp[1] << ";" << K_chaine_cp[2] << std::endl ;
    fichier_data << "time ; " << "strain ; " << "stress " << std::endl ;
    fichier_data.close();

    //Chargement
    while(F.getCurrentTime() < 30)
    {
	time_step *= 1.1 ;
	F.setDeltaTime(time_step) ;

    F.step();

    std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::app);
    fichier_data << F.getCurrentTime() << " ; " << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << " ; " << F.getAverageField( REAL_STRESS_FIELD, -1, 1.)[1] << std::endl ;
    fichier_data.close();

    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;
    }
    //Visualisation des champs meca
    std::string test = "../../../../amie_data/test_visco_CASH_CSH1" + str;
    TriangleWriter writer(test.c_str(), &F, 1.) ;
    writer.getField(STRAIN_FIELD) ;
    writer.getField(PRINCIPAL_REAL_STRESS_FIELD) ;
//	writer.getField(SCALAR_DAMAGE_FIELD) ;
    writer.getField(TWFT_STIFFNESS) ;
    writer.write() ;
    //Dechargement
    /*
    F.addBoundaryCondition(load) ;
    time_step = 0.4 ;
    while(F.getCurrentTime() < 60)
    {
	i++ ;
    time_step *= 1.1 ;
	F.setDeltaTime(time_step) ;
//	load->setData(0.0001*F.getCurrentTime()) ;
    F.step();
//	load->setData(0.0001*F.getCurrentTime()) ;

    std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::app);
    fichier_data << F.getCurrentTime() << " ; " << F.getAverageField( STRAIN_FIELD, -1, 1.)[1] << " ; " << F.getAverageField( REAL_STRESS_FIELD, -1, 1.)[1] << std::endl ;
    fichier_data.close();
    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;
    }*/
    }

    return 0 ;

}

