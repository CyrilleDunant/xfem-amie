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
#include "../features/microstructuregenerator.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../utilities/itoa.h"
#include "../physics/damagemodels/spacetimefiberbasedisotropiclineardamage.h"
#include "../physics/fracturecriteria/maxstrain.h"
#include "../geometry/geometry_base.h"
#include "../features/microstructuregenerator.h"

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
  int nbr_test = 10 ;
  for (int k = 0; k < (nbr_test-7); k++)
  {
        std::stringstream ss;
	  ss << k;
        std::string strfile = ss.str();
  ///////////////////////////////////////////////////
  // Discretisation temporelle///////////////////////
  ///////////////////////////////////////////////////
  //donnees sorties
  double poreArea = 0 ;
  double incArea = 0 ;
  double C3SArea = 0 ;
  double CHArea = 0 ;
  double AftArea = 0 ;
  double MsArea = 0 ;
  double cpArea = 0 ;
  double tol_param = 0.1e9;
  ///////////////////////////////////////////////////
  // Geometrie du maillage///////////////////////////
  ///////////////////////////////////////////////////
        RectangularFeature samplers(nullptr, 1.,1.,0.,0.) ;
  FeatureTree F(&samplers) ;
  F.setSamplingNumber(60) ;
  F.setOrder(LINEAR_TIME_LINEAR) ;
  F.setSamplingRestriction ( 0.01 ) ;
  ///////////////////////////////////////////////////
  // Modele et Materiaux      ///////////////////////
  ///////////////////////////////////////////////////
  // Elastique
  ///////////////////////////////////////////////////
  //Valeur cible du béton
  double kconcrete = 40.4e9 ;
  // CSH + Water
//   Value obtained with the scheme proposed by Roberts and Garboczi Elastic Properties of Model Porous Ceramics
//   overlapping sphere
//   double k_elas = 6.149e9 ;
//   double nu_elas = 0.162;
//   overlapping pore
//   double k_elas = 11.183e9 ;
//   double nu_elas = 0.21;  
//   overlapping ellipsoidal pores
  double k_elas = 8.2654e9 + 4*0.5e9 + 2*0.5e9;
  double nu_elas = 0.18;  
  Matrix E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //C3S+C2S+C3A+C4AF
  double E_C3S = 130.0e9;
  double nu_C3S = 0.3;
  double rmax_C3S = 0.04;
  double rmin_C3S = rmax_C3S*0.1;
  double f_C3S = 0.163;
  Matrix E_C3S_elas = Tensor::cauchyGreen( E_C3S, nu_C3S, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //CH
  double E_CH = 38.0e9;
  double nu_CH = 0.305;
  double rmax_CH = 0.02;
  double rmin_CH = rmax_CH*0.01;
  double f_CH = 0.151;
  Matrix E_CH_elas = Tensor::cauchyGreen( E_CH, nu_CH, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;
  //Aft    
  double E_Aft = 22.4e9;
  double nu_Aft = 0.255;
  double rmax_Aft = 0.02;
  double rmin_Aft = rmax_Aft*0.01;
  double f_Aft = 0.083;
  Matrix E_Aft_elas = Tensor::cauchyGreen( E_Aft, nu_Aft, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;    
  //Ms+Autre    
  double E_Ms = 42.3e9;
  double nu_Ms = 0.324;
  double rmax_Ms = 0.02;
  double rmin_Ms = rmax_Ms*0.01;
  double f_Ms = 0.044;
  Matrix E_Ms_elas = Tensor::cauchyGreen( E_Ms, nu_Ms, SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON) ;      
  double f_inc;
  f_inc = f_Aft+f_Ms+f_CH+f_C3S;
  // Viscoelastique
  //////////////////////////////////////////////////
  // CSH + Water
  std::vector<std::pair<Matrix, Matrix> > branches ;
  double factor_k = 0.38 + 4*(0.1)/10;
  std::vector<double> K_chaine_cp = {2.e11*factor_k,1.2e11*factor_k, 0.97e11*factor_k,7.4e10*factor_k, 2.6e10*factor_k,  2.2e10*factor_k} ;
  for(size_t i = 0 ; i < K_chaine_cp.size() ; i++)
  {
    double tau = 0.02*std::pow(10., (double) i - 1 );
    Matrix K_i = Tensor::cauchyGreen(K_chaine_cp[i], nu_elas,  SPACE_TWO_DIMENSIONAL , PLANE_STRESS, YOUNG_POISSON)  ;
    Matrix Am_i = Tensor::cauchyGreen( K_chaine_cp[i]*tau, nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
    branches.push_back(std::make_pair(K_i, Am_i)) ;
  }
  Viscoelasticity * paste = new Viscoelasticity( GENERALIZED_KELVIN_VOIGT, E_cp_elas, branches) ;
  Viscoelasticity * elaspaste = new Viscoelasticity( PURE_ELASTICITY, E_cp_elas,0,0) ; 
  //C3S+C2S+C3A+C4AF
  Viscoelasticity * vC3S = new Viscoelasticity( PURE_ELASTICITY, E_C3S_elas,3) ;
  //CH
  Viscoelasticity * vCH = new Viscoelasticity( PURE_ELASTICITY, E_CH_elas,3) ;  
  //Aft
  Viscoelasticity * vAft = new Viscoelasticity( PURE_ELASTICITY, E_Aft_elas,3) ;      
  //Ms+Autre
  Viscoelasticity * vMs = new Viscoelasticity( PURE_ELASTICITY, E_Ms_elas,3) ;          
  //Voids
//   double rmax_pore = 0.00625;
//   double rmin_pore = rmax_pore*0.1 ;
//   double f_pore = 0.1;
  //
    samplers.setBehaviour( elaspaste ); /*Elasticite  */
//       samplers.setBehaviour( paste ); /*Creep  */
  //
  //      ConstantSizeDistribution * psd_b = new ConstantSizeDistribution() ;
 GranuloFromCumulativePSD * psd_b = new GranuloFromCumulativePSD("../examples/data/irradiation/inclusion_psd", CUMULATIVE_PERCENT) ;
//GranuloFromFile * psd_b=new GranuloFromFile("../examples/data/irradiation/inclusion_ellips_psd",6);   
  //ConstantSizeDistribution * psd_b = new ConstantSizeDistribution();
 //GravelPolygonalInclusionGenerator * inc_elli = new GravelPolygonalInclusionGenerator(1.5, 0.2, 2, 5) ; 
  
  //     std::vector<Inclusion *> inc_pores = PSDGenerator::get2DInclusions(rmax_pore,f_pore*0.1*0.1, psd_b, PSDEndCriteria(rmin_pore, 0.01,20)) ;
  std::vector<Inclusion *> inc_C3S = PSDGenerator::get2DInclusions(rmax_C3S, f_C3S, psd_b, PSDEndCriteria(rmin_C3S, 1e-6,120000 )) ;
  std::vector<Inclusion *> inc_CH = PSDGenerator::get2DInclusions(rmax_CH, f_CH, psd_b, PSDEndCriteria(rmin_CH, 1e-6,120000 )) ;
  std::vector<Inclusion *> inc_Aft = PSDGenerator::get2DInclusions(rmax_Aft, f_Aft, psd_b, PSDEndCriteria(rmin_Aft, 1e-6,120000 )) ;
  std::vector<Inclusion *> inc_Ms = PSDGenerator::get2DInclusions(rmax_Ms, f_Ms, psd_b, PSDEndCriteria(rmin_Ms, 1e-6,120000 )) ;
  //     GeometryType geo_tr = TRIANGLE ;
  //     RandomDistribution * aire = new ConstantDistribution(1.0);
  //     RandomDistribution * ar = new ConstantDistribution(1.0);
  //     RandomDistribution * ori = new ConstantDistribution(0.0);
  //     InclusionConverter * inc_tri = new InclusionConverter(geo_tr, aire,  ar, ori ) ;
  //     std::vector<Feature *>  inc_ag = inc_tri->convert(inc_agg) ;
  //Circle
  //      EllipsoidalInclusionGenerator * inc_elli = new EllipsoidalInclusionGenerator(1.0) ; 
  //Ellipse
  //EllipsoidalInclusionGenerator * inc_elli = new EllipsoidalInclusionGenerator(1.0/2.0) ;
  //GravelPolygon
  
  std::vector<Feature *> inc_tot ;
  int i_C3S = 0;
  std::cout << "nombre de C3S " << inc_C3S.size() << std::endl ;
  std::cout << "nombre de CH " << inc_CH.size() << std::endl ;
  std::cout << "nombre de Aft " << inc_Aft.size() << std::endl ;
  std::cout << "nombre de Ms " << inc_Ms.size() << std::endl ;
  for(size_t i = 0 ; i < (inc_C3S.size()); i++)
  {
    inc_tot.push_back(inc_C3S[i]);
    //         inc_tot.push_back(inc_elli->convert(inc_agg[i]));
    inc_tot[i]->setBehaviour(vC3S);
    incArea += inc_tot[i]->area() ;
    i_C3S += 1 ;
  }
  int i_CH = 0;    
  for(size_t i = 0 ; i < (inc_CH.size()); i++)
  {
    inc_tot.push_back(inc_CH[i]);
    inc_tot[i + i_C3S]->setBehaviour(vCH) ;
    incArea += inc_tot[i + i_C3S]->area() ;
    i_CH += 1 ;
  }
  int i_Aft = 0;        
  for(size_t i = 0 ; i < (inc_Aft.size()); i++)
  {
    inc_tot.push_back(inc_Aft[i]);
    inc_tot[i + i_C3S+ i_CH]->setBehaviour(vAft) ;
    incArea += inc_tot[i + i_C3S+ i_CH]->area() ;
    i_Aft += 1 ;	
  }
  int i_Ms = 0;            
  for(size_t i = 0 ; i < (inc_Ms.size()); i++)
  {
    inc_tot.push_back(inc_Ms[i]);
    inc_tot[i + i_C3S + i_CH + i_Aft]->setBehaviour(vMs) ;
    incArea += inc_tot[i + i_C3S + i_CH + i_Aft]->area() ;
    i_Ms += 1 ;	
  }    
  //Si on une distrib uniforme il faut les merger randomiquement
  //std::random_shuffle(pores.begin(), pores.end());2000000
  //pour la distance entre granulats prendre rmin/20
  inc_tot = placement2D(samplers.getPrimitive(), inc_tot, rmin_C3S, 0, 60000) ;
  for(size_t i = 0 ; i < inc_tot.size() ; i++)
    F.addFeature(&samplers, inc_tot[i]);
  F.step();
  //std::vector<Feature *> finc = PSDGenerator::get2DConcrete( &F, vaggregates, 10000, 2e-1, 0.0000000001, psd_b, CIRCLE, 1.0, M_PI, 1.0e5, 0.4) ;
  //LogarithmicCreepWithExternalParameters aggregates("young_modulus = 60e9, poisson_ratio = 0.2") ;
  ///////////////////////////////////////////////////
  // CONDITIONS AUX LIMITES      ////////////////////
  // Prendre garde à bien mettre AFTER_ /////////////
  ///////////////////////////////////////////////////
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA, TOP_AFTER, -10.e6)) ;
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT_AFTER)) ;
  F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM_AFTER)) ;
  
  F.setMaxIterationsPerStep(256) ;
  
  
  
  //Chargement
  int i_kcp = 0 ;
  double kcp_inc = 1e9 ;
  double k_error = 1.;
  double k_error_0 = 1.;
  double k_homo ;
  double nu_homo ;
  
//   int cache_cp = F.get2DMesh()->generateCache((Feature *) &samplers) ;
  std::vector<bool> i_cp ;
  std::vector<bool> iC3S ;
  std::vector<bool> iCH ;
  std::vector<bool> iAft ;
  std::vector<bool> iMs ;
  double lame_C3S = E_C3S/(1+nu_C3S);
  double lame_CH = E_CH/(1+nu_CH);
  double lame_Aft = E_Aft/(1+nu_Aft);
  double lame_Ms = E_Ms/(1+nu_Ms);
  for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
  {      
    if((i->getBehaviour()->param[2][2] > (lame_C3S- tol_param)) && (i->getBehaviour()->param[2][2] < (lame_C3S + tol_param)))
    {
      C3SArea += i->area() ;
      iC3S.push_back(true) ;
      iCH.push_back(false);
      iAft.push_back(false);
      iMs.push_back(false);
      i_cp.push_back(false);
    }
    else if((i->getBehaviour()->param[2][2] > (lame_CH- tol_param)) && (i->getBehaviour()->param[2][2] < (lame_CH + tol_param)))
    {
      CHArea += i->area() ;
      iC3S.push_back(false) ;
      iCH.push_back(true);
      iAft.push_back(false);
      iMs.push_back(false);
      i_cp.push_back(false);
    }
    else if((i->getBehaviour()->param[2][2] > (lame_Aft- tol_param)) && (i->getBehaviour()->param[2][2] < (lame_Aft+ tol_param)))
    {
      AftArea += i->area() ;
      iC3S.push_back(false) ;
      iCH.push_back(false);
      iAft.push_back(true);
      iMs.push_back(false);
      i_cp.push_back(false);
    }
    else if((i->getBehaviour()->param[2][2] > (lame_Ms- tol_param)) && (i->getBehaviour()->param[2][2] < (lame_Ms+ tol_param)))
    {
      MsArea += i->area() ;
      iC3S.push_back(false) ;
      iCH.push_back(false);
      iAft.push_back(false);
      iMs.push_back(true);
      i_cp.push_back(false);
    }
    else
    {
      cpArea += i->area() ;
      iC3S.push_back(false) ;
      iCH.push_back(false);
      iAft.push_back(false);
      iMs.push_back(false);	
      i_cp.push_back(true);
    }
  }
  std::cout << " Area_C3s "<< " =============== "  << C3SArea  << std::endl ;
  std::cout << " Area_CH "<< " =============== "  << CHArea  << std::endl ;
  std::cout << " Area_Aft "<< " =============== "  << AftArea  << std::endl ;
  std::cout << " Area_Ms "<< " =============== "  << MsArea  << std::endl ;
  std::cout << " Area_matrix "<< " =============== "  << cpArea  << std::endl ; 
  
  std::string test_data = "homogenization_elas"  + strfile ;
  std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::trunc);
  fichier_data << "f_C3S ; " << C3SArea  << " ; nombre C3S ;  " << inc_C3S.size() 
               << "; f_CH ; " << CHArea  << " ; nombre CH ;  " << inc_CH.size() 
	       << " ; f_Aft ; " << AftArea<< " ; nombre Aft ;  " << inc_Aft.size() 
	       << " ; f_Ms ; " << MsArea<< " ; nombre Ms ;  " << inc_Ms.size() 
      << "; Kv_param ; " << K_chaine_cp[0] << " ; " <<  K_chaine_cp[1] << " ; "<<   K_chaine_cp[2] << " ; " <<   K_chaine_cp[3] << " ; " << std::endl ;  
  fichier_data<< "Time ;"  << "S11 ; " << "S22 ; "  << " S12 ;  " << " epsi11 ;  " << " epsi22 ;  " << " epsi12 ;  " << " dep2 ;  "  << " nbre points calcul dep2;  "  
       << " dep1 ;  "  << " nbre points calcul dep1;  "     << " dep3 ;  "  << " nbre points calcul dep3;  "  << "S11 C3S ; " << "S22 C3S ; "  << " S12 C3S;  "
       << "S11 Aft ; " << "S22 Aft ; "  << " S12 Aft;  " << "S11 Ms ; " << "S22 Ms ; "  << " S12 Ms;  "<< "S11 CSH ; " << "S22 CSH ; "  << " S12 CSH;  " 
       << " epsi11 C3S;  " << " epsi22 C3S;  " << " epsi12 C3S;" << "epsi11 Aft ; " << "epsi22 Aft ; "  << " epsi12 Aft;  "
       << "epsi11 Ms ; " << "epsi22 Ms ; "  << " epsi12 Ms;  "<< " epsi11 CSH;  " << " epsi22 CSH;  " << " epsi12 CSH;" <<  std::endl ;
  //     fichier_data << "frac agg ; " << f_agg << " ; frac pore ;  " << 1 - f_agg  << " ; nombre agg ;  " << inc_agg.size() << " ; surface agg effective ; " << aggArea  <<  std::endl ;
  //     fichier_data << "E agg ; " << E_agg << " ; nu_agg ; " << nu_agg  << std::endl ;
  //     fichier_data << "k_homo ; " << "nu_elas ; " << "k_elas " << std::endl ;
  
  double time_step = 0.01 ;    
  while(F.getCurrentTime() < 0.0126)
  {
    time_step *= 1.25 ;
    F.setDeltaTime(time_step) ;
    F.step();
    
    Vector strainC3S(3) ;
    Vector strainCH(3) ;  
    Vector strainAft(3) ;     
    Vector strainMs(3) ;
    Vector straincp(3) ;
    Vector stressC3S(3) ;
    Vector stressCH(3) ;  
    Vector stressAft(3) ;     
    Vector stressMs(3) ;
    Vector stresscp(3) ;
    Vector strainMaxC3S(3) ;
    Vector strainMaxCH(3) ;  
    Vector strainMaxAft(3) ;     
    Vector strainMaxMs(3) ;
    Vector strainMaxcp(3) ;
    Vector stressMaxC3S(3) ;
    Vector stressMaxCH(3) ;  
    Vector stressMaxAft(3) ;     
    Vector stressMaxMs(3) ;
    Vector stressMaxcp(3) ;   
    Vector tmp(3) ; 
    for(auto i = F.get2DMesh()->begin() ; i != F.get2DMesh()->end() ; i++)
    {
      if(iC3S[i.getPosition()])
      {
	double a = i->area() ;
	i->getState().getAverageField(TOTAL_STRAIN_FIELD, tmp, nullptr, 1.) ;
	strainC3S += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > strainMaxC3S[j])
	    strainMaxC3S[j] = tmp[j] ;
	}
	i->getState().getAverageField(REAL_STRESS_FIELD, tmp, nullptr, 1.) ;
	stressC3S += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > stressMaxC3S[j])
	    stressMaxC3S[j] = tmp[j] ;
	}
      }
      else if(iCH[i.getPosition()])
      {
	double a = i->area() ;
	i->getState().getAverageField(TOTAL_STRAIN_FIELD, tmp, nullptr, 1.) ;
	strainCH += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > strainMaxCH[j])
	    strainMaxCH[j] = tmp[j] ;
	}
	i->getState().getAverageField(REAL_STRESS_FIELD, tmp, nullptr, 1.) ;
	stressCH += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > stressMaxCH[j])
	    stressMaxCH[j] = tmp[j] ;
	}
      }
      else if(iAft[i.getPosition()])
      {
	double a = i->area() ;
	i->getState().getAverageField(TOTAL_STRAIN_FIELD, tmp, nullptr, 1.) ;
	strainAft += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > strainMaxAft[j])
	    strainMaxAft[j] = tmp[j] ;
	}
	i->getState().getAverageField(REAL_STRESS_FIELD, tmp, nullptr, 1.) ;
	stressAft += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > stressMaxAft[j])
	    stressMaxAft[j] = tmp[j] ;
	}
      }
      else if(iMs[i.getPosition()])
      {
	double a = i->area() ;
	i->getState().getAverageField(TOTAL_STRAIN_FIELD, tmp, nullptr, 1.) ;
	strainMs += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > strainMaxMs[j])
	    strainMaxMs[j] = tmp[j] ;
	}
	i->getState().getAverageField(REAL_STRESS_FIELD, tmp, nullptr, 1.) ;
	stressMs += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > stressMaxMs[j])
	    stressMaxMs[j] = tmp[j] ;
	}
      }
      else if(i_cp[i.getPosition()])
      {
	double a = i->area() ;
	i->getState().getAverageField(TOTAL_STRAIN_FIELD, tmp, nullptr, 1.) ;
	straincp += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > strainMaxcp[j])
	    strainMaxcp[j] = tmp[j] ;
	}
	i->getState().getAverageField(REAL_STRESS_FIELD, tmp, nullptr, 1.) ;
	stresscp += tmp*a ;
	for(size_t j = 0 ; j < 3 ; j++)
	{
	  if(tmp[j] > stressMaxcp[j])
	    stressMaxcp[j] = tmp[j] ;
	}
      }
    }
    
    
    strainC3S /= C3SArea ;
    stressC3S /= C3SArea ; 
    strainCH /= CHArea ;
    stressCH /= CHArea ; 
    strainAft /= AftArea ;
    stressAft /= AftArea ; 
    strainMs /= MsArea ;
    stressMs /= MsArea ; 	
    straincp /= cpArea ;
    stresscp /= cpArea ; 		
    
    int npoints = F.get2DMesh()->begin()->getBoundingPoints().size() ;

        double deltaE = 0;
	double deltanu = 0;
	double deltanu1 = 0;
        double deltacountE = 0;
	double deltacountnu = 0;
	double deltacountnu1 = 0;
        double volume = 0 ;
	

        for ( auto k = F.get2DMesh()->begin() ; k != F.get2DMesh()->end() ; k++ )
        {
            if ( k->getBehaviour()->type != VOID_BEHAVIOUR )
            {
                volume += k->area() ;
                for ( int p = npoints/2 ; p < npoints ; p++ )
                {

                    if (  (k->getBoundingPoint ( p ).getX()  > -.51) && ((k->getBoundingPoint ( p ).getY()) > 0.4985 ))
                    {
                        deltacountE++ ;
                        Vector dummy(k->getBehaviour()->getNumberOfDegreesOfFreedom()) ;
                        k->getState().getField(DISPLACEMENT_FIELD, k->getBoundingPoint ( p ), dummy, false); 
                        deltaE += dummy[1] ;
                    }
                    else if ( ( k->getBoundingPoint ( p ).getX()  > 0.4985) & ((k->getBoundingPoint ( p ).getY()) > -0.505 ))
                    {
                        deltacountnu++ ;
                        Vector dummy(k->getBehaviour()->getNumberOfDegreesOfFreedom()) ;
                        k->getState().getField(DISPLACEMENT_FIELD, k->getBoundingPoint ( p ), dummy, false); 
                        deltanu += dummy[0] ;
                    }
                    else if ( ( k->getBoundingPoint ( p ).getX()  < -0.4985) & ((k->getBoundingPoint ( p ).getY()) > -0.505 ))
                    {
                        deltacountnu1++ ;
                        Vector dummy(k->getBehaviour()->getNumberOfDegreesOfFreedom()) ;
                        k->getState().getField(DISPLACEMENT_FIELD, k->getBoundingPoint ( p ), dummy, false); 
                        deltanu1 += dummy[0] ;
                    }
                }
            }
        }


        deltaE /= deltacountE ;
	deltanu /= deltacountnu ;
	deltanu1 /= deltacountnu1 ;
    
//     k_homo = F.getAverageField(REAL_STRESS_FIELD)[1] / F.getAverageField( STRAIN_FIELD)[1];
//     nu_homo = -F.getAverageField(STRAIN_FIELD)[0] / F.getAverageField(STRAIN_FIELD)[1];
//     k_error = (k_homo - kconcrete )/kconcrete ;
/*    
    std::cout << " =============== " << stressC3S[0]  << " =============== " << stressC3S[1]  << " =============== " 	<< stressC3S[2]    << " =============== " 
    << " Strain " << strainC3S[0] << " =============== " <<  strainC3S[1]	  << " =============== " << strainC3S[2]  << std::endl ;  
    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;	
    
    std::cout << " =============== " << stressCH[0]  << " =============== " << stressCH[1]  << " =============== " 	<< stressCH[2]    << " =============== " 
    << " Strain " << strainCH[0] << " =============== " <<  strainCH[1]	  << " =============== " << strainCH[2]  << std::endl ;  
    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;
    
    std::cout << " =============== " << stressAft[0]  << " =============== " << stressAft[1]  << " =============== " 	<< stressAft[2]    << " =============== " 
    << " Strain " << strainAft[0] << " =============== " <<  strainAft[1]	  << " =============== " << strainAft[2]  << std::endl ;  
    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;
    
    std::cout << " =============== " << stressMs[0]  << " =============== " << stressMs[1]  << " =============== " 	<< stressMs[2]    << " =============== " 
    << " Strain " << strainMs[0] << " =============== " <<  strainMs[1]	  << " =============== " << strainMs[2]  << std::endl ;  
    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;	
    
    std::cout << " =============== " << stresscp[0]  << " =============== " << stresscp[1]  << " =============== " 	<< stresscp[2]    << " =============== " 
    << " Strain " << straincp[0] << " =============== " <<  straincp[1]	  << " =============== " << straincp[2]  << std::endl ;  
    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;		*/
    
    std::cout << " =============== " << F.getAverageField( REAL_STRESS_FIELD)[0]  << " =============== " << F.getAverageField( REAL_STRESS_FIELD)[1]  << " =============== " 
    << F.getAverageField( REAL_STRESS_FIELD)[2]    << " =============== " 
    << " Strain " << F.getAverageField( TOTAL_STRAIN_FIELD)[0]  << " =============== " << F.getAverageField( TOTAL_STRAIN_FIELD)[1]  
    << " =============== " << F.getAverageField( TOTAL_STRAIN_FIELD)[2]  << std::endl ;  
    std::cout << "temps " <<  F.getCurrentTime() << std::endl ;
    
    
    std::ofstream fichier_data(test_data.c_str(), std::ios::out | std::ios::app);
    fichier_data << F.getCurrentTime()<< " ; " <<F.getAverageField( REAL_STRESS_FIELD)[0]<< " ; " << F.getAverageField( REAL_STRESS_FIELD)[1]<< " ; "<< F.getAverageField( REAL_STRESS_FIELD)[2] 
    << " ; " << F.getAverageField( TOTAL_STRAIN_FIELD)[0] << " ; " << F.getAverageField( TOTAL_STRAIN_FIELD)[1] <<  " ; " << F.getAverageField( TOTAL_STRAIN_FIELD)[2] 
    << " ; " << deltaE <<  " ; " << deltacountE  << " ; " << deltanu <<  " ; " << deltacountnu <<  " ; " << deltanu1 <<  " ; " << deltacountnu1 <<  " ; "
    <<    stressC3S[0] <<  " ; " << stressC3S[1]	<<  " ; " << stressC3S[2]	<<  " ; " 
    <<    stressAft[0] <<  " ; " << stressAft[1]	<<  " ; " << stressAft[2]	<<  " ; "     
    <<    stressMs[0] <<  " ; " << stressMs[1]	<<  " ; " << stressMs[2]	<<  " ; "      
    << stresscp[0] << " ; "<<     stresscp[1] <<  " ; " << stresscp[2]	        <<  " ; " 
    << strainC3S[0]	<<  " ; " << strainC3S[1] << " ; "   <<    strainC3S[2] <<  " ; " 
    <<    strainAft[0] <<  " ; " << strainAft[1]	<<  " ; " << strainAft[2]	 << " ; "     
    <<    strainMs[0] <<  " ; " << strainMs[1]	<<  " ; " << strainMs[2]	<<  " ; "     
        << straincp[0]	<<  " ; " << straincp[1]	<<  " ; " << straincp[2] << " ; " 
    << std::endl ;
    fichier_data.close();
  }
  //      while( std::abs(k_error) > 0.1)
  //Visualisation des champs meca
  // 	if(k_error <= 0 && std::abs(k_error_0) >= std::abs(k_error) )
  // 	    {k_elas += kcp_inc;iii=0;kcp_inc *= 0.5 ;}
  // 	if(k_error >= 0 && std::abs(k_error_0) >= std::abs(k_error) )
  // 	    {k_elas -= kcp_inc;iii=1;kcp_inc *= 0.5 ;}
  // 	if(k_error <= 0 && std::abs(k_error_0) <= std::abs(k_error)) 
  // 	    {k_elas += kcp_inc;iii=2;}
  // 	if(k_error >= 0 && std::abs(k_error_0) <= std::abs(k_error)) 
  // 	    {k_elas -= kcp_inc;iii=3;}	
  // 	if(std::abs(F.getAverageField( STRAIN_FIELD)[1]) < 1e-6) 
  // 	    {k_elas = kconcrete; k_error = 1. ;}	    
  // 	k_error_0 = k_error;
  //         E_cp_elas = Tensor::cauchyGreen( k_elas, nu_elas,  SPACE_TWO_DIMENSIONAL, PLANE_STRESS, YOUNG_POISSON ) ;
  //             for (auto elem = F.get2DMesh()->begin(cacheID);
  //                     elem != F.get2DMesh()->end(cacheID) ; elem ++)
  //             {
  //                 elem->setBehaviour(F.get2DMesh(), new Stiffness(E_cp_elas));
  //             }
  
  //      }
  std::string test = "creep_homo_elas" + strfile ;
  TriangleWriter writer(test.c_str(), &F, 1.) ;
  writer.getField(TOTAL_STRAIN_FIELD) ;
  writer.getField(TWFT_STIFFNESS) ;	 
  writer.getField(REAL_STRESS_FIELD) ;
  writer.write() ;
  }
  
  return 0 ;
  
}


