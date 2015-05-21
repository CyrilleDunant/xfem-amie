
//
// C++ Implementation: mazars
//
// Description:
//
//
// Author: Adrien Hilaire <adrien.hilaire@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "mazars.h"
#include "../../mesher/delaunay.h"
#include "../damagemodels/damagemodel.h"
namespace Amie {

NonLocalMazars::NonLocalMazars(double thresh, double E, double nu, double Gf, double cstress, double cstrain, double radius, planeType pt, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
    , threshold(std::abs(thresh)), E(E), Gf(Gf), nu(nu), cstrain(cstrain), cstress(cstress), pt(pt)
{
    setMaterialCharacteristicRadius(radius);
    met = false ;
    B_t = (0.1*E*thresh) / (Gf - E*0.1*thresh*thresh*0.5);
    B_c = -1./(sqrt(2.)*nu*cstrain);
    A_c = -(E*thresh + sqrt(2.)*cstress*nu)/((E*std::exp(B_c*thresh - 1.)/B_c) - E*thresh);
}


NonLocalMazars::~NonLocalMazars()
{
}

double NonLocalMazars::gradeAtTime(ElementState &s, double t)
{
    std::pair<Vector, Vector> sstrain = getSmoothedFields(PRINCIPAL_REAL_STRESS_FIELD, PRINCIPAL_MECHANICAL_STRAIN_FIELD, s,  t) ;
    Vector stress = sstrain.first ;
    Vector strain = sstrain.second ;
    double strainzz =  -nu*(strain[0] + strain[1])/(1. - nu) ;
    double dama_predict = s.getParent()->getBehaviour()->getDamageModel()->getState().max();
    std::vector<double> posstrain ;
    double maxStrain = 0 ;
    double pseudo_dama = 0 ;
    double true_threshold = 0 ;
    double talpha = 0 ;
    double calpha = 0 ;
    //introduction of gamma to optimize the bi-compressive (tri-compressive) behavior
    double gamma = 1.0 ;
    double Trcsig = 0.0;
    double Trtsig = 0.0;
    met = false ;
    if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL && pt==PLANE_STRESS)
    {
        posstrain.push_back( (0.5*( std::abs(stress[0])  + stress[0] ) - nu*(0.5*( std::abs(stress[1])  + stress[1]) ))/(E*(1 - dama_predict))) ;
        posstrain.push_back( (0.5*( std::abs(stress[1])  + stress[1] ) - nu*(0.5*( std::abs(stress[0])  + stress[0])))/(E*(1 - dama_predict))) ;
        posstrain.push_back( ( std::max(0.0, stress[0]) + std::max(0.0,stress[1]) )*-nu/(E*(1 - dama_predict)) ) ;
        Trtsig = 0.5*( std::abs(stress[0])  + stress[0] ) + 0.5*( std::abs(stress[1])  + stress[1]) ;
        Trcsig =  ( std::min(0.0, stress[0]) + std::min(0.0,stress[1]));
        if( Trcsig <= -POINT_TOLERANCE &&  std::abs(Trtsig/Trcsig) <= 1e-6)
        {
            gamma =  std::sqrt( std::min(0.0, stress[0])*std::min(0.0, stress[0]) + std::min(0.0,stress[1])*std::min(0.0,stress[1]) ) / (-Trcsig) ;
        }
        maxStrain = gamma*std::max(0.0,sqrt ( pow( (0.5*( std::abs(strain[0])  + strain[0] )) ,2.0) + pow((0.5*( std::abs(strain[1])  + strain[1] ) ),2.0) +   pow((0.5*( std::abs(strainzz)  + strainzz ) ),2.0))) ;
        talpha = (posstrain[0]*(0.5*( std::abs(strain[0])  + strain[0] )) + posstrain[1]*(0.5*( std::abs(strain[1])  + strain[1] ))  + posstrain[2]*(0.5*( std::abs(strainzz)  + strainzz ))) / (maxStrain*maxStrain);
         talpha = 1.0;
	calpha = 1.0 - talpha ;
    }
    
    else if( s.getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL && pt==PLANE_STRAIN)
    {
        std::cout << "\n" << "MODEL MISSING FOR PLANE STRAIN CONDITIONS" << std::endl;
    }

    else if( s.getParent()->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
    {
        posstrain.push_back( (0.5*( std::abs(stress[0])  + stress[0] ) - nu*(0.5*( std::abs(stress[1])  + stress[1])  +  0.5*( std::abs(stress[2])  + stress[2])))/(E*(1 - dama_predict))) ;
        posstrain.push_back( (0.5*( std::abs(stress[1])  + stress[1] ) - nu*(0.5*( std::abs(stress[0])  + stress[0]) +  0.5*( std::abs(stress[2])  + stress[2])))/(E*(1 - dama_predict))) ;
        posstrain.push_back( (0.5*( std::abs(stress[2])  + stress[2] ) - nu*(0.5*( std::abs(stress[0])  + stress[0]) +  0.5*( std::abs(stress[1])  + stress[1])))/(E*(1 - dama_predict))) ;
        Trtsig = 0.5*( std::abs(stress[0])  + stress[0] ) + 0.5*( std::abs(stress[1])  + stress[1]) + 0.5*( std::abs(stress[2])  + stress[2]);
        Trcsig =  ( std::min(0.0, stress[0]) + std::min(0.0,stress[1]) + std::min(0.0,stress[2]));
        if( Trcsig <= -POINT_TOLERANCE &&  std::abs(Trtsig/Trcsig) <= 1e-4)
        {
            gamma =  std::sqrt( std::min(0.0, stress[0])*std::min(0.0, stress[0]) + std::min(0.0,stress[1])*std::min(0.0,stress[1]) + std::min(0.0,stress[2])*std::min(0.0,stress[2])) / (-Trcsig) ;
        }
        maxStrain = gamma*std::max(0. , std::sqrt( pow((0.5*(std::abs(strain[0])  + strain[0])),2.0) + pow((0.5*(std::abs(strain[1])  + strain[1])),2.0) + pow((0.5*(std::abs(strain[2])  + strain[2])),2.0) ));
        talpha = (posstrain[0]*(0.5*( std::abs(strain[0])  + strain[0] )) + posstrain[1]*(0.5*( std::abs(strain[1])  + strain[1] ))  + posstrain[2]*(0.5*( std::abs(strain[2])  + strain[2]  ))) / (maxStrain*maxStrain);
        if(talpha < 1.0e-3) {
            talpha= 0.;
        }
        calpha = 1.0 - talpha ;
    }

    //Critere seuil sur l'endommagement
    true_threshold = std::max(threshold, maxStrain);
    pseudo_dama =  std::min(talpha*(1. - (threshold/true_threshold) * (std::exp(- B_t*(true_threshold - threshold)))) + calpha*(1. - (threshold*(1 - A_c)/true_threshold)  -  A_c*(std::exp(- B_c*(true_threshold - threshold)))), 1.0 - POINT_TOLERANCE) ;
    if(pseudo_dama < 1.0e-6 || std::isnan(pseudo_dama)) {
        pseudo_dama = -1.0;
    }
    if( (maxStrain >= threshold) &&  (pseudo_dama >= dama_predict))
    {
        met = true ;
        double un_dama = pseudo_dama - dama_predict;
        return un_dama ;
    }
    return pseudo_dama - dama_predict;

}

double NonLocalMazars::grade( ElementState &s )
{
    return gradeAtTime(s, 0.0) ;
}

FractureCriterion * NonLocalMazars::getCopy() const
{
    return new NonLocalMazars(threshold,E, nu , Gf, cstress, cstrain,  getMaterialCharacteristicRadius(), pt) ;
}


NonLocalSpaceTimeMazars::NonLocalSpaceTimeMazars(double thresh, double E, double nu, double Gf, double cstress, double cstrain, double radius, planeType pt, MirrorState mirroring, double delta_x, double delta_y, double delta_z) : NonLocalMazars(thresh,  E,  nu,  Gf,  cstress,  cstrain,   radius, pt,  mirroring,  delta_x,  delta_y,  delta_z) { }

NonLocalSpaceTimeMazars::~NonLocalSpaceTimeMazars() { }

double NonLocalSpaceTimeMazars::grade(ElementState &s)
{
    double gradeBefore = gradeAtTime(s, -1) ;
    double gradeAfter = gradeAtTime(s, 1) ;

    if(gradeAfter < 0)
        return gradeAfter ;
    if(gradeBefore > 0)
    {
//         std::cout << gradeBefore << " vs " << gradeAfter << std::endl ;
//         exit(0) ;
//         return .1 ;
        return .999 ;
    }

    double upTime = 1 ;
    double downTime = -1 ;
    double testTime = 0 ;
    
    while(std::abs(upTime-downTime) > 1e-6)
    {
        double gradeTest = gradeAtTime(s, testTime) ;
        if(gradeTest < 0)
            downTime = testTime ;
        else if(gradeTest > 0)
            upTime = testTime ;
        else
            return testTime ;
        
        testTime = 0.5*(downTime+upTime) ;
    }
    return 1.-(testTime*.5+.5) ;
}

FractureCriterion *NonLocalSpaceTimeMazars::getCopy() const
{
    return new NonLocalSpaceTimeMazars( threshold,E, nu , Gf, cstress, cstrain,  getMaterialCharacteristicRadius(), pt ) ;
}




}
