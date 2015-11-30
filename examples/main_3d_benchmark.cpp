// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/laplacian.h"
#include "../features/pore.h"
#include "../features/sample3d.h"
#include "../features/inclusion3d.h"
#include "../features/expansiveZone3d.h"
#include "../utilities/granulo.h"
#include "../physics/stiffness.h"
#include <sys/time.h>

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include "../utilities/writer/voxel_writer.h"
#define DEBUG
using namespace Amie ;

int viewangle = 0 ;
int viewangle2 = 0 ;

FeatureTree * featureTree ;
std::vector<bool> cracked ;

double E_min = 10;
double E_max = 0;

double prop ;

double x_max = 0 ;
double y_max = 0 ;
double z_max = 0 ;

double x_min = 0 ;
double y_min = 0 ;
double z_min = 0 ;


timeval t0, t1 ;

double timepos = 0.1 ;

bool firstRun = true ;

std::vector<DelaunayTriangle *> tris__ ;

Vector b(0) ;
Vector x(0) ;
Vector sigma(0) ;
Vector sigma11(0) ;
Vector sigma22(0) ;
Vector sigma33(0) ;
Vector sigma12(0) ;
Vector sigma13(0) ;
Vector sigma23(0) ;
Vector epsilon(0) ;
Vector epsilon11(0) ;
Vector epsilon22(0) ;
Vector epsilon33(0) ;
Vector epsilon12(0) ;
Vector epsilon13(0) ;
Vector epsilon23(0) ;
Vector vonMises(0) ;
Vector stiffness(0) ;
Vector angle(0) ;
Vector damage(0) ;

Vector epsilon_bcx(6) ;
Vector sigma_bcx(6) ;
Vector epsilon_bcy(6) ;
Vector sigma_bcy(6) ;
Vector epsilon_bcz(6) ;
Vector sigma_bcz(6) ;

double nu = 0.2 ;
double E_agg = 100 ;//softest
double E_paste = 1 ;//stiff
double E_stiff = E_agg*10 ;//stiffer
double E_soft = E_agg/10; //stiffest

double factor = 0.3 ;
MinimumAngle cri(M_PI/6.) ;
bool nothingToAdd = false ;
bool dlist = false ;
int count = 0 ;
double aggregateArea = 0;

typedef enum
{
    XS1,
    S1,
    S2024,
    S3200,
    O1,
} BenchmarkMicrostructure ;

typedef enum
{
    DIFFUSION,
    ELASTICITY,
} BenchmarkPhenomenon ;

BenchmarkMicrostructure getMicrostructure(std::string micro)
{
    if(micro == std::string("XS1"))
        return XS1 ;
    if(micro == std::string("S1"))
        return S1 ;
    if(micro == std::string("S2024"))
        return S2024 ;
    if(micro == std::string("S3200"))
        return S3200 ;
    if(micro == std::string("O1"))
        return O1 ;
    return S1 ;
}

BenchmarkPhenomenon getPhenomenon(std::string pheno)
{
    if(pheno == std::string("diffusion"))
        return DIFFUSION ;
    if(pheno == std::string("elasticity"))
        return ELASTICITY ;
    return DIFFUSION ;
}

double getLength(BenchmarkMicrostructure micro)
{
    switch(micro)
    {
    case(S1) :
        return 0.15 ;
    case(S2024) :
        return 0.15 ;
    case(S3200) :
        return 400. ;
    case(O1) :
        return 1. ;
    default:
        return 1. ;
    }
           return 1. ;
}

int main(int argc, char *argv[])
{
    gettimeofday(&t0, NULL) ;
    std::cerr << "usage: benchmark <microstructure> <phenomenon> <order> <scale> <properties> <sampling>" << std::endl ;
    std::cerr << "all fields are required" << std::endl ;
    std::cerr << "<microstructure>\tstring among <XS1, S1, S2024, S3200, O1>" << std::endl ;
    std::cerr << "<phenomenon>\t\tstring among <diffusion, elasticity>" << std::endl ;
    std::cerr << "<order>\t\t\tinteger : order of the elements (1 for linear, 2 for quadratic)" << std::endl ;
    std::cerr << "<scale>\t\t\tdouble : scale factor" << std::endl ;
    std::cerr << "<properties>\t\tdouble : value of the Young's Modulus or Diffusion coefficient of the inclusion(s)" << std::endl ;
    std::cerr << "<sampling>\t\tinteger : number of points at the surface of the REV" << std::endl ;

//    if(argc != 7)
//        return 1 ;

    BenchmarkMicrostructure micro = S1 ;//getMicrostructure(std::string(argv[1])) ;
    BenchmarkPhenomenon pheno = ELASTICITY ;//getPhenomenon(std::string(argv[2])) ;
    int order = 1; //atof(argv[3]) ;
    if(order < 1 || order > 2) {
        order = 1 ;
    }
    double scale = 100.;//atof(argv[4]) ;
    prop = 3;//atof(argv[5]) ;
    int sampling = 8000;//atof(argv[6]) ;
    double length = getLength(micro) ;

    double size = scale*length ;
    double halfSize = size/2 ;

    Sample3D sample(nullptr, size, size, size, halfSize, halfSize, halfSize) ;
    FeatureTree F(&sample) ;
    featureTree = &F ;

    F.setProjectionOnBoundaries(false) ;

    Form* behaviour = nullptr ;
    Matrix m1(6,6) ;
    switch(pheno)
    {
        case DIFFUSION:
        {
            Matrix d0(3,3) ;
            double lambda = 1 ;
            d0[0][0] = lambda ;
            d0[1][1] = lambda ;
            d0[2][2] = lambda ;
            sample.setBehaviour(new Laplacian(d0)) ;

            Matrix d1(3,3) ;
            lambda = prop ;
            d1[0][0] = lambda ;
            d1[1][1] = lambda ;
            d1[2][2] = lambda ;
            behaviour = new Laplacian(d1) ;
            break ;
        }
        case ELASTICITY:
        {

            double nu = 0.2 ;
            double E = 1. ;
            sample.setBehaviour(new Stiffness(Tensor::cauchyGreen(E, nu, true,  SPACE_THREE_DIMENSIONAL, PLANE_STRESS ))) ;

//             std::cout << prop << std::endl ;

            Vector alpha(0., 6) ;
            behaviour = new StiffnessWithImposedDeformation(Tensor::cauchyGreen(10., .2, true,  SPACE_THREE_DIMENSIONAL, PLANE_STRESS ),alpha) ;
            break ;
        }
    }

    std::string str_micro = "S1" ;
//    if(micro != O1)
//    {
    std::vector<Feature * > inclusions ;
//        if(micro == S1)
//        {*/
    inclusions.push_back(new Inclusion3D(0.0623*scale, sample.getCenter().getX(), sample.getCenter().getY(), sample.getCenter().getZ())) ;
    Vector str(0., 6) ;
//     inclusions.push_back(new ExpansiveZone3D(&sample , 0.0623*scale, sample.getCenter().getX(), sample.getCenter().getY(), sample.getCenter().getZ(),Tensor::cauchyGreen(10., .2, true,  SPACE_THREE_DIMENSIONAL ), str )) ;
    inclusions[0]->setBehaviour(behaviour) ;
    F.addFeature(&sample, inclusions[0]) ;

    F.setSamplingNumber(sampling) ;
    F.setMaxIterationsPerStep(2);
    F.setDeltaTime(0.001);
    F.setElementGenerationMethod(0,true) ;

    if(order == 2)
        F.setOrder(QUADRATIC) ;
    else
        F.setOrder(LINEAR) ;

    Function pos("x") ;
    Function grad = pos*1./(length*scale) ;


    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_XI, RIGHT, 1.)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI, LEFT)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ZETA, BACK)) ;

    F.step() ;

    x.resize(F.getDisplacements().size()) ;
    x = featureTree->getDisplacements() ;

    switch(pheno)
    {
    case DIFFUSION:
    {
        Vector gradient11(4*F.get3DMesh()->begin().size()) ;
        Vector gradient22(4*F.get3DMesh()->begin().size()) ;
        Vector gradient33(4*F.get3DMesh()->begin().size()) ;
        Vector flux11(4*F.get3DMesh()->begin().size()) ;
        Vector flux22(4*F.get3DMesh()->begin().size()) ;
        Vector flux33(4*F.get3DMesh()->begin().size()) ;
        {
            Vector gradient(12*F.get3DMesh()->begin().size()) ;
            Vector flux(12*F.get3DMesh()->begin().size()) ;
            {
                std::pair<Vector,Vector> gradient_flux ;
                gradient_flux.first.resize(12*F.get3DMesh()->begin().size()) ;
                gradient_flux.second.resize(12*F.get3DMesh()->begin().size()) ;
                gradient_flux = featureTree->getGradientAndFlux(F.get3DMesh()->begin(),F.get3DMesh()->end()) ;
                gradient = gradient_flux.first ;
                flux = gradient_flux.second ;
            }

            std::cout << "get gradient and flux..." << std::endl ;
            for(size_t i = 0 ; i < F.get3DMesh()->begin().size() ; i++)
            {
                for(size_t j = 0 ; j < 4 ; j++)
                {
                    gradient11[4*i+j] = gradient[4*3*i+3*j+0] ;
                    gradient22[4*i+j] = gradient[4*3*i+3*j+1] ;
                    gradient33[4*i+j] = gradient[4*3*i+3*j+2] ;
                    flux11[4*i+j] = flux[4*3*i+3*j+0] ;
                    flux22[4*i+j] = flux[4*3*i+3*j+1] ;
                    flux33[4*i+j] = flux[4*3*i+3*j+2] ;
                }
            }
        }

        std::cout << "averaging gradient and flux..." << std::endl ;
        Vector average_gradient(3) ;
        Vector average_flux(3) ;
        double total_volume = 0. ;
        for(auto i = F.get3DMesh()->begin() ; i != F.get3DMesh()->end() ; i++)
        {
            double volume = i->volume() ;
            total_volume += volume ;
            for(size_t j = 0 ; j < 4 ; j++)
            {
                average_gradient[0] += gradient11[4*i.getPosition()+j]*volume/4 ;
                average_gradient[1] += gradient22[4*i.getPosition()+j]*volume/4 ;
                average_gradient[2] += gradient33[4*i.getPosition()+j]*volume/4 ;
                average_flux[0] += flux11[4*i.getPosition()+j]*volume/4 ;
                average_flux[1] += flux22[4*i.getPosition()+j]*volume/4 ;
                average_flux[2] += flux33[4*i.getPosition()+j]*volume/4 ;
            }
        }

        std::cout << std::endl ;
        std::cout << "max value :" << x_max << std::endl ;
        std::cout << "min value :" << x_min << std::endl ;
        std::cout << "max flux11 :" << flux11.max() << std::endl ;
        std::cout << "min flux11 :" << flux11.min() << std::endl ;
        std::cout << "max flux22 :" << flux22.max() << std::endl ;
        std::cout << "min flux22 :" << flux22.min() << std::endl ;
        std::cout << "max flux33 :" << flux33.max() << std::endl ;
        std::cout << "min flux33 :" << flux33.min() << std::endl ;

        std::cout << "max gradient11 :" << gradient11.max() << std::endl ;
        std::cout << "min gradient11 :" << gradient11.min() << std::endl ;
        std::cout << "max gradient22 :" << gradient22.max() << std::endl ;
        std::cout << "min gradient22 :" << gradient22.min() << std::endl ;
        std::cout << "max gradient33 :" << gradient33.max() << std::endl ;
        std::cout << "min gradient33 :" << gradient33.min() << std::endl ;

        std::cout << "average flux11 : " << average_flux[0]/total_volume << std::endl ;
        std::cout << "average flux22 : " << average_flux[1]/total_volume << std::endl ;
        std::cout << "average flux33 : " << average_flux[2]/total_volume << std::endl ;
        std::cout << "average gradient11 : " << average_gradient[0]/total_volume << std::endl ;
        std::cout << "average gradient22 : " << average_gradient[1]/total_volume << std::endl ;
        std::cout << "average gradient33 : " << average_gradient[2]/total_volume << std::endl ;


        std::string filebench("benchmark.txt") ;
        std::fstream out ;
        out.open(filebench.c_str(), std::ios::out|std::ios::app) ;
        out << "DIFFUSION\t" << str_micro << "\t" << "D_inc = " << prop << "\t" ;
        if(order==2)
            out << "QUAD" ;
        out << "dof = " << x.size() << "\t"
            << "D11 = " << average_flux[0]/average_gradient[0] << std::endl ;
        out.close() ;

        /*        VoxelWriter vw("fem_s1", 100) ;
                vw.getField(featureTree, VWFT_FLUX) ;
                vw.write();*/

        break ;
    }
    case ELASTICITY:
    {
        Vector strain11(4*F.get3DMesh()->begin().size()) ;
        Vector strain22(4*F.get3DMesh()->begin().size()) ;
        Vector strain33(4*F.get3DMesh()->begin().size()) ;
        Vector stress11(4*F.get3DMesh()->begin().size()) ;
        Vector stress22(4*F.get3DMesh()->begin().size()) ;
        Vector stress33(4*F.get3DMesh()->begin().size()) ;
        {
            Vector stress(24*F.get3DMesh()->begin().size()) ;
            Vector strain(24*F.get3DMesh()->begin().size()) ;
            {
                std::pair<Vector,Vector> stress_strain ;
                stress_strain.first.resize(24*F.get3DMesh()->begin().size()) ;
                stress_strain.second.resize(24*F.get3DMesh()->begin().size()) ;
                stress_strain = featureTree->getStressAndStrain(F.get3DMesh()->begin(),F.get3DMesh()->end()) ;
                stress = stress_strain.first ;
                strain = stress_strain.second ;
            }

            std::cout << "get stress and strain..." << std::endl ;
            for(auto i = F.get3DMesh()->begin() ; i != F.get3DMesh()->end() ; i++)
            {
                for(size_t j = 0 ; j < 4 ; j++)
                {
                    stress11[4*i.getPosition()+j] = stress[4*6*i.getPosition()+6*j+0] ;
                    stress22[4*i.getPosition()+j] = stress[4*6*i.getPosition()+6*j+1] ;
                    stress33[4*i.getPosition()+j] = stress[4*6*i.getPosition()+6*j+2] ;
                    strain11[4*i.getPosition()+j] = strain[4*6*i.getPosition()+6*j+0] ;
                    strain22[4*i.getPosition()+j] = strain[4*6*i.getPosition()+6*j+1] ;
                    strain33[4*i.getPosition()+j] = strain[4*6*i.getPosition()+6*j+2] ;
                }
            }
        }

        std::cout << "averaging stress and strain..." << std::endl ;
        Vector average_stress(3) ;
        Vector average_strain(3) ;
        double total_volume = 0. ;
        for(auto i = F.get3DMesh()->begin() ; i != F.get3DMesh()->end() ; i++)
        {
            double volume = i->volume() ;
            total_volume += volume ;
            for(size_t j = 0 ; j < 4 ; j++)
            {
                average_stress[0] += stress11[4*i.getPosition()+j]*volume/4 ;
                average_stress[1] += stress22[4*i.getPosition()+j]*volume/4 ;
                average_stress[2] += stress33[4*i.getPosition()+j]*volume/4 ;
                average_strain[0] += strain11[4*i.getPosition()+j]*volume/4 ;
                average_strain[1] += strain22[4*i.getPosition()+j]*volume/4 ;
                average_strain[2] += strain33[4*i.getPosition()+j]*volume/4 ;
            }
        }

        std::cout << std::endl ;
        std::cout << "max value :" << x.max() << std::endl ;
        std::cout << "min value :" << x.min() << std::endl ;
        std::cout << "max stress11 :" << std::setprecision(16) << stress11.max() << std::endl ;
        std::cout << "min stress11 :" <<  std::setprecision(16) << stress11.min() << std::endl ;
        std::cout << "noise stress11 :" <<  std::setprecision(16) << stress11.max()-stress11.min() << std::endl ;
        std::cout << "max stress22 :" << std::setprecision(16) <<  stress22.max() << std::endl ;
        std::cout << "min stress22 :" <<  std::setprecision(16) << stress22.min() << std::endl ;
        std::cout << "noise stress22 :" <<  std::setprecision(16) << stress22.max()-stress22.min() << std::endl ;
        std::cout << "max stress33 :" <<  std::setprecision(16) << stress33.max() << std::endl ;
        std::cout << "min stress33 :" <<  std::setprecision(16) << stress33.min() << std::endl ;
        std::cout << "noise stress33 :" <<  std::setprecision(16) << stress33.max()-stress33.min() << std::endl ;

        std::cout << "max strain11 :" <<  std::setprecision(16) << strain11.max() << std::endl ;
        std::cout << "min strain11 :" <<  std::setprecision(16) << strain11.min() << std::endl ;
        std::cout << "noise strain11 :" <<  std::setprecision(16) << strain11.max()-strain11.min() << std::endl ;
        std::cout << "max strain22 :" << std::setprecision(16) <<  strain22.max() << std::endl ;
        std::cout << "min strain22 :" <<  std::setprecision(16) << strain22.min() << std::endl ;
        std::cout << "noise strain22 :" <<  std::setprecision(16) << strain22.max()-strain22.min() << std::endl ;
        std::cout << "max strain33 :" <<  std::setprecision(16) << strain33.max() << std::endl ;
        std::cout << "min strain33 :" <<  std::setprecision(16) << strain33.min() << std::endl ;
        std::cout << "noise strain33 :" <<  std::setprecision(16) << strain33.max()-strain33.min() << std::endl ;

        std::cout << "average stress11 : " <<  std::setprecision(16) << average_stress[0]/total_volume << std::endl ;
        std::cout << "average stress22 : " <<  std::setprecision(16) << average_stress[1]/total_volume << std::endl ;
        std::cout << "average stress33 : " << std::setprecision(16) <<  average_stress[2]/total_volume << std::endl ;
        std::cout << "average strain11 : " <<  std::setprecision(16) << average_strain[0]/total_volume << std::endl ;
        std::cout << "average strain22 : " <<  std::setprecision(16) << average_strain[1]/total_volume << std::endl ;
        std::cout << "average strain33 : " <<  std::setprecision(16) << average_strain[2]/total_volume << std::endl ;

        /*        average_strain[1] = (average_strain[1]+average_strain[2])*.5 ;
                average_stress[1] = (average_stress[1]+average_stress[2])*.5 ;*/

        Matrix K(2,2) ;
        K[0][0] = average_strain[0] ;
        K[0][1] = average_strain[1]+average_strain[2] ;
        K[1][0] = average_strain[1] ;
        K[1][1] = average_strain[0]+average_strain[2] ;
        invert2x2Matrix(K) ;
        K.print() ;

        Vector s(2) ;
        s[0] = average_stress[0] ;
        s[1] = average_stress[1] ;

        Vector c = K*s ;

        gettimeofday(&t1, NULL) ;

//         double delta = t1.tv_sec*1000000 - t0.tv_sec*1000000 + t1.tv_usec - t0.tv_usec ;

        std::string filebench("benchmark.txt") ;
        /*        std::fstream out ;
                out.open(filebench.c_str(), std::ios::out|std::ios::app) ;
                out << "ELASTICITY\t" << str_micro ;
                    if(order==2)
                        out << "QUAD" ;
                    out << "\t" << "E_inc = " << prop << "\t"
                    << "dof = " << x.size() << "\t"
                    << "C1111 = " << c[0] << "\t"
                    << "C1122 = " << c[1] << std::endl ;
                out << "Time to solve (s) =" << delta/1e6 << std::endl ;
                out.close() ;*/

        /*        VoxelWriter vw("simulation_out_stiffness", 200) ;
                vw.getField(featureTree, VWFT_STIFFNESS) ;
                vw.write();*/

        Vector tmp(6) ;
       for(double x = 0 ;  x <= size ; x += size*.005)
       {
           for(double y = 0 ;  y <= size ; y += size*.005)
           {
               featureTree->get3DMesh()->getField(STRAIN_FIELD, Point(x,y,halfSize),tmp) ;
                std::cout <<  tmp[0]<< "  "<< std::flush ;
           }
           std::cout << std::endl ;
       }
       for(double x = 0 ;  x <= size ; x += size*.005)
       {
           for(double y = 0 ;  y <= size ; y += size*.005)
           { 
               featureTree->get3DMesh()->getField(REAL_STRESS_FIELD, Point(x,y,halfSize),tmp) ;
                std::cout << tmp[0]<< "  "<< std::flush ;
           }
           std::cout << std::endl ;
       }
        
        
        VoxelWriter vw0("simulation_out_stress", 50) ;
        vw0.getField(featureTree, VWFT_STRESS) ;
        vw0.write();

        VoxelWriter vw1("simulation_out_strain", 50) ;
        vw1.getField(featureTree, VWFT_STRAIN) ;
        vw1.write();

        break ;
    }
    }

    return 0 ;
}
