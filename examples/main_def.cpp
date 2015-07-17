// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../utilities/samplingcriterion.h"
#include "../features/features.h"
#include "../physics/physics_base.h"
#include "../physics/fracturecriteria/mohrcoulomb.h"
#include "../physics/fracturecriteria/ruptureenergy.h"
#include "../physics/stiffness_with_variable_imposed_deformation.h"
#include "../physics/stiffness_with_variable_imposed_deformation_and_fracture.h"
#include "../physics/weibull_distributed_stiffness.h"
#include "../physics/weibull_distributed_stiffness_with_variable_imposed_deformation_and_fracture.h"
#include "../features/pore.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../features/crack.h"
#include "../features/enrichmentInclusion.h"
#include "../mesher/delaunay_3d.h"
#include "../solvers/assembly.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"

#include <fstream>

#include <cmath>
#include <typeinfo>
#include <limits>

#include <time.h>


using namespace Amie ;

void step(FeatureTree * featureTree )
{
    int nsteps = 1;
    featureTree->setDeltaTime(0.0004);
    featureTree->setMaxIterationsPerStep(128) ;
    for(int i = 0 ; i < nsteps ; i++)
    {
        featureTree->step() ;
        featureTree->printReport(i == 0);
    }
}



int main(int argc, char *argv[])
{

    double nu = 0.3 ;
    double E_agg = 58.9e9 ;
    double E_paste = 12e9 ;

    Matrix m0_agg(3,3) ;
    m0_agg[0][0] = E_agg/(1-nu*nu) ;
    m0_agg[0][1] =E_agg/(1-nu*nu)*nu ;
    m0_agg[0][2] = 0 ;
    m0_agg[1][0] = E_agg/(1-nu*nu)*nu ;
    m0_agg[1][1] = E_agg/(1-nu*nu) ;
    m0_agg[1][2] = 0 ;
    m0_agg[2][0] = 0 ;
    m0_agg[2][1] = 0 ;
    m0_agg[2][2] = E_agg/(1-nu*nu)*(1.-nu)/2. ;

    Matrix m0_paste(3,3) ;
    m0_paste[0][0] = E_paste/(1-nu*nu) ;
    m0_paste[0][1] =E_paste/(1-nu*nu)*nu ;
    m0_paste[0][2] = 0 ;
    m0_paste[1][0] = E_paste/(1-nu*nu)*nu ;
    m0_paste[1][1] = E_paste/(1-nu*nu) ;
    m0_paste[1][2] = 0 ;
    m0_paste[2][0] = 0 ;
    m0_paste[2][1] = 0 ;
    m0_paste[2][2] = E_paste/(1-nu*nu)*(1.-nu)/2. ;

    Sample sample(nullptr, 0.04, 0.04, 0, 0) ;

    FeatureTree F(&sample) ;


    double itzSize = 0.00005;
    int inclusionNumber = 4000 ; // 10 100 500 1000 2000 4000

//     double masseInitiale = .00000743;
//     double densite = 1.;
// 	inclusions = GranuloBolome(masseInitiale, densite, BOLOME_A)(.008, .0001, inclusionNumber, itzSize);
//	inclusions = GranuloBolome(4.79263e-07, 1, BOLOME_D)(.0025, .0001, inclusionNumber, itzSize);
    auto inclusions = PSDGenerator::get2DInclusions(.0025, 4.79263e-07, new PSDBolomeD(), PSDEndCriteria(-1, 0.001, inclusionNumber)) ;

    std::vector<Feature *> feats ;
    for(size_t i = 0; i < inclusions.size() ; i++)
        feats.push_back(inclusions[i]) ;


    inclusions.clear() ;
    for(size_t i = 0; i < feats.size() ; i++)
        inclusions.push_back(static_cast<Inclusion *>(feats[i])) ;


    feats=placement2D(sample.getPrimitive(), feats, 0.000001, 0, 6400);
    double volume = 0 ;
    for(size_t i = 0 ; i < feats.size() ; i++)
        volume += feats[i]->area() ;
    if(!feats.empty())
        std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius()
                  << ", smallest r =" << feats.back()->getRadius()
                  << ", filling = " << volume/sample.area()*100.<< "%"<< std::endl ;

    Vector a(3) ;
    a[0] =  0.0008+.002;
    a[1] =  0.0008+.002 ;
    sample.setBehaviour(new WeibullStiffnessWithVariableImposedDeformationAndFracture(m0_paste,a, new MohrCoulomb(12500000, -12500000*8))) ;
// 	sample.setBehaviour(new Stiffness(m0_paste)) ;
    for(size_t i = 0 ; i < inclusions.size() ; i++)
    {
        inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
        inclusions[i]->setBehaviour(new WeibullDistributedStiffness(E_agg, nu, SPACE_TWO_DIMENSIONAL,-8.*57000000, 57000000)) ;
        F.addFeature(&sample,inclusions[i]) ;
    }
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , TOP_LEFT));
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , BOTTOM_LEFT));
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , BOTTOM_LEFT));
    F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA , BOTTOM_RIGHT));

    std::cout << "largest inclusion with r = " << (*inclusions.begin())->getRadius() << std::endl ;
    std::cout << "smallest inclusion with r = " << (*inclusions.rbegin())->getRadius() << std::endl ;

    F.setSamplingNumber(1400) ;

    F.setOrder(LINEAR) ;

    step(&F) ;

// 	delete dt ;

    return 0 ;
}
