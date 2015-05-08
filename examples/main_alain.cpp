// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../features/microstructuregenerator.h"
#include "../features/polygonSample.h"
#include "../physics/physics_base.h"
#include "../physics/stiffness.h"
#include "../physics/dual_behaviour.h"
#include "../physics/finite_difference_viscoelasticity.h"
#include "../physics/logarithmic_creep.h"
#include "../physics/homogenization/phase.h"
#include "../physics/homogenization/composite.h"
#include "../physics/logarithmic_creep_with_external_parameters.h"
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
#include "../features/polygonSample.h"
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
#include "../physics/fracturecriteria/spacetimemultilinearsofteningfracturecriterion.h"

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
 	omp_set_num_threads(1) ;
        std::valarray<Point *> pts(4) ;
        pts[0] =  new Point(0,0) ;
        pts[1] =  new Point(0.08,0) ;
        pts[2] =  new Point(0.08,0.08) ;
        pts[3] =  new Point(0,0.04) ;
        
	PolygonalSample s(nullptr, pts) ;
        Sample rect(nullptr, 0.07,0.07,0,0) ;
	Inclusion * inc = new Inclusion( 0.02,0.,0. ) ;
        Inclusion * son = new Inclusion( 0.02, -0.01, 0 ) ;
	rect.setBehaviour(new ViscoElasticOnlyPasteBehaviour() ) ;
	s.setBehaviour(new ElasticOnlyPasteBehaviour() ) ;
	inc->setBehaviour(new ElasticOnlyAggregateBehaviour() ) ;
	son->setBehaviour(new ElasticOnlyAggregateBehaviour(44e9) ) ;

	FeatureTree f(&s) ;
	f.setSamplingNumber(1) ;
//        f.addFeature( &s, inc ) ;
//        f.addFeature( inc, son ) ;

//	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_FORCE_XI, BOTTOM_RIGHT,1e6 ) ) ;
//	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( SET_PROPORTIONAL_DISPLACEMENT_XI_ETA, TOP, -1 ) ) ; // ux = 0.5 u_y
        Point n(-0.004,0.008) ;
	f.addBoundaryCondition( new GeometryAndFaceDefinedSurfaceBoundaryCondition( SET_TANGENT_DISPLACEMENT, dynamic_cast<Polygon*>(&s), n, 0.0001 ) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI, BOTTOM_LEFT) ) ;
	f.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA, BOTTOM ) ) ;

	f.step() ;

//        f.getAssembly()->print() ;

        std::vector<Point *> nodes = f.getNodes() ;
        Vector disp = f.getDisplacements() ;
        for(size_t i = 0 ; i < nodes.size() ; i++)
        {
            if(nodes[i]->getY() > 0.039999+0.5*nodes[i]->getX() )
            {
                double nx = (0.5*disp[i*2] - disp[i*2+1])/2.5 ;
                double vx = disp[i*2] - nx ;
                double vy = 0.5*vx ;
                double ny = disp[i*2+1] - vy ;
                std::cout << nodes[i]->getId() << "\t" << nodes[i]->getX() << "\t" << nodes[i]->getY() << "\t" << disp[i*2] << "\t" << disp[i*2+1] << "\t" << vx << "\t" << vy << "\t" << nx << "\t" << ny << "\t" << std::sqrt( vx*vx+vy*vy ) << std::endl ;
            }
        }

//        f.getAssembly()->print() ;

	TriangleWriter trg( "toto_prop", &f, 1.) ;
	trg.getField( TWFT_STIFFNESS ) ;
	trg.write() ;




	return 0 ;
}

