// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../features/features.h"
#include "../features/sample.h"
#include "../features/crack.h"

#include "../physics/materials/paste_behaviour.h"

#include "../utilities/writer/triangle_writer.h"

using namespace Mu ;

FeatureTree * featureTree ;


int main(int argc, char *argv[])
{
	int nsteps = 2 ;
	double nu = 0.2 ;
	double E_paste = 30e9 ;

	double width = 0.12;
	double height = 0.12;
	Sample sample(width, height , 0., 0.) ;

	featureTree = new FeatureTree(&sample) ;
	
	Point * a = new Point(-0.025,0) ;
	Point * b = new Point(0.025,0) ;
	BranchedCrack * Crack1 = new BranchedCrack( a,  b);
	Crack1->setEnrichementRadius(0.01);

	featureTree->addFeature(&sample, Crack1);//MY

// 	sample.setBehaviour(new OrthotropicStiffness(E_paste, E_paste*.5,  E_paste*.5/(2.*1-nu*0.5),  nu, M_PI*.15)) ;
	sample.setBehaviour(new ElasticOnlyPasteBehaviour(E_paste, nu)) ;

 	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_ETA , TOP, .001)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_ALONG_XI , TOP, 0)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_XI , BOTTOM_RIGHT)) ;
	featureTree->addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(FIX_ALONG_ETA, BOTTOM)) ;

	featureTree->setSamplingNumber(atof(argv[1])) ;
	featureTree->setOrder(QUADRATIC) ;
	featureTree->setMaxIterationsPerStep(1600);

	
	MultiTriangleWriter writerm( "displacements_enrichment", "displacements_enrichment", nullptr ) ;
	for(size_t v = 0 ; v < nsteps ; v++)
	{

		Crack1->print() ;
		bool go_on = featureTree->step() ;

		writerm.reset( featureTree ) ;
		writerm.setGeometry(Crack1->getPrimitive());
		writerm.getField( REAL_STRESS_FIELD ) ;
		writerm.getField( TWFT_INTERSECTION ) ;
		writerm.getField( TWFT_ENRICHMENT ) ;
		writerm.append() ;
		writerm.writeSvg(5, true) ;
		
		if(!go_on)
			break ;
	}
	return 0 ;
}
