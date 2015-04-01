
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "main.h"
#include "../features/features.h"
#include "../physics/orthotropicstiffness.h"
#include "../features/sample.h"
#include "../features/inclusion.h"
#include "../features/expansiveZone.h"
#include "../utilities/granulo.h"
#include "../utilities/placement.h"
#include "../physics/materials/gel_behaviour.h"
#include "../physics/materials/paste_behaviour.h"
#include "../physics/materials/aggregate_behaviour.h"
#include "../utilities/writer/triangle_writer.h"


#include <fstream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>


using namespace Amie ;

FeatureTree *featureTree ;

double placed_area = 0 ;

double basesize = 0.04 ;

Rectangle baseGeometry( basesize, basesize, 0, 0 ) ;

std::vector<std::pair<ExpansiveZone *, Inclusion *> > zones ;

std::vector<std::pair<double, double> > expansion_reaction ;
std::vector<std::pair<double, double> > expansion_stress_xx ;
std::vector<std::pair<double, double> > expansion_stress_yy ;
std::vector<std::pair<double, double> > apparent_extension ;
std::vector<double> cracked_volume ;
std::vector<double> damaged_volume ;


double nu = 0.3 ;
double E_agg = 58.9e9 ;
double E_paste = 12e9 ;

int totit = 1 ;

double aggregateArea = 0;

GelBehaviour * gel = new GelBehaviour() ;


void step()
{
	int nsteps = 90;
	int nstepstot = 90;
	featureTree->setMaxIterationsPerStep( 400 ) ;

	for( int i = 0 ; i < nsteps ; i++ )
	{
		std::cout << "\r iteration " << i << "/" << nsteps << std::flush ;
		bool go_on = featureTree->step() ;

		if( featureTree->solverConverged() )
		{
			cracked_volume.push_back( featureTree->crackedVolume ) ;
			damaged_volume.push_back( featureTree->damagedVolume ) ;
		}


		Vector x = featureTree->getDisplacements() ;

		std::cout << "unknowns :" << x.size() << std::endl ;


		int npoints = featureTree->get2DMesh()->begin()->getBoundingPoints().size() ;

		std::string filename( "triangles" ) ;

		if( !go_on )
			filename = std::string( "intermediate-triangles" ) ;

		if( !featureTree->solverConverged() )
			filename = std::string( "failed-triangles" ) ;
                
		filename.append( itoa( totit++, 10 ) ) ;
		std::cout << filename << std::endl ;

		TriangleWriter writer(filename.c_str(), featureTree) ;
// 		writer.reset( featureTree ) ;
// 		writer.getField(TWFT_STRAIN_AND_STRESS ) ;
		writer.getField(TWFT_PRINCIPAL_STRESS ) ;
		writer.getField(TWFT_PRINCIPAL_STRAIN ) ;
		writer.getField( TWFT_VON_MISES ) ;
		writer.getField( TWFT_STIFFNESS ) ;
		writer.getField( TWFT_DAMAGE ) ;

		writer.write() ;

		double volume = 0 ;
		double xavg = 0 ;
		double yavg = 0 ;
		for(auto k = featureTree->get2DMesh()->begin() ; k != featureTree->get2DMesh()->end() ; k++)
		{
			if(k->getBehaviour()->type != VOID_BEHAVIOUR )
			{

				if(baseGeometry.in(k->getCenter()))
				{				
					double ar = k->area() ;
					volume += ar ;
					for(int l = 0 ; l < npoints ;l++)
					{
						xavg += x[k->getBoundingPoint(l).getId()*2]*ar/npoints ;
						yavg += x[k->getBoundingPoint(l).getId()*2+1]*ar/npoints ;
					}
				}
			}
		}
			
		xavg /= volume ;
		yavg /= volume ;
		std::pair<Vector, Vector> stempm = featureTree->getFieldMinMax(REAL_STRESS_FIELD) ;
		std::pair<Vector, Vector> etempm = featureTree->getFieldMinMax(STRAIN_FIELD) ;
		std::pair<Vector, Vector> vmm = featureTree->getFieldMinMax(VON_MISES_REAL_STRESS_FIELD) ;
		Vector stemp = featureTree->getAverageField(REAL_STRESS_FIELD) ;
		Vector etemp = featureTree->getAverageField(STRAIN_FIELD) ;
		
		std::cout << std::endl ;
		std::cout << "max value :" << x.max() << std::endl ;
		std::cout << "min value :" << x.min() << std::endl ;
		std::cout << "avg x value :" << xavg << std::endl ;
		std::cout << "avg y value :" << xavg << std::endl ;

		std::cout << "max sigma11 :" << stempm.second[0]  << std::endl ;
		std::cout << "min sigma11 :" << stempm.first[0]   << std::endl ;
		std::cout << "max sigma12 :" << stempm.second[2]  << std::endl ;
		std::cout << "min sigma12 :" << stempm.first[2]   << std::endl ;
		std::cout << "max sigma22 :" << stempm.second[1]  << std::endl ;
		std::cout << "min sigma22 :" << stempm.first[1]   << std::endl ;
		
		std::cout << "max epsilon11 :" << etempm.second[0] << std::endl ;
		std::cout << "min epsilon11 :" << etempm.first[0]  << std::endl ;
		std::cout << "max epsilon12 :" << etempm.second[2] << std::endl ;
		std::cout << "min epsilon12 :" << etempm.first[2]  << std::endl ;
		std::cout << "max epsilon22 :" << etempm.second[1] << std::endl ;
		std::cout << "min epsilon22 :" << etempm.first[1]  << std::endl ;
		
		std::cout << "max von Mises :" << vmm.second[0] << std::endl ;
		std::cout << "min von Mises :" << vmm.first[0] << std::endl ;
		
		std::cout << "average sigma11 : " << stemp[0] << std::endl ;
		std::cout << "average sigma22 : " << stemp[1] << std::endl ;
		std::cout << "average sigma12 : " << stemp[2] << std::endl ;
		std::cout << "average epsilon11 : " << etemp[0] << std::endl ;
		std::cout << "average epsilon22 : " << etemp[1] << std::endl ;
		std::cout << "average epsilon12 : " << etemp[2] << std::endl ;
		
		std::cout << std::endl ;


		if( go_on )
		{
			featureTree->forceEnrichmentChange();
			double delta_r = sqrt( aggregateArea * 0.03 / ( ( double )zones.size() * M_PI ) ) / ( double )nstepstot ;
			double reactedArea = 0 ;

			Inclusion *current = nullptr ;

			if( !zones.empty() )
				current = zones[0].second ;

			double current_area = 0 ;
			int current_number = 0 ;
			int stopped_reaction = 0 ;

			for( size_t z = 0 ; z < zones.size() ; z++ )
			{
				zones[z].first->setRadius( zones[z].first->getRadius() + delta_r ) ;

				if( zones[z].second == current )
				{
					current_area += zones[z].first->area() ;
					current_number++ ;
				}
				else
				{
					if( current_area / zones[z - 1].second->area() > 0.03 )
					{
						stopped_reaction++ ;

						for( int m = 0 ; m < current_number ; m++ )
						{
							reactedArea -= zones[z - 1 - m].first->area() ;
							zones[z - 1 - m].first->setRadius( zones[z].first->getRadius() - delta_r ) ;
							reactedArea += zones[z - 1 - m].first->area() ;
						}
					}

					current_area = zones[z].first->area() ;
					current_number = 1 ;
					current = zones[z].second ;
				}

				reactedArea += zones[z].first->area() ;
			}

			std::cout << "reacted Area : " << reactedArea << ", reaction stopped in " << stopped_reaction << " aggs." << std::endl ;


			std::vector<double> macro_strain = featureTree->getMedianMacroscopicStrain(&baseGeometry) ;
			if( go_on )
			{
				expansion_reaction.push_back( std::make_pair( reactedArea / placed_area, macro_strain[0]) ) ;
				expansion_stress_xx.push_back( std::make_pair( etemp[0], stemp[0] ) ) ;
				expansion_stress_yy.push_back( std::make_pair( etemp[1], stemp[1] ) ) ;
				apparent_extension.push_back( std::make_pair(macro_strain[0], macro_strain[1]) ) ;
			}

		}
			std::cout << "reaction" << "   "
// 			          << expansion_reaction[i].second << "   "
			          << "eps xx" << "\t"
			          << "eps yy" << "\t" 
								<< "sig xx" << "\t"
			          << "sig yy" << "\t"
			          << "   dx "  << "\t"
			          << "   dy "  << "\t"
			          << "cracks"  << "\t"
			          << "damage"  << "\t"
			          << std::endl ;
		for( size_t i = 0 ; i < expansion_reaction.size() ; i++ )
			std::cout << expansion_reaction[i].first << "\t"
// 			          << expansion_reaction[i].second << "\t"
			          << expansion_stress_xx[i].first << "\t"
			          << expansion_stress_yy[i].first << "\t" 
								<< expansion_stress_xx[i].second << "\t"
			          << expansion_stress_yy[i].second << "\t"
			          << apparent_extension[i].first  << "\t"
			          << apparent_extension[i].second  << "\t"
			          << cracked_volume[i]  << "\t"
			          << damaged_volume[i]  << "\t"
			          << std::endl ;

	}

	for( size_t i = 0 ; i < expansion_reaction.size() ; i++ )
		std::cout << expansion_reaction[i].first << "   "
		          << expansion_reaction[i].second << "   "
		          << expansion_stress_xx[i].first << "   "
		          << expansion_stress_xx[i].second << "   "
		          << expansion_stress_yy[i].first << "   "
		          << expansion_stress_yy[i].second << "   "
		          << apparent_extension[i].first  << "   "
		          << apparent_extension[i].second  << "   "
		          << cracked_volume[i]  << "   "
		          << damaged_volume[i]  << "   "
		          << std::endl ;
}

std::vector<std::pair<ExpansiveZone *, Inclusion *> > generateExpansiveZonesHomogeneously( int n, std::vector<Inclusion * > & incs , FeatureTree &F, const Rectangle & sample )
{
	double radiusFraction = 10 ;
	double radius = 0.00001 ;
	std::vector<std::pair<ExpansiveZone *, Inclusion *> > ret ;
	aggregateArea = 0 ;

	std::vector<ExpansiveZone *> zonesToPlace ;

	srand(1);
	for( int i = 0 ; i < n ; i++ )
	{
		Point pos( ( ( double )rand() / RAND_MAX - .5 ) * ( sample.width() - radius * radiusFraction ), ( ( double )rand() / RAND_MAX - .5 ) * ( sample.height() - radius * radiusFraction ) ) ;
		bool alone  = true ;

		for( size_t j = 0 ; j < zonesToPlace.size() ; j++ )
		{
			if( squareDist( pos, zonesToPlace[j]->Circle::getCenter() ) < ( radius * radiusFraction + radius * radiusFraction ) * ( radius * radiusFraction + radius * radiusFraction ) )
			{
				alone = false ;
				break ;
			}
		}

		if( alone )
			zonesToPlace.push_back( new ExpansiveZone( nullptr, radius/15., pos.getX(), pos.getY(), gel ) ) ;
		else
			i-- ;
	}

	std::map<Inclusion *, int> zonesPerIncs ;

	for( size_t i = 0 ; i < zonesToPlace.size() ; i++ )
	{
		bool placed = false ;

		for( size_t j = 0 ; j < incs.size() ; j++ )
		{
			if( dist( zonesToPlace[i]->getCenter(), incs[j]->getCenter() ) < incs[j]->getRadius() - radius * radiusFraction /*&& incs[j]->getRadius() <= 0.008 && incs[j]->getRadius() > 0.004*/ && baseGeometry.in( zonesToPlace[i]->getCenter() ) )
			{
				zonesPerIncs[incs[j]]++ ; ;
				F.addFeature( incs[j], zonesToPlace[i] ) ;
				ret.push_back( std::make_pair( zonesToPlace[i], incs[j] ) ) ;
				placed = true ;
				break ;
			}
		}

		if( !placed )
			delete zonesToPlace[i] ;
	}

	int count = 0 ;

	for( auto i = zonesPerIncs.begin() ; i != zonesPerIncs.end() ; ++i )
	{
		aggregateArea += i->first->area() ;
		count += i->second ;
// 		std::cout << aggregateArea << "  " << count << std::endl ;
	}

	std::cout << "initial Reacted Area = " << M_PI *radius *radius *ret.size() << " in " << ret.size() << " zones" << std::endl ;
	std::cout << "Reactive aggregate Area = " << aggregateArea << std::endl ;
	return ret ;
}

int main( int argc, char *argv[] )
{

	int nzones = atof( argv[1] ) ;
	double dmax = atof( argv[2] ) ;
	double fact = atof( argv[3] ) ;
	double fact0 = atof(argv[4]) ;
	
	double restraintDepth = 0.01 ;
	if(fact0 < 10 && fact < 10)
		restraintDepth = 0 ;
	Sample sample( nullptr, basesize + restraintDepth, basesize + restraintDepth, 0, 0 ) ;

	Matrix m0_agg( 3, 3 ) ;
	m0_agg[0][0] = E_agg / ( 1 - nu * nu ) ;
	m0_agg[0][1] = E_agg / ( 1 - nu * nu ) * nu ;
	m0_agg[0][2] = 0 ;
	m0_agg[1][0] = E_agg / ( 1 - nu * nu ) * nu ;
	m0_agg[1][1] = E_agg / ( 1 - nu * nu ) ;
	m0_agg[1][2] = 0 ;
	m0_agg[2][0] = 0 ;
	m0_agg[2][1] = 0 ;
	m0_agg[2][2] = E_agg / ( 1 - nu * nu ) * ( 1. - nu ) / 2. ;

	Matrix m0_paste( 3, 3 ) ;
	m0_paste[0][0] = E_paste / ( 1 - nu * nu ) ;
	m0_paste[0][1] = E_paste / ( 1 - nu * nu ) * nu ;
	m0_paste[0][2] = 0 ;
	m0_paste[1][0] = E_paste / ( 1 - nu * nu ) * nu ;
	m0_paste[1][1] = E_paste / ( 1 - nu * nu ) ;
	m0_paste[1][2] = 0 ;
	m0_paste[2][0] = 0 ;
	m0_paste[2][1] = 0 ;
	m0_paste[2][2] = E_paste / ( 1 - nu * nu ) * ( 1. - nu ) / 2. ;

	Matrix m0_support( 3, 3 ) ;
	m0_support[0][0] = 1. / ( 1 - 0.2 * 0.2 ) ;
	m0_support[0][1] = 1. / ( 1 - 0.2 * 0.2 ) * 0.2 ;
	m0_support[0][2] = 0 ;
	m0_support[1][0] = 1. / ( 1 - 0.2 * 0.2 ) * 0.2 ;
	m0_support[1][1] = 1. / ( 1 - 0.2 * 0.2 ) ;
	m0_support[1][2] = 0 ;
	m0_support[2][0] = 0 ;
	m0_support[2][1] = 0 ;
	m0_support[2][2] = 1. / ( 1 - 0.2 * 0.2 ) * ( 1. - 0.2 ) / 2. ;

// 	Sample reinforcement0(nullptr, 8,.15,0,.5) ;
// 	reinforcement0.setBehaviour(new Stiffness(m0*5)) ;
//
// 	Sample reinforcement1(nullptr, 8,.15,0,-.5) ;
// 	reinforcement1.setBehaviour(new Stiffness(m0*5)) ;

	FeatureTree F( &sample ) ;
	featureTree = &F ;


	double itzSize = 0.00002;
//	int inclusionNumber = 10 ;
 	int inclusionNumber = 8172 ;

	std::vector<Feature *> feats  = PSDGenerator::get2DConcrete(&F, nullptr,  inclusionNumber, dmax*0.5, itzSize, new PSDBolomeA(), CIRCLE, 1., M_PI, 100000, 0.8, &baseGeometry) ;
	std::vector<Inclusion *> inclusions ;
	
	for( size_t i = 0; i < feats.size() ; i++ )
		inclusions.push_back( static_cast<Inclusion *>( feats[i] ) ) ;

	Rectangle placeGeometry( basesize, basesize, 0, 0 ) ;
	int nAgg = 1 ;
	feats = placement( &placeGeometry, feats, &nAgg, 1, 6400 );
	double volume = 0 ;

	for( size_t i = 0 ; i < feats.size() ; i++ )
		volume += feats[i]->area() ;

	if( !feats.empty() )
		std::cout << "n = " << feats.size() << ", largest r = " << feats.front()->getRadius() - itzSize
		          << ", smallest r =" << feats.back()->getRadius() - itzSize << ", ratio = " << (feats.front()->getRadius() - itzSize) / (feats.back()->getRadius() - itzSize) <<std::endl ;

	sample.setBehaviour( new PasteBehaviour()) ;
	

// 		sample.setBehaviour(new Stiffness(m0_paste)) ;
	Vector setExpansion(0., 3) ; setExpansion[0] = -0.0015 ; setExpansion[1] = -0.0015 ; setExpansion[2] = 0 ;
// 	sample.setBehaviour(new StiffnessWithImposedDeformation(m0_paste, setExpansion)) ;
	if( restraintDepth > 0 )
	{
		Sample *voidtop = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().getX() - ( sample.width() - restraintDepth )*.5 - restraintDepth * .25, sample.getCenter().getY() + ( sample.height() - restraintDepth )*.5 + 0.0025 ) ;
		voidtop->isVirtualFeature = true ;
		voidtop->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidtop );
		F.setSamplingFactor(voidtop, .5);
		Sample *voidbottom = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().getX() - ( sample.width() - restraintDepth )*.5 - restraintDepth * .25, sample.getCenter().getY() - ( sample.height() - restraintDepth )*.5 - 0.0025 ) ;
		voidbottom->isVirtualFeature = true ;
		voidbottom->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidbottom );
		F.setSamplingFactor(voidbottom, .5);
		Sample *voidleft = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().getX() + ( sample.width() - restraintDepth )*.5 + restraintDepth * .25, sample.getCenter().getY() + ( sample.height() - restraintDepth )*.5 + 0.0025 ) ;
		voidleft->isVirtualFeature = true ;
		voidleft->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidleft );
		F.setSamplingFactor(voidleft, .5);
		Sample *voidright = new Sample( nullptr, restraintDepth * .5, restraintDepth * .5, sample.getCenter().getX() + ( sample.width() - restraintDepth )*.5 + restraintDepth * .25, sample.getCenter().getY() - ( sample.height() - restraintDepth )*.5 - 0.0025 ) ;
		voidright->isVirtualFeature = true ;
		voidright->setBehaviour( new VoidForm() );
		F.addFeature( &sample, voidright );
		F.setSamplingFactor(voidright, .5);
                
		//width are  6544984695	10226538586	14726215564      done: 11 13 10 20 12
		//length are 5113269293	26179938780	40906154344      next: 12

		Sample *blocktop = new Sample( nullptr, sample.width() - restraintDepth, restraintDepth * .5, sample.getCenter().getX(), sample.getCenter().getY() + ( sample.height() - restraintDepth )*.5 + restraintDepth * .25 ) ;
		if(fact0 > 10)
			blocktop->setBehaviour(new OrthotropicStiffness(fact0*1e-4, fact0, fact0*1e-4*fact0/(fact0+fact0),  0.1, 0.) ) ;
		else
			blocktop->setBehaviour(new VoidForm()) ;
		
		F.addFeature( nullptr, blocktop );
		F.setSamplingFactor(blocktop, 0.5);

		Sample *blockbottom = new Sample( nullptr, sample.width() - restraintDepth, restraintDepth * .5, sample.getCenter().getX(), sample.getCenter().getY() - ( sample.height() - restraintDepth )*.5 - restraintDepth * .25 ) ;
		if(fact0 > 10)
			blockbottom->setBehaviour(new OrthotropicStiffness(fact0*1e-4, fact0, fact0*1e-4*fact0/(fact0+fact0),  0.1, 0.) ) ;
		else
			blockbottom->setBehaviour(new VoidForm()) ;
		
		F.addFeature( nullptr, blockbottom );
		F.setSamplingFactor(blockbottom, 0.5);

		Sample *blockleft = new Sample( nullptr, restraintDepth * .5, sample.height() - restraintDepth, sample.getCenter().getX() - ( sample.width() - restraintDepth )*.5 - restraintDepth * .25, sample.getCenter().getY() ) ;
		if(fact > 10)
			blockleft->setBehaviour(new OrthotropicStiffness(fact*1e-4, fact, fact*1e-4*fact/(fact+fact),  0.1, 0.) ) ;
		else
			blockleft->setBehaviour(new VoidForm()) ;
		
		F.addFeature( nullptr, blockleft );
		F.setSamplingFactor(blockleft, 0.5);

		Sample *blockright = new Sample( nullptr, restraintDepth * .5, sample.height() - restraintDepth, sample.getCenter().getX() + ( sample.width() - restraintDepth )*.5 + restraintDepth * .25, sample.getCenter().getY() ) ;
		if(fact > 10)
			blockright->setBehaviour(new OrthotropicStiffness(fact*1e-4, fact, fact*1e-4*fact/(fact+fact),  0.1, 0.) ) ;
		else
			blockright->setBehaviour(new VoidForm()) ;
		F.addFeature( nullptr, blockright );
		F.setSamplingFactor(blockright, 0.5);
	}

	std::vector<Inclusion *> placedinclusions ;

	for( size_t i = 0 ; i < feats.size() ; i++ )
	{
		Point a( inclusions[i]->getCenter() + Point( 0, inclusions[i]->getRadius() ) ) ;
		Point b( inclusions[i]->getCenter() + Point( 0, -inclusions[i]->getRadius() ) ) ;
		Point c( inclusions[i]->getCenter() + Point( inclusions[i]->getRadius(), 0 ) ) ;
		Point d( inclusions[i]->getCenter() + Point( -inclusions[i]->getRadius(), 0 ) ) ;

		if( !( !baseGeometry.in( a ) && !baseGeometry.in( b ) && !baseGeometry.in( c ) && !baseGeometry.in( d ) ) )
		{
			inclusions[i]->setRadius(inclusions[i]->getRadius()-itzSize) ;
			AggregateBehaviour *stiff = new AggregateBehaviour() ;
// 			StiffnessWithImposedDeformation * stiff = new StiffnessWithImposedDeformation(m0_agg, setExpansion) ;
// 			Stiffness * stiff = new Stiffness(m0_agg) ;
			// 		stiff->variability = .5 ;
			inclusions[i]->setBehaviour( stiff ) ;
			F.addFeature( &sample, inclusions[i] ) ;
			F.setSamplingFactor(inclusions[i], 4.);
			placed_area += inclusions[i]->area() ;
			placedinclusions.push_back( inclusions[i] );
		}
	}

	inclusions.erase( inclusions.begin() + feats.size(), inclusions.end() ) ;
	std::cout << ", filling = " << placed_area / baseGeometry.area() * 100. << "%" << std::endl ;

	if( !inclusions.empty() )
	{
		std::cout << "largest inclusion with r = " << ( *inclusions.begin() )->getRadius() << std::endl ;
		std::cout << "smallest inclusion with r = " << ( *inclusions.rbegin() )->getRadius() << std::endl ;
		std::cout << "placed area = " <<  placed_area << std::endl ;
	}

	zones = generateExpansiveZonesHomogeneously(nzones, placedinclusions, F , sample) ;
	F.setSamplingNumber( 80 ) ;

	if( restraintDepth > 0 )
	{
// 		F.addBoundaryCondition(new GeometryDefinedBoundaryCondition(FIX_ALONG_XI, new Rectangle(0.035+restraintDepth*.5, basesize+restraintDepth*1.1, -(0.035+restraintDepth*.5)*.5, 0))) ;
// 		F.addBoundaryCondition(new GeometryDefinedBoundaryCondition(FIX_ALONG_ETA, new Rectangle(basesize+restraintDepth*1.1,0.035+restraintDepth*.5, 0,-(0.035+restraintDepth*.5)*.5))) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI , LEFT ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI , RIGHT ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA , BOTTOM ) ) ;
 		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA , TOP ) ) ;
	}
	else
	{
// 		F.addBoundaryCondition(new BoundingBoxDefinedBoundaryCondition(SET_STRESS_ETA , TOP, -15e6)) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_XI , BOTTOM_LEFT ) ) ;
		F.addBoundaryCondition( new BoundingBoxDefinedBoundaryCondition( FIX_ALONG_ETA , BOTTOM ) ) ;
	}
	F.setSamplingFactor(&sample, 2.);
	F.setOrder( LINEAR ) ;
//
	step() ;



	return 0 ;
}
