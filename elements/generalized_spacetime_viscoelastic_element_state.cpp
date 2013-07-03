#include "generalized_spacetime_viscoelastic_element_state.h"
#include "../physics/viscoelasticity.h"
#include "../physics/damagemodels/damagemodel.h"
#include "../mesher/delaunay.h"
#include "../utilities/random.h"
#include "../features/inclusion.h"

using namespace Mu ;

GeneralizedSpaceTimeViscoElasticElementState::GeneralizedSpaceTimeViscoElasticElementState(IntegrableEntity * e) : ElementState(e)
{

  
}

GeneralizedSpaceTimeViscoElasticElementState::GeneralizedSpaceTimeViscoElasticElementState(const GeneralizedSpaceTimeViscoElasticElementState &s) : ElementState(s)
{

}


GeneralizedSpaceTimeViscoElasticElementState & GeneralizedSpaceTimeViscoElasticElementState::operator =(const GeneralizedSpaceTimeViscoElasticElementState & s) 
{
	ElementState::operator =(s) ;
	
	return *this ;
}

GaussPointArray genEquivalentGaussPointArray( TriElement * trg, double time)
{
	GaussPointArray gp ;
	switch(trg->getOrder())
	{
		case LINEAR_TIME_LINEAR:
		case LINEAR_TIME_QUADRATIC:
			gp = gaussPointSet(LINEAR, trg) ;
			break ;
		case QUADRATIC_TIME_LINEAR:
		case QUADRATIC_TIME_QUADRATIC:
			gp = gaussPointSet(QUADRATIC, trg) ;
			break ;
	}
	if(trg->getEnrichmentFunctions().size() > 0)
	{
		std::vector<std::pair<Point, double> > gp_alternative ;
		GaussPointArray original = trg->getGaussPoints() ;
		for(size_t i = 0 ; i < original.gaussPoints.size()/3 ; i++)
		{
			gp_alternative.push_back(original.gaussPoints[i*3]);
			gp_alternative.back().first.t = time ;
			gp_alternative.back().second /= 5./9 ;
		}
		
// 			Point A(0,1) ;
// 			Point B(0,0) ;
// 			Point C(1,0) ;
// 			double a = std::sqrt(0.6) ;
// 			
// 			TriangularInclusion triangular(A,B,C) ;
// 			triangular.sample(16) ;
// 			DelaunayTree * dt = new DelaunayTree(&A,&B,&C) ;
// // 			std::cout << trg.getInPoints().size() << std::endl ;
// // 			std::cout << trg.getBoundingPoints().size() << std::endl ;
// 			for(size_t i = 0 ; i < triangular.getBoundingPoints().size() ; i++)
// 			{
// 				dt->insert(&triangular.getBoundingPoint(i));
// 			}
// 			for(size_t i = 0 ; i < triangular.getInPoints().size() ; i++)
// 			{
// 				dt->insert(&triangular.getInPoint(i));
// 			}
// 			std::vector<DelaunayTriangle *> tris = dt->getElements() ;
// 			for(size_t i = 0 ; i < tris.size() ; i++)
// 			{
// 				Point c = tris[i]->getCenter() ;
// 				double ar = 2.*tris[i]->area() ;
// 				gp_alternative.push_back(std::make_pair(Point(c.x,c.y,0,-a), ar*5./18));
// 				gp_alternative.push_back(std::make_pair(Point(c.x,c.y,0,0), ar*8./18));
// 				gp_alternative.push_back(std::make_pair(Point(c.x,c.y,0,a), ar*5./18));
// 				
// 			}
// /*			
// 			
// 			while(npoints > 0)
// 			{
// 				double x = (double)random()/RAND_MAX ;
// 				double y = (double)random()/RAND_MAX ;
// 				if(x+y < 1)
// 				{
// 					gp_alternative.push_back(std::make_pair(Point(x,y,0,-a), 5./18));
// 					gp_alternative.push_back(std::make_pair(Point(x,y,0,0), 8./18));
// 					gp_alternative.push_back(std::make_pair(Point(x,y,0,a), 5./18));
// 					npoints-- ;
// 				}
// 			}*/
// 			double jac = 2.*trg->area() ;
// 			
// 			for(size_t i = 0 ; i < gp_alternative.size() ; i++)
// 			{
// 				gp_alternative[i].second *= jac ;
// 			}			
		gp.gaussPoints.resize(gp_alternative.size()) ;
		std::copy(gp_alternative.begin(), gp_alternative.end(), &gp.gaussPoints[0]);
// 		gp.id = -1 ;
// 		delete dt ;

// 			double w1 = 0 ;
// 			for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
// 				w1 += gp.gaussPoints[i].second ;
// 			double w2 = 0 ;
// 			for(size_t i = 0 ; i < gp_alternative.size() ; i++)
// 				w2 += gp_alternative[i].second ;
// 			
// 			std::cout << w1 << "\t" << w2 << std::endl ;

	}
	
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		gp.gaussPoints[i].first.t = time ;
	return gp ;
}

void GeneralizedSpaceTimeViscoElasticElementState::getAverageField( FieldType f, Vector & ret, int dummy , double t) 
{
	GaussPointArray gp = parent->getGaussPoints() ;
	ret = 0 ;
	double total = 0 ;
	if(dummy<0)
	{
		gp = genEquivalentGaussPointArray( dynamic_cast<TriElement *>(parent), t) ;
	}
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Point p_ = gp.gaussPoints[i].first ;
		double w = gp.gaussPoints[i].second ;
		Vector tmp(0., ret.size()) ;
		getField(f, p_, tmp, true, dummy) ;
		ret += tmp * w ;
		total += w ;
	}	
	
	double shapes = 0. ;
	for(size_t i = 0 ; i < getParent()->getShapeFunctions().size() ; i++)
	{
		shapes += VirtualMachine().ieval( getParent()->getShapeFunction(i), gp) ;
	}
	double enriched = shapes ;
	for(size_t i = 0 ; i < getParent()->getEnrichmentFunctions().size() ; i++)
	{
		enriched += VirtualMachine().ieval( getParent()->getEnrichmentFunction(i), gp) ;
		//std::cout << "(" << gp.gaussPoints.size() << ")" << enriched << "\t" ;
	}
	
// 	if(parent->getEnrichmentFunctions().size())
// 	{
// 	for(size_t i = 0 ; i < parent->getEnrichmentFunctions().size() ; i++)
// 	  std::cout << parent->getEnrichmentFunction(i).getDofID() << ";" ;
// 	std::cout << std::endl ;
// 	}
//	std::cout << gp.gaussPoints.size() << "\t" << ret.size() << "\t" << getParent()->area() << std::endl ;
	ret /=  (total);//(shapes/enriched) ;
	
// 	if(f == STRAIN_FIELD)
// 	{
// 		Vector toto(3) ; 
// 		getAverageField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, toto, dummy, t) ;
// 		std::cout << parent->getEnrichmentFunctions().size() << ";" << ret[0]-toto[0] << "\t" ;
// 	}

}

void GeneralizedSpaceTimeViscoElasticElementState::getAverageField( FieldType f1, FieldType f2, Vector & r1, Vector & r2, int dummy , double t) 
{	
	getAverageField(f1, r1, dummy, t) ;
	getAverageField(f2, r2, dummy, t) ;
}

void GeneralizedSpaceTimeViscoElasticElementState::getField( FieldType f, const Point & p, Vector & ret, bool local, int )  const 
{
	ret = 0. ; 
	
	VirtualMachine vm ;
	int totaldof = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;	
	int realdof = parent->spaceDimensions() ;
	int blocks = totaldof / realdof ;
	Point p_ = p ;
	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;
	
	Form * visco = (parent->getBehaviour()) ;
	
	switch(f)
	{
		case DISPLACEMENT_FIELD:
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f =  vm.eval( parent->getShapeFunction( j ) , p_) ;
				for(size_t k = 0 ; k < realdof ; k++)
					ret[k] += f * displacements[j*totaldof+k] ;
			}
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm.eval( parent->getEnrichmentFunction( j ) , p_) ;
				for(size_t k = 0 ; k < realdof ; k++)
					ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
			}
			return ;
		case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD:
//			std::cout << ret.size() << "\t" << totaldof << std::endl ;
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f =  vm.eval( parent->getShapeFunction( j ) , p_) ;
				for(size_t k = 0 ; k < totaldof ; k++)
					ret[k] += f * displacements[j*totaldof+k] ;
			}
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm.eval( parent->getEnrichmentFunction( j ) , p_) ;
				for(size_t k = 0 ; k < totaldof ; k++)
					ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
			}
			return ;
		case ENRICHED_DISPLACEMENT_FIELD:
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm.eval( parent->getEnrichmentFunction( j ) , p_) ;
				for(size_t k = 0 ; k < realdof ; k++)
					ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
			}
			return ;
		case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD:
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm.eval( parent->getEnrichmentFunction( j ) , p_) ;
				for(size_t k = 0 ; k < totaldof ; k++)
					ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
			}
			return ;
		case SPEED_FIELD:
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f =  vm.deval( parent->getShapeFunction( j ) , TIME_VARIABLE, p_) ;
				for(size_t k = 0 ; k < realdof ; k++)
					ret[k] += f * displacements[j*totaldof+k] ;
			}
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm.deval( parent->getEnrichmentFunction( j ) , TIME_VARIABLE,  p_) ;
				for(size_t k = 0 ; k < realdof ; k++)
					ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
			}
			return ;
		case GENERALIZED_VISCOELASTIC_SPEED_FIELD:
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f =  vm.deval( parent->getShapeFunction( j ) , TIME_VARIABLE, p_) ;
				for(size_t k = 0 ; k < totaldof ; k++)
					ret[k] += f * displacements[j*totaldof+k] ;
			}
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm.deval( parent->getEnrichmentFunction( j ) , TIME_VARIABLE, p_) ;
				for(size_t k = 0 ; k < totaldof ; k++)
					ret[k] += f * enrichedDisplacements[j*totaldof+k] ;
			}
			return ;
		case STRAIN_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			  
				double x_xi = 0;
				double x_eta = 0;
				double x_tau = 0 ;
				double y_xi = 0;
				double y_eta = 0;
				double y_tau = 0 ;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_tau = vm.deval( parent->getShapeFunction( j ), TIME_VARIABLE, p_ ) ;
					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					x_tau += f_tau * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;
					y_tau += f_tau * displacements[j * totaldof + 1] ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;

					x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
					x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
					y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
					y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2] ;
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] ) ;// + x_tau * Jinv[1][2]  + y_tau * Jinv[0][2] );
								
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getShapeFunction( j ), ZETA, p_ ) ;
					double x = displacements[j * totaldof] ;
					double y = displacements[j * totaldof + 1] ;
					double z = displacements[j * totaldof + 2] ;

					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getEnrichmentFunction( j ), ZETA, p_ ) ;
					double x = enrichedDisplacements[j * totaldof] ;
					double y = enrichedDisplacements[j * totaldof + 1] ;
					double z = enrichedDisplacements[j * totaldof + 2] ;

					x_xi += f_xi * x;
					x_eta += f_eta * x ;
					x_zeta += f_zeta * x ;
					y_xi += f_xi * y ;
					y_eta += f_eta * y ;
					y_zeta += f_zeta * y ;
					z_xi += f_xi * z ;
					z_eta += f_eta * z ;
					z_zeta += f_zeta * z ;
				}

				Matrix Jinv( 3, 3 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
				ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

				ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
			}
			return ;
		case GENERALIZED_VISCOELASTIC_STRAIN_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			 
				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
				Vector dt(0., totaldof) ;
				
				double x_xi = 0;
				double x_eta = 0;
				double x_tau = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_tau = 0;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_tau = vm.deval( parent->getShapeFunction( j ), TIME_VARIABLE, p_ ) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * displacements[j * totaldof + i] ;
						dy[i] += f_eta * displacements[j * totaldof + i] ;
						dt[i] += f_tau * displacements[j * totaldof + i] ;
					}
					
/*					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					double f_tau = vm.deval( parent->getEnrichmentFunction( j ), TIME_VARIABLE, p_ ) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * enrichedDisplacements[j * totaldof + i] ;
						dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
						dt[i] += f_tau * enrichedDisplacements[j * totaldof + i] ;
					}

/*					x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
					x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
					y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
					y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;*/
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi  = dx[ i * realdof + 0 ] ;
					x_eta = dy[ i * realdof + 0 ] ;
					x_tau = dt[ i * realdof + 0 ] ;
					y_xi  = dx[ i * realdof + 1 ] ;
					y_eta = dy[ i * realdof + 1 ] ;
					y_tau = dt[ i * realdof + 1 ] ;
					ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2];
					ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
					ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
//					std::cout << ret.size() << std::endl ;
				}
				
// 				Vector strainns(3) ;
// 				this->getField(STRAIN_FIELD, p_, strainns, true);
// 				
// 				std::cout << parent->getEnrichmentFunctions().size()<<"-" << ret[0] - strainns[0] << " ; " ;
				
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
				Vector dz(0., totaldof) ;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getShapeFunction( j ), ZETA, p_ ) ;

					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi   * displacements[j * totaldof + i] ;
						dy[i] += f_eta  * displacements[j * totaldof + i] ;
						dz[i] += f_zeta * displacements[j * totaldof + i] ;
					}

/*					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;*/
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getEnrichmentFunction( j ), ZETA, p_ ) ;
					
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi   * enrichedDisplacements[j * totaldof + i] ;
						dy[i] += f_eta  * enrichedDisplacements[j * totaldof + i] ;
						dz[i] += f_zeta * enrichedDisplacements[j * totaldof + i] ;
					}

/*					x_xi += f_xi * x;
					x_eta += f_eta * x ;
					x_zeta += f_zeta * x ;
					y_xi += f_xi * y ;
					y_eta += f_eta * y ;
					y_zeta += f_zeta * y ;
					z_xi += f_xi * z ;
					z_eta += f_eta * z ;
					z_zeta += f_zeta * z ;*/
				}

				Matrix Jinv( 4, 4 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi   = dx[ i * realdof + 0 ] ;
					x_eta  = dy[ i * realdof + 0 ] ;
					x_zeta = dz[ i * realdof + 0 ] ;
					y_xi   = dx[ i * realdof + 1 ] ;
					y_eta  = dy[ i * realdof + 1 ] ;
					y_zeta = dz[ i * realdof + 1 ] ;
					z_xi   = dx[ i * realdof + 2 ] ;
					z_eta  = dy[ i * realdof + 2 ] ;
					z_zeta = dz[ i * realdof + 2 ] ;
					ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
					ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
					ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

					ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
							    ( y_eta ) * Jinv[2][1] +
							    ( y_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[1][0] +
							    ( z_eta ) * Jinv[1][1] +
							    ( z_zeta ) * Jinv[1][2] );

					ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
							    ( x_eta ) * Jinv[2][1] +
							    ( x_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[0][0] +
							    ( z_eta ) * Jinv[0][1] +
							    ( z_zeta ) * Jinv[0][2] );

					ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
							    ( y_eta )  * Jinv[0][1] +
							    ( y_zeta ) * Jinv[0][2] +
							    ( x_xi )   * Jinv[1][0] +
							    ( x_eta )  * Jinv[1][1] +
							    ( x_zeta ) * Jinv[1][2] );
				}
				
			}
			return ;
		case PRINCIPAL_STRAIN_FIELD:
		{
			Vector strains(0.,3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			this->getField(STRAIN_FIELD, p_, strains, true) ;
			ret = toPrincipal(strains) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_PRINCIPAL_STRAIN_FIELD:
		{
			Vector strains(0.,blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true) ;
			Vector tmp(0., 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			Vector ptmp(0., 2+(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			for(size_t i = 0 ; i < blocks ; i++)
			{
				for(size_t j = 0 ; j < tmp.size() ; j++)
					tmp[j] = strains[ i*tmp.size() + j ] ;
				ptmp = toPrincipal(tmp) ;
				for(size_t j = 0 ; j < tmp.size() ; j++)
					ret[ i*ptmp.size() + j] = ptmp[j] ;				
			}
			return ;
		}
		case NON_ENRICHED_STRAIN_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			  
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getShapeFunction( j ), ZETA, p_ ) ;
					double x = displacements[j * totaldof] ;
					double y = displacements[j * totaldof + 1] ;
					double z = displacements[j * totaldof + 2] ;

					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;
				}

				Matrix Jinv( 3, 3 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
				ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

				ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
			}
			return ;
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			 
				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
			  
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * displacements[j * totaldof + i] ;
						dy[i] += f_eta * displacements[j * totaldof + i] ;
					}
					
/*					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi  = dx[ i * realdof + 0 ] ;
					x_eta = dy[ i * realdof + 0 ] ;
					y_xi  = dx[ i * realdof + 1 ] ;
					y_eta = dy[ i * realdof + 1 ] ;
					ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
					ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
					ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
				}
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
				Vector dz(0., totaldof) ;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getShapeFunction( j ), ZETA, p_ ) ;

					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi   * displacements[j * totaldof + i] ;
						dy[i] += f_eta  * displacements[j * totaldof + i] ;
						dz[i] += f_zeta * displacements[j * totaldof + i] ;
					}

/*					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;*/
				}


				Matrix Jinv( 4, 4 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi   = dx[ i * realdof + 0 ] ;
					x_eta  = dy[ i * realdof + 0 ] ;
					x_zeta = dz[ i * realdof + 0 ] ;
					y_xi   = dx[ i * realdof + 1 ] ;
					y_eta  = dy[ i * realdof + 1 ] ;
					y_zeta = dz[ i * realdof + 1 ] ;
					z_xi   = dx[ i * realdof + 2 ] ;
					z_eta  = dy[ i * realdof + 2 ] ;
					z_zeta = dz[ i * realdof + 2 ] ;
					ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
					ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
					ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

					ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
							    ( y_eta ) * Jinv[2][1] +
							    ( y_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[1][0] +
							    ( z_eta ) * Jinv[1][1] +
							    ( z_zeta ) * Jinv[1][2] );

					ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
							    ( x_eta ) * Jinv[2][1] +
							    ( x_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[0][0] +
							    ( z_eta ) * Jinv[0][1] +
							    ( z_zeta ) * Jinv[0][2] );

					ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
							    ( y_eta )  * Jinv[0][1] +
							    ( y_zeta ) * Jinv[0][2] +
							    ( x_xi )   * Jinv[1][0] +
							    ( x_eta )  * Jinv[1][1] +
							    ( x_zeta ) * Jinv[1][2] );
				}
				
			}
			return ;
		case VON_MISES_STRAIN_FIELD:
		{
			Vector eps(0., (size_t) parent->spaceDimensions()) ;
			this->getField( PRINCIPAL_STRAIN_FIELD, p_, eps, true ) ;
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				ret[0] = ( 2. / 3. * ( eps[0] * eps[0] + eps[1] * eps[1] ) ) ;
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				ret[0] = sqrt( 2. / 3. * ( eps[0] * eps[0] + eps[1] * eps[1] + eps[2] * eps[2] ) ) ;
			}
			return ;
		}
		case STRAIN_RATE_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			  
				double x_xi = 0;
				double x_eta = 0;
				
				double y_xi = 0;
				double y_eta = 0;

				Vector Xdx(3) ; Xdx = 0 ;
				Vector Xdxx(6) ; Xdxx = 0 ;
				Vector Ydx(3) ; Ydx = 0 ;
				Vector Ydxx(6) ; Ydxx = 0 ;
				
				Vector xdx(3) ;
				Vector xdxx(6) ;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ ) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_) ;

					xdx[0] = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					xdx[1] = vm.deval( parent->getShapeFunction( j ), ETA, p_) ;
					xdx[2] = vm.deval( parent->getShapeFunction( j ), TIME_VARIABLE, p_ ) ;
					
					xdxx[0] = vm.ddeval( parent->getShapeFunction( j ), XI, XI, p_ , 1e-5) ;
					xdxx[1] = vm.ddeval( parent->getShapeFunction( j ), ETA, ETA, p_ , 1e-5) ;
					xdxx[2] = vm.ddeval( parent->getShapeFunction( j ), TIME_VARIABLE, TIME_VARIABLE, p_ , 1e-5) ;
					xdxx[3] = vm.ddeval( parent->getShapeFunction( j ), XI, ETA, p_ , 1e-5) ;
					xdxx[4] = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					xdxx[5] = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					
					for(int x = 0 ; x < 3 ; x++)
					{
						Xdx[x] += xdx[x] * displacements[j * totaldof] ;
						Xdxx[x] += xdxx[x] * displacements[j * totaldof] ;
						Xdxx[x+3] += xdxx[x+3] * displacements[j * totaldof] ;

					  
						Ydx[x] += xdx[x] * displacements[j * totaldof + 1] ;
						Ydxx[x] += xdxx[x] * displacements[j * totaldof + 1] ;
						Ydxx[x+3] += xdxx[x+3] * displacements[j * totaldof + 1] ;
					}
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm.ddeval( parent->getEnrichmentFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getEnrichmentFunction( j ), ETA, TIME_VARIABLE,  p_ , 1e-5) ;

					x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
					x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
					y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
					y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;
				}

				Matrix Jinv ;
				Matrix T1, T2 ;
				parent->getSecondJacobianMatrix( p_, T1, T2 ) ;
				
				Vector dXX = (Vector) (T1 * Xdx) + (Vector) (T2 * Xdxx) ;
				Vector dYY = (Vector) (T1 * Ydx) + (Vector) (T2 * Ydxx) ;
				
				ret[0] = dXX[5] ;
				ret[1] = dYY[4] ;
				ret[2] = 0.5 * (dXX[4] + dYY[5]) ;
				
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
// 				
// 				ret *= Jinv[2][2] ;
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					double f_zeta = vm.ddeval( parent->getShapeFunction( j ), ZETA, TIME_VARIABLE, p_ , 1e-5) ;
					double x = displacements[j * totaldof] ;
					double y = displacements[j * totaldof + 1] ;
					double z = displacements[j * totaldof + 2] ;

					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
				{
					double f_xi = vm.ddeval( parent->getEnrichmentFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getEnrichmentFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					double f_zeta = vm.ddeval( parent->getEnrichmentFunction( j ), ZETA, TIME_VARIABLE, p_ , 1e-5) ;
					double x = enrichedDisplacements[j * totaldof] ;
					double y = enrichedDisplacements[j * totaldof + 1] ;
					double z = enrichedDisplacements[j * totaldof + 2] ;

					x_xi += f_xi * x;
					x_eta += f_eta * x ;
					x_zeta += f_zeta * x ;
					y_xi += f_xi * y ;
					y_eta += f_eta * y ;
					y_zeta += f_zeta * y ;
					z_xi += f_xi * z ;
					z_eta += f_eta * z ;
					z_zeta += f_zeta * z ;
				}

				Matrix Jinv( 4, 4 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
				ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

				ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
				
				ret *= Jinv[3][3] ;
			}
			return ;
		case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			 
// 				if(totaldof == 2)
// 				{
// 					this->getField( STRAIN_RATE_FIELD, p_, ret, true) ;
// 					return ;
// 				}
			  
			  
				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
				Vector dt(0., totaldof) ;
			  
				double x_xi = 0;
				double x_eta = 0;
				double x_tau = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_tau = 0;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 10.*default_derivation_delta) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE,p_ , 10.*default_derivation_delta) ;
					double f_tau = vm.ddeval( parent->getShapeFunction( j ), TIME_VARIABLE, TIME_VARIABLE,p_ , 10.*default_derivation_delta) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * displacements[j * totaldof + i] ;
						dy[i] += f_eta * displacements[j * totaldof + i] ;
						dt[i] += f_tau * displacements[j * totaldof + i] ;
					}
					
/*					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm.ddeval( parent->getEnrichmentFunction( j ), XI, TIME_VARIABLE,p_ , 10.*default_derivation_delta) ;
					double f_eta = vm.ddeval( parent->getEnrichmentFunction( j ), ETA, TIME_VARIABLE,p_ , 10.*default_derivation_delta) ;
					double f_tau = vm.ddeval( parent->getEnrichmentFunction( j ), TIME_VARIABLE,TIME_VARIABLE, p_ , 10.*default_derivation_delta) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * enrichedDisplacements[j * totaldof + i] ;
						dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
						dt[i] += f_tau * enrichedDisplacements[j * totaldof + i] ;
					}

/*					x_xi += f_xi * enrichedDisplacements[j * totaldof] ;
					x_eta += f_eta * enrichedDisplacements[j * totaldof] ;
					y_xi += f_xi * enrichedDisplacements[j * totaldof + 1] ;
					y_eta += f_eta * enrichedDisplacements[j * totaldof + 1] ;*/
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi  = dx[ i * realdof + 0 ] ;
					x_eta = dy[ i * realdof + 0 ] ;
					x_tau = dt[ i * realdof + 0 ] ;
					y_xi  = dx[ i * realdof + 1 ] ;
					y_eta = dy[ i * realdof + 1 ] ;
					y_tau = dt[ i * realdof + 1 ] ;
					ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;//+ x_tau * Jinv[0][2];
					ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;//+ y_tau * Jinv[1][2] ;
					ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );//+ x_tau * Jinv[1][2]  + y_tau * Jinv[0][2]);
					//					std::cout << ret.size() << std::endl ;
				}
				ret *= Jinv[2][2] ;
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
				Vector dz(0., totaldof) ;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
					double f_zeta = vm.ddeval( parent->getShapeFunction( j ), ZETA, TIME_VARIABLE,  p_, 1e-5 ) ;

					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi   * displacements[j * totaldof + i] ;
						dy[i] += f_eta  * displacements[j * totaldof + i] ;
						dz[i] += f_zeta * displacements[j * totaldof + i] ;
					}

/*					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;*/
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
				{
					double f_xi = vm.ddeval( parent->getEnrichmentFunction( j ), XI, TIME_VARIABLE, p_, 1e-5 ) ;
					double f_eta = vm.ddeval( parent->getEnrichmentFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					double f_zeta = vm.ddeval( parent->getEnrichmentFunction( j ), ZETA, TIME_VARIABLE, p_, 1e-5 ) ;
					
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi   * enrichedDisplacements[j * totaldof + i] ;
						dy[i] += f_eta  * enrichedDisplacements[j * totaldof + i] ;
						dz[i] += f_zeta * enrichedDisplacements[j * totaldof + i] ;
					}

/*					x_xi += f_xi * x;
					x_eta += f_eta * x ;
					x_zeta += f_zeta * x ;
					y_xi += f_xi * y ;
					y_eta += f_eta * y ;
					y_zeta += f_zeta * y ;
					z_xi += f_xi * z ;
					z_eta += f_eta * z ;
					z_zeta += f_zeta * z ;*/
				}

				Matrix Jinv( 4, 4 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi   = dx[ i * realdof + 0 ] ;
					x_eta  = dy[ i * realdof + 0 ] ;
					x_zeta = dz[ i * realdof + 0 ] ;
					y_xi   = dx[ i * realdof + 1 ] ;
					y_eta  = dy[ i * realdof + 1 ] ;
					y_zeta = dz[ i * realdof + 1 ] ;
					z_xi   = dx[ i * realdof + 2 ] ;
					z_eta  = dy[ i * realdof + 2 ] ;
					z_zeta = dz[ i * realdof + 2 ] ;
					ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
					ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
					ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

					ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
							    ( y_eta ) * Jinv[2][1] +
							    ( y_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[1][0] +
							    ( z_eta ) * Jinv[1][1] +
							    ( z_zeta ) * Jinv[1][2] );

					ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
							    ( x_eta ) * Jinv[2][1] +
							    ( x_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[0][0] +
							    ( z_eta ) * Jinv[0][1] +
							    ( z_zeta ) * Jinv[0][2] );

					ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
							    ( y_eta )  * Jinv[0][1] +
							    ( y_zeta ) * Jinv[0][2] +
							    ( x_xi )   * Jinv[1][0] +
							    ( x_eta )  * Jinv[1][1] +
							    ( x_zeta ) * Jinv[1][2] );
				}
				ret *= Jinv[3][3] ;
				
			}
			return ;
		case NON_ENRICHED_STRAIN_RATE_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			  
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;
				}
				
				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
				
				ret *= Jinv[2][2] ;
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					double f_zeta = vm.ddeval( parent->getShapeFunction( j ), ZETA, TIME_VARIABLE, p_ , 1e-5) ;
					double x = displacements[j * totaldof] ;
					double y = displacements[j * totaldof + 1] ;
					double z = displacements[j * totaldof + 2] ;

					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;
				}

				Matrix Jinv( 4, 4 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
				ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

				ret[3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
				
				ret *= Jinv[3][3] ;
			}
			return ;
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{			 
				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
			  
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * displacements[j * totaldof + i] ;
						dy[i] += f_eta * displacements[j * totaldof + i] ;
					}
					
/*					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;*/
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi  = dx[ i * realdof + 0 ] ;
					x_eta = dy[ i * realdof + 0 ] ;
					y_xi  = dx[ i * realdof + 1 ] ;
					y_eta = dy[ i * realdof + 1 ] ;
					ret[i*3+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
					ret[i*3+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
					ret[i*3+2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
				}
				ret *= Jinv[2][2] ;
// 				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
// 				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
// 				ret[2] = 0.5 * ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;
				double y_xi = 0;
				double y_eta = 0;
				double y_zeta = 0;
				double z_xi = 0;
				double z_eta = 0;
				double z_zeta = 0;

				Vector dx(0., totaldof) ;
				Vector dy(0., totaldof) ;
				Vector dz(0., totaldof) ;
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5 ) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5 ) ;
					double f_zeta = vm.ddeval( parent->getShapeFunction( j ), ZETA, TIME_VARIABLE,  p_, 1e-5 ) ;

					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi   * displacements[j * totaldof + i] ;
						dy[i] += f_eta  * displacements[j * totaldof + i] ;
						dz[i] += f_zeta * displacements[j * totaldof + i] ;
					}

/*					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
					y_xi   += f_xi   * y ;
					y_eta  += f_eta  * y ;
					y_zeta += f_zeta * y ;
					z_xi   += f_xi   * z ;
					z_eta  += f_eta  * z ;
					z_zeta += f_zeta * z ;*/
				}

				Matrix Jinv( 4, 4 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				for(size_t i = 0 ; i < blocks ; i++)
				{
					x_xi   = dx[ i * realdof + 0 ] ;
					x_eta  = dy[ i * realdof + 0 ] ;
					x_zeta = dz[ i * realdof + 0 ] ;
					y_xi   = dx[ i * realdof + 1 ] ;
					y_eta  = dy[ i * realdof + 1 ] ;
					y_zeta = dz[ i * realdof + 1 ] ;
					z_xi   = dx[ i * realdof + 2 ] ;
					z_eta  = dy[ i * realdof + 2 ] ;
					z_zeta = dz[ i * realdof + 2 ] ;
					ret[i*6+0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
					ret[i*6+1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
					ret[i*6+2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

					ret[i*6+3] = 0.5 * ( ( y_xi ) * Jinv[2][0] +
							    ( y_eta ) * Jinv[2][1] +
							    ( y_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[1][0] +
							    ( z_eta ) * Jinv[1][1] +
							    ( z_zeta ) * Jinv[1][2] );

					ret[i*6+4] = 0.5 * ( ( x_xi ) * Jinv[2][0] +
							    ( x_eta ) * Jinv[2][1] +
							    ( x_zeta ) * Jinv[2][2] +
							    ( z_xi ) * Jinv[0][0] +
							    ( z_eta ) * Jinv[0][1] +
							    ( z_zeta ) * Jinv[0][2] );

					ret[i*6+5] = 0.5 * ( ( y_xi )   * Jinv[0][0] +
							    ( y_eta )  * Jinv[0][1] +
							    ( y_zeta ) * Jinv[0][2] +
							    ( x_xi )   * Jinv[1][0] +
							    ( x_eta )  * Jinv[1][1] +
							    ( x_zeta ) * Jinv[1][2] );
				}
				ret *= Jinv[3][3] ;
				
			}
			return ;
		case REAL_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			Vector stresses = (Vector) (visco->getTensor(p_, parent) * strains) 
					+ (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			for(size_t i = 0 ; i < 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL) ; i++)
				ret[i] = stresses[i] ;
/*			if(parent->getEnrichmentFunctions().size() > 0)
			{
				Vector toto = getParent()->getBehaviour()->getImposedStress(p_, parent) ;
				Matrix tata = getParent()->getBehaviour()->getTensor( p_ ) ;
				std::cout << toto[0] << "\t" << tata[0][0] << "\t" << tata[0][1] << std::endl ;
			}*/
			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD:
		{
// 			visco->getTensor(p_,parent).print() ;
// 			visco->getViscousTensor(p_, parent).print() ;
// 			exit(0) ;
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			ret = (Vector) (visco->getTensor(p_, parent) * strains) 
			    + (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			Vector stresses = visco->getImposedStress(p_, parent) ;
			for(size_t i = 0 ; i < stresses.size() ; i++)
				ret[i] -= stresses[i] ;
			return ;
		}
		case PRINCIPAL_REAL_STRESS_FIELD:
		{
			Vector stress(0.,3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			this->getField(REAL_STRESS_FIELD, p_, stress, true) ;
			ret = toPrincipal(stress) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_PRINCIPAL_REAL_STRESS_FIELD:
		{
			Vector stress(0.,blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD, p_, stress, true) ;
			Vector tmp(0., 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			Vector ptmp(0., 2+(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			for(size_t i = 0 ; i < blocks ; i++)
			{
				for(size_t j = 0 ; j < tmp.size() ; j++)
					tmp[j] = stress[ i*tmp.size() + j ] ;
				ptmp = toPrincipal(tmp) ;
				for(size_t j = 0 ; j < tmp.size() ; j++)
					ret[ i*ptmp.size() + j] = ptmp[j] ;				
			}
			return ;
		}
		case NON_ENRICHED_REAL_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			Vector stresses = (Vector) (visco->getTensor(p_, parent) * strains) 
					+ (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			for(size_t i = 0 ; i < 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL) ; i++)
				ret[i] = stresses[i] ;
			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			ret = (Vector) (visco->getTensor(p_, parent) * strains) 
			    + (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			Vector stresses = visco->getImposedStress(p_, parent) ;
			for(size_t i = 0 ; i < stresses.size() ; i++)
				ret[i] -= stresses[i] ;
			return ;
		}
		case VON_MISES_REAL_STRESS_FIELD:
			if( parent->getOrder() == LINEAR )
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector sigma(0., 2) ;
					Point c(1./3., 1./3.) ;
					this->getField(PRINCIPAL_REAL_STRESS_FIELD, c, sigma, true) ;
					ret[0] = sqrt( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
					return ;
				}
				else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
				{
					Mu::PointArray pts( 4 ) ;
					pts[0] = &parent->getBoundingPoint( 0 ) ;
					pts[1] = &parent->getBoundingPoint( 1 ) ;
					pts[2] = &parent->getBoundingPoint( 2 ) ;
					pts[3] = &parent->getBoundingPoint( 3 ) ;
					Vector sigma(0., 24) ;
//					this->getField( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false) ;
					sigma[0] = sqrt( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
					sigma[1] = sqrt( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
					sigma[2] = sqrt( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
					sigma[3] = sqrt( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

					ret[0] = std::max( std::max( std::max( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
					
					return ;
				}
				return ;
			}
			else
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector principalStresses(0., parent->getBoundingPoints().size()*2) ;
//					this->getField(PRINCIPAL_REAL_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false) ;
					double maxS = 0 ;

					for( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
					{
						maxS = std::max( maxS,
								sqrt( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
					}

					ret[0] = maxS ;
					return ;

				}
				else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
				{
					Mu::PointArray pts( 4 ) ;
					pts[0] = &parent->getBoundingPoint( 0 ) ;
					pts[1] = &parent->getBoundingPoint( 2 ) ;
					pts[2] = &parent->getBoundingPoint( 4 ) ;
					pts[3] = &parent->getBoundingPoint( 6 ) ;
					Vector sigma(0., 24) ;
//					this->getField( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false) ;
					sigma[0] = sqrt( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6 ;
					sigma[1] = sqrt( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6 ;
					sigma[2] = sqrt( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6 ;
					sigma[3] = sqrt( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6 ;

					ret[0] = std::max( std::max( std::max( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
					return ;
				}
			}
			return ;
		case EFFECTIVE_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			Vector stresses = (Vector) (visco->param * strains) 
					+ (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			for(size_t i = 0 ; i < 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL) ; i++)
				ret[i] = stresses[i] ;
			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			ret = (Vector) (visco->param * strains) 
			    + (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			Vector stresses = visco->getImposedStress(p_, parent) ;
			for(size_t i = 0 ; i < stresses.size() ; i++)
				ret[i] -= stresses[i] ;
			return ;
		}
		case PRINCIPAL_EFFECTIVE_STRESS_FIELD:
		{
			Vector stress(0.,3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			this->getField(EFFECTIVE_STRESS_FIELD, p_, stress, true) ;
			ret = toPrincipal(stress) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_PRINCIPAL_EFFECTIVE_STRESS_FIELD:
		{
			Vector stress(0.,blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD, p_, stress, true) ;
			Vector tmp(0., 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			Vector ptmp(0., 2+(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			for(size_t i = 0 ; i < blocks ; i++)
			{
				for(size_t j = 0 ; j < tmp.size() ; j++)
					tmp[j] = stress[ i*tmp.size() + j ] ;
				ptmp = toPrincipal(tmp) ;
				for(size_t j = 0 ; j < tmp.size() ; j++)
					ret[ i*ptmp.size() + j] = ptmp[j] ;				
			}
			return ;
		}
		case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			Vector stresses = (Vector) (visco->param * strains) 
					+ (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			for(size_t i = 0 ; i < 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL) ; i++)
				ret[i] = stresses[i] ;
			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			ret = (Vector) (visco->param * strains) 
			    + (Vector) (visco->getViscousTensor(p_, parent) * speeds) ;
			Vector stresses = visco->getImposedStress(p_, parent) ;
			for(size_t i = 0 ; i < stresses.size() ; i++)
				ret[i] -= stresses[i] ;
			return ;
		}
		case VON_MISES_EFFECTIVE_STRESS_FIELD:
			if( parent->getOrder() == LINEAR )
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector sigma(0., 2) ;
					Point c(1./3., 1./3.) ;
					this->getField(PRINCIPAL_EFFECTIVE_STRESS_FIELD, c, sigma, true) ;
					ret[0] = sqrt( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
					return ;
				}
				else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
				{
					Mu::PointArray pts( 4 ) ;
					pts[0] = &parent->getBoundingPoint( 0 ) ;
					pts[1] = &parent->getBoundingPoint( 1 ) ;
					pts[2] = &parent->getBoundingPoint( 2 ) ;
					pts[3] = &parent->getBoundingPoint( 3 ) ;
					Vector sigma(0., 24) ;
//					this->getField( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false) ;
					sigma[0] = sqrt( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
					sigma[1] = sqrt( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
					sigma[2] = sqrt( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
					sigma[3] = sqrt( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

					ret[0] = std::max( std::max( std::max( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
					
					return ;
				}
				return ;
			}
			else
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector principalStresses(0., parent->getBoundingPoints().size()*2) ;
//					this->getField(PRINCIPAL_EFFECTIVE_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false) ;
					double maxS = 0 ;

					for( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
					{
						maxS = std::max( maxS,
								sqrt( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
					}

					ret[0] = maxS ;
					return ;

				}
				else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
				{
					Mu::PointArray pts( 4 ) ;
					pts[0] = &parent->getBoundingPoint( 0 ) ;
					pts[1] = &parent->getBoundingPoint( 2 ) ;
					pts[2] = &parent->getBoundingPoint( 4 ) ;
					pts[3] = &parent->getBoundingPoint( 6 ) ;
					Vector sigma(0., 24) ;
//					this->getField( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false) ;
					sigma[0] = sqrt( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6 ;
					sigma[1] = sqrt( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6 ;
					sigma[2] = sqrt( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6 ;
					sigma[3] = sqrt( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6 ;

					ret[0] = std::max( std::max( std::max( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
					return ;
				}
			}
			return ;
		case PRINCIPAL_ANGLE_FIELD:
		{
			Vector strains(0., 3+3*(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
			this->getField(STRAIN_FIELD,  p_, strains, true ) ;
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				ret[0] =  0.5 * atan2( strains[2], strains[0] - strains[1] ) ;
			else
			{
				ret[0] =  0.5 * atan2(strains[3] , strains[0] - strains[1] ) ;
				ret[1] =  0.5 * atan2(strains[4] , strains[0] - strains[2] ) ;
				ret[2] =  0.5 * atan2(strains[5] , strains[1] - strains[2] ) ;
			}
			return ;
		}
/*		case GRADIENT_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 1 )
			{
				double x_xi = 0;
				double x_eta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					x_xi += f_xi * displacements[j] ;
					x_eta += f_eta * displacements[j] ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size(); j++ )
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					x_xi += f_xi * enrichedDisplacements[j] ;
					x_eta += f_eta * enrichedDisplacements[j] ;

				}

				Matrix Jinv( 2, 2 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
				ret[1] = ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1] ;
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 1 )
			{
				double x_xi = 0;
				double x_eta = 0;
				double x_zeta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getShapeFunction( j ), ZETA, p_ ) ;
					double x = displacements[j] ;

					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					double f_zeta = vm.deval( parent->getEnrichmentFunction( j ), ZETA, p_ ) ;
					double x = enrichedDisplacements[j] ;

					x_xi += f_xi * x;
					x_eta += f_eta * x ;
					x_zeta += f_zeta * x ;
				}

				Matrix Jinv( 3, 3 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
				ret[1] = ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( x_zeta ) * Jinv[1][2];
				ret[2] = ( x_xi ) * Jinv[2][0] + ( x_eta ) * Jinv[2][1]  + ( x_zeta ) * Jinv[2][2];
			}
			return ;
		case FLUX_FIELD:
			this->getField(GRADIENT_FIELD, p_, ret, true) ;
			ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret) ;
			return ;*/
	}
}

void GeneralizedSpaceTimeViscoElasticElementState::getFieldAtNodes( FieldType f, Vector & ret, int ) 
{
	int totaldof = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	int realdof = parent->spaceDimensions() ;
	int blocks = totaldof / realdof ;
	switch(f)
	{
		case DISPLACEMENT_FIELD:
			for(int j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				for(int i = 0 ; i < realdof ; i++)
					ret[j*realdof+i] = displacements[ j*totaldof + i ] ;
			}
			return ;
		case ENRICHED_DISPLACEMENT_FIELD:
			for(int j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				for(int i = 0 ; i < realdof ; i++)
					ret[j*realdof+i] = enrichedDisplacements[ j*totaldof + i ] ;
			}
			return ;
		case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD:
			ret = displacements ;
			return ;
		case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD:
			ret = enrichedDisplacements ;
			return ;
		case STRAIN_FIELD :

			ElementState::getField(f, parent->getBoundingPoints(), ret, false) ;
			return ;
		case PRINCIPAL_STRAIN_FIELD :
			ElementState::getField(f, parent->getBoundingPoints(), ret, false) ;
			return ;
		case NON_ENRICHED_STRAIN_FIELD :
			ElementState::getField(f, parent->getBoundingPoints(), ret, false) ;
			return ;
		case REAL_STRESS_FIELD:
			ElementState::getField(f,  parent->getBoundingPoints(), ret, false) ;
			return ;
		case PRINCIPAL_REAL_STRESS_FIELD :
			ElementState::getField(f, parent->getBoundingPoints(), ret, false) ;
			return ;
		case NON_ENRICHED_REAL_STRESS_FIELD:
			ElementState::getField(f,  parent->getBoundingPoints(), ret, false) ;
			return ;
		case EFFECTIVE_STRESS_FIELD:
			ElementState::getField(f,  parent->getBoundingPoints(), ret, false) ;
			return ;
		case PRINCIPAL_EFFECTIVE_STRESS_FIELD :
			ElementState::getField(f, parent->getBoundingPoints(), ret, false) ;
			return ;
		case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
			ElementState::getField(f,  parent->getBoundingPoints(), ret, false) ;
			return ;
	}
	ElementState::getField(f, parent->getBoundingPoints(), ret, false) ;
}

void GeneralizedSpaceTimeViscoElasticElementState::getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, int , int)  const 
{
	this->getField(f1, p, ret1, local) ;
	this->getField(f2, p, ret2, local) ;
}

void GeneralizedSpaceTimeViscoElasticElementState::getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, int , int)   
{
	this->getFieldAtNodes(f1, ret1) ;
	this->getFieldAtNodes(f2, ret2) ;
}


