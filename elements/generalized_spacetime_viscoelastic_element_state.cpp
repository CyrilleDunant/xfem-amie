#include "generalized_spacetime_viscoelastic_element_state.h"
#include "../physics/generalized_spacetime_viscoelasticity.h"

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

void GeneralizedSpaceTimeViscoElasticElementState::getAverageField( FieldType f, Vector & ret, int dummy ) 
{
	GaussPointArray gp = parent->getGaussPoints() ;
	ret = 0 ;
	double total = 0 ;
	for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
	{
		Point p_ = gp.gaussPoints[i].first ;
		Vector tmp = ret ;
		getField(f, p_, tmp, dummy) ;
		ret += tmp*gp.gaussPoints[i].second ;
		total += gp.gaussPoints[i].second ;
	}
	ret /= total ;
}

void GeneralizedSpaceTimeViscoElasticElementState::getField( FieldType f, const Point & p, Vector & ret, bool local, int )  const 
{
	VirtualMachine vm ;
	int totaldof = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	int realdof = parent->spaceDimensions() ;
	int blocks = totaldof / realdof ;
	Point p_ = p ;
	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;
	
	GeneralizedSpaceTimeViscoelasticity * visco = dynamic_cast<GeneralizedSpaceTimeViscoelasticity *>(parent->getBehaviour()) ;
	
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

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm.deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm.deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * enrichedDisplacements[j * totaldof + i] ;
						dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
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
				
				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm.ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					x_xi += f_xi * displacements[j * totaldof] ;
					x_eta += f_eta * displacements[j * totaldof] ;
					y_xi += f_xi * displacements[j * totaldof + 1] ;
					y_eta += f_eta * displacements[j * totaldof + 1] ;
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

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm.ddeval( parent->getEnrichmentFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm.ddeval( parent->getEnrichmentFunction( j ), ETA, TIME_VARIABLE,  p_ , 1e-5 ) ;
					for(size_t i = 0 ; i < totaldof ; i++)
					{
						dx[i] += f_xi  * enrichedDisplacements[j * totaldof + i] ;
						dy[i] += f_eta * enrichedDisplacements[j * totaldof + i] ;
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
					+ (Vector) (visco->eta * speeds) ;
			for(size_t i = 0 ; i < 3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL) ; i++)
				ret[i] = stresses[i] ;
			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			return ;
		}
		case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD:
		{
			Vector strains(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			Vector speeds(0., blocks*(3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL))) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_FIELD, p_, strains, true) ;
			this->getField(GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD, p_, speeds, true) ;
			ret = (Vector) (visco->getTensor(p_, parent) * strains) 
			    + (Vector) (visco->eta * speeds) ;
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
					+ (Vector) (visco->eta * speeds) ;
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
			    + (Vector) (visco->eta * speeds) ;
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
					+ (Vector) (visco->eta * speeds) ;
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
			    + (Vector) (visco->eta * speeds) ;
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
					+ (Vector) (visco->eta * speeds) ;
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
			    + (Vector) (visco->eta * speeds) ;
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
			if( strainAtNodes.size() == 0 )
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					strainAtNodes.resize( 3 * parent->getBoundingPoints().size() ) ;
				else
					strainAtNodes.resize( 6 * parent->getBoundingPoints().size() ) ;

				ElementState::getField(f, parent->getBoundingPoints(), strainAtNodes, false) ;
			}
			ret = strainAtNodes ;
			return ;
		case PRINCIPAL_STRAIN_FIELD :
			if( strainAtNodes.size() == 0 )
			{
				strainAtNodes.resize( parent->spaceDimensions() * parent->getBoundingPoints().size() ) ;
				ElementState::getField(f, parent->getBoundingPoints(), strainAtNodes, false) ;
			}
			ret = strainAtNodes ;
			return ;
		case NON_ENRICHED_STRAIN_FIELD :
			if( strainAtNodes.size() == 0 )
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					strainAtNodes.resize( 3 * parent->getBoundingPoints().size() ) ;
				else
					strainAtNodes.resize( 6 * parent->getBoundingPoints().size() ) ;

				ElementState::getField(f, parent->getBoundingPoints(), strainAtNodes, false) ;
			}
			ret = strainAtNodes ;
			return ;
		case REAL_STRESS_FIELD:
			if( stressAtNodes.size() == 0 )
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					stressAtNodes.resize( 3 * parent->getBoundingPoints().size() ) ;
				else
					stressAtNodes.resize( 6 * parent->getBoundingPoints().size() ) ;

				ElementState::getField(f,  parent->getBoundingPoints(), stressAtNodes, false) ;
			}
			ret = stressAtNodes ;
			return ;
		case PRINCIPAL_REAL_STRESS_FIELD :
			if( stressAtNodes.size() == 0 )
			{
				stressAtNodes.resize( parent->spaceDimensions() * parent->getBoundingPoints().size() ) ;
				ElementState::getField(f, parent->getBoundingPoints(), stressAtNodes, false) ;
			}
			ret = stressAtNodes ;
			return ;
		case NON_ENRICHED_REAL_STRESS_FIELD:
			if( stressAtNodes.size() == 0 )
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					stressAtNodes.resize( 3 * parent->getBoundingPoints().size() ) ;
				else
					stressAtNodes.resize( 6 * parent->getBoundingPoints().size() ) ;

				ElementState::getField(f,  parent->getBoundingPoints(), stressAtNodes, false) ;
			}
			ret = stressAtNodes ;
			return ;
		case EFFECTIVE_STRESS_FIELD:
			if( stressAtNodes.size() == 0 )
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					stressAtNodes.resize( 3 * parent->getBoundingPoints().size() ) ;
				else
					stressAtNodes.resize( 6 * parent->getBoundingPoints().size() ) ;

				ElementState::getField(f,  parent->getBoundingPoints(), stressAtNodes, false) ;
			}
			ret = stressAtNodes ;
			return ;
		case PRINCIPAL_EFFECTIVE_STRESS_FIELD :
			if( stressAtNodes.size() == 0 )
			{
				stressAtNodes.resize( parent->spaceDimensions() * parent->getBoundingPoints().size() ) ;
				ElementState::getField(f, parent->getBoundingPoints(), stressAtNodes, false) ;
			}
			ret = stressAtNodes ;
			return ;
		case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
			if( stressAtNodes.size() == 0 )
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					stressAtNodes.resize( 3 * parent->getBoundingPoints().size() ) ;
				else
					stressAtNodes.resize( 6 * parent->getBoundingPoints().size() ) ;

				ElementState::getField(f,  parent->getBoundingPoints(), stressAtNodes, false) ;
			}
			ret = stressAtNodes ;
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
