
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution
//

#include "integrable_entity.h"
#include "../physics/fracturecriteria/fracturecriterion.h"
#include "../physics/kelvinvoight.h"
#include "../solvers/assembly.h"
#include "../features/boundarycondition.h"
#include "../physics/damagemodels/damagemodel.h"
using namespace Mu ;


size_t Mu::fieldTypeElementarySize(FieldType f, SpaceDimensionality dim, size_t blocks) 
{
	switch(f)
	{
		case VON_MISES_STRAIN_FIELD :
		case VON_MISES_REAL_STRESS_FIELD :
		case VON_MISES_EFFECTIVE_STRESS_FIELD :
		case PRINCIPAL_ANGLE_FIELD :
			return 1 ;
	
		case DISPLACEMENT_FIELD :
		case ENRICHED_DISPLACEMENT_FIELD :
		case SPEED_FIELD :
		case FLUX_FIELD:
		case GRADIENT_FIELD:
		case PRINCIPAL_STRAIN_FIELD:
		case PRINCIPAL_EFFECTIVE_STRESS_FIELD:
		case PRINCIPAL_REAL_STRESS_FIELD:
			return dim == SPACE_THREE_DIMENSIONAL ? 3 : 2 ;
			
		case STRAIN_FIELD:
		case STRAIN_RATE_FIELD:
		case EFFECTIVE_STRESS_FIELD:
		case REAL_STRESS_FIELD:
		case NON_ENRICHED_STRAIN_FIELD:
		case NON_ENRICHED_STRAIN_RATE_FIELD:
		case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
		case NON_ENRICHED_REAL_STRESS_FIELD:
			return dim == SPACE_THREE_DIMENSIONAL ? 6 : 3 ;
		  

		case GENERALIZED_VISCOELASTIC_DISPLACEMENT_FIELD :
		case GENERALIZED_VISCOELASTIC_ENRICHED_DISPLACEMENT_FIELD :
		case GENERALIZED_VISCOELASTIC_SPEED_FIELD :
			return dim == SPACE_THREE_DIMENSIONAL ? 3*blocks : 2*blocks ;
			
		case GENERALIZED_VISCOELASTIC_STRAIN_FIELD:
		case GENERALIZED_VISCOELASTIC_STRAIN_RATE_FIELD:
		case GENERALIZED_VISCOELASTIC_EFFECTIVE_STRESS_FIELD:
		case GENERALIZED_VISCOELASTIC_REAL_STRESS_FIELD:
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_FIELD:
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_STRAIN_RATE_FIELD:
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
		case GENERALIZED_VISCOELASTIC_NON_ENRICHED_REAL_STRESS_FIELD:
			return dim == SPACE_THREE_DIMENSIONAL ? 6*blocks : 3*blocks ;
	  
	}
	std::cout << f << " non-existing field !" << std::endl ;
	exit(0) ;
	return 0 ;
}


ElementState * Form::createElementState( IntegrableEntity * e) 
{
	return new ElementState(e) ;
}

IntegrableEntity::IntegrableEntity() : boundaryConditionCache( nullptr ), cachedGps( nullptr )
{
	state = new ElementState( this ) ;
	enrichmentUpdated = false ;
	behaviourUpdated = false ;
}


Function IntegrableEntity::getZTransform() const { return Function("1") ;};

Function IntegrableEntity::getTTransform() const { return Function("1") ;};

Vector Form::getImposedStress(const Point & p, IntegrableEntity * e, int g) const
{
	if(getDamageModel() && getDamageModel()->hasInducedForces())
		return getDamageModel()->getImposedStress(p) ;
	
	return Vector(double(0), getTensor(p, e).numCols()) ;
}

Vector Form::getImposedStrain(const Point & p, IntegrableEntity * e, int g) const
{
	if(getDamageModel() && getDamageModel()->hasInducedForces())
		return getDamageModel()->getImposedStrain(p) ;
	
	return Vector(double(0), getTensor(p).numCols()) ;
}

bool Form::hasInducedForces() const
{
	return false || getDamageModel() && getDamageModel()->hasInducedForces();
} ;

std::vector<BoundaryCondition * > Form::getBoundaryConditions(const ElementState & s, size_t id,  const Function & p_i, const GaussPointArray &gp, const std::valarray<Matrix> &Jinv) const 
{ 
	std::vector<BoundaryCondition * > ret ;
	if(getDamageModel() && getDamageModel()->hasInducedBoundaryConditions())
	{
		return getDamageModel()->getBoundaryConditions(s, id,  p_i, gp, Jinv) ;
	}
		
	return  ret ;
} ;

// void Form::setFractureCriterion(FractureCriterion * frac) 
// {
// 	if(frac)
// 		this->getFractureCriterion() = frac ;
// }


void IntegrableEntity::applyBoundaryCondition( Assembly *a )
{
	if( !getBehaviour())
		return ;
	std::vector<BoundaryCondition *> bcCache ;
	if( getBehaviour()->type != VOID_BEHAVIOUR )
	{
		if(boundaryConditionCache)
		{
			for(size_t i = 0 ; i < boundaryConditionCache->size() ; i++)
			{
				if ((*boundaryConditionCache)[i])
				{
					if( get2DMesh() )
						(*boundaryConditionCache)[i]->apply( a, get2DMesh() ) ;
					else
						(*boundaryConditionCache)[i]->apply( a, get3DMesh() ) ;			  
				}
			}
		}
		
		std::valarray<Matrix> Jinv( getGaussPoints().gaussPoints.size() ) ;
		
		for( size_t i = 0 ; i < getGaussPoints().gaussPoints.size() ;  i++ )
		{
			getInverseJacobianMatrix( getGaussPoints().gaussPoints[i].first, Jinv[i] ) ;
		}

		size_t start = 0 ;
		if(timePlanes() > 1)
		{
			start = getBoundingPoints().size() - getBoundingPoints().size()/timePlanes() ;
		}
		for( size_t i = start ; i < getBoundingPoints().size() ; i++ )
		{
			std::vector<BoundaryCondition *> boundaryConditionCachetmp = getBehaviour()->getBoundaryConditions( getState(), getBoundingPoint( i ).id,  getShapeFunction( i ), getGaussPoints(), Jinv ) ;
			for(size_t j = 0 ; j < boundaryConditionCachetmp.size() ; j++)
			{
				if( get2DMesh() )
					boundaryConditionCachetmp[j]->apply( a, get2DMesh() ) ;
				else
					boundaryConditionCachetmp[j]->apply( a, get3DMesh() ) ;
				delete boundaryConditionCachetmp[j] ;
			}
		}

// 		start = 0 ;
// 		if(timePlanes() > 1)
// 		{
// 			start = getEnrichmentFunctions().size() - getEnrichmentFunctions().size()/timePlanes() ;
// 		}
// 		for( size_t i = start ; i < getEnrichmentFunctions().size() ; i++ )
// 		{
// //			std::cout << getGaussPoints().gaussPoints.size() << std::endl ;
// 			std::vector<BoundaryCondition *> boundaryConditionCachetmp = getBehaviour()->getBoundaryConditions( getState(), getEnrichmentFunction( i ).getDofID(),  getEnrichmentFunction( i ), getGaussPoints(), Jinv ) ;
// 			for(size_t j = 0 ; j < boundaryConditionCachetmp.size() ; j++)
// 			{
// 				if( get2DMesh() )
// 					boundaryConditionCachetmp[j]->apply( a, get2DMesh() ) ;
// 				else
// 					boundaryConditionCachetmp[j]->apply( a, get3DMesh() ) ;
// 				delete boundaryConditionCachetmp[j] ;
// 			}
// 		}
		
	}
}

IntegrableEntity::~IntegrableEntity()
{
	if( boundaryConditionCache )
	{
		for( size_t i = 0 ; i < boundaryConditionCache->size() ; i++ )
			delete( *boundaryConditionCache )[i] ;
	}
	delete state ;
	delete boundaryConditionCache ;
	delete cachedGps ;
}


void IntegrableEntity::addBoundaryCondition( BoundaryCondition * bc) 
{
	if( !boundaryConditionCache )
	{
		boundaryConditionCache = new std::vector<BoundaryCondition *>(0, nullptr) ;
	}
	boundaryConditionCache->push_back(bc);
}

void IntegrableEntity::clearBoundaryConditions() 
{
	if(boundaryConditionCache)
	{
		for( size_t i = 0 ; i < boundaryConditionCache->size() ; i++ )
			delete( *boundaryConditionCache )[i] ;
		delete boundaryConditionCache ;
	}
	boundaryConditionCache = nullptr ;
}


void IntegrableEntity::setState( ElementState * s)
{
	delete state ;
	state = s ;
}

const ElementState &IntegrableEntity::getState() const
{
	return *state ;
}

ElementState &IntegrableEntity::getState()
{
	return *state ;
}


const Vector &ElementState::getDisplacements() const
{
	return this->displacements ;
}

Vector &ElementState::getDisplacements()
{
	return this->displacements ;
}

const Vector &ElementState::getEnrichedDisplacements() const
{
	return this->enrichedDisplacements ;
}

Vector &ElementState::getEnrichedDisplacements()
{
	return this->enrichedDisplacements ;
}

const Vector &ElementState::getPreviousDisplacements() const
{
	return this->previousDisplacements ;
}

Vector &ElementState::getPreviousDisplacements()
{
	return this->previousDisplacements ;
}

const Vector &ElementState::getPreviousEnrichedDisplacements() const
{
	return this->previousEnrichedDisplacements ;
}

Vector &ElementState::getPreviousEnrichedDisplacements()
{
	return this->previousEnrichedDisplacements ;
}



ElementState &ElementState::operator =( const ElementState &s )
{
	displacements.resize( s.getDisplacements().size() ) ;
	displacements = s.getDisplacements() ;
	enrichedDisplacements.resize( s.getEnrichedDisplacements().size() ) ;
	enrichedDisplacements  = s.getEnrichedDisplacements();

	previousDisplacements.resize( s.getPreviousDisplacements().size() ) ;
	previousDisplacements = s.getPreviousDisplacements() ;
	previousEnrichedDisplacements.resize( s.getPreviousEnrichedDisplacements().size() ) ;
	previousEnrichedDisplacements = s. getPreviousEnrichedDisplacements();


	buffer.resize( s.getBuffer().size() ) ;
	buffer = s.getBuffer() ;

	timePos = s.getTime();
	previousTimePos = s.getTime() - s.getDeltaTime();

	parent = s.getParent();
	return *this ;
}

ElementState::ElementState( const ElementState &s )
{
	strainAtGaussPoints.resize(0) ;
	stressAtGaussPoints.resize(0) ;
	effectivePStressAtGaussPoints.resize(0) ;
	pstrainAtGaussPoints.resize(0) ;
	pstressAtGaussPoints.resize(0) ;
	
	displacements.resize( s.getDisplacements().size() ) ;
	displacements = s.getDisplacements() ;
	enrichedDisplacements.resize( s.getEnrichedDisplacements().size() ) ;
	enrichedDisplacements  = s.getEnrichedDisplacements();

	previousDisplacements.resize( s.getPreviousDisplacements().size() ) ;
	previousDisplacements = s.getPreviousDisplacements() ;
	previousEnrichedDisplacements.resize( s.getPreviousEnrichedDisplacements().size() ) ;
	previousEnrichedDisplacements = s. getPreviousEnrichedDisplacements();


	buffer.resize( s.getBuffer().size() ) ;
	buffer = s.getBuffer() ;

	timePos = s.getTime();
	previousTimePos = s.getTime() - s.getDeltaTime();

	parent = s.getParent();
}

Vector ElementState::getAverageDisplacement(VirtualMachine * vm) const
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	FunctionMatrix disps = getDisplacementFunction() ;
	Matrix ret = vm->ieval( disps, getParent() ) ;
	Vector v( ret.numRows() ) ;

	if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		for( size_t i = 0 ;  i < v.size() ; i++ )
			v[i] = ret[i][0] / parent->area() ;
	}
	else
	{
		for( size_t i = 0 ;  i < v.size() ; i++ )
			v[i] = ret[i][0] / parent->volume() ;
	}
	if (cleanup) delete vm ;
	return v ;
}

ElementState::ElementState( IntegrableEntity *s )
{

	parent = s ;
	this->timePos = 0 ;
	this->previousTimePos = 0 ;
// 	size_t ndof = 2 ;
// // 	if(s->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
// // 		ndof = 3 ;
// 	this->displacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->previousDisplacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->previousPreviousDisplacements.resize(s->getBoundingPoints()->size()*ndof) ;
// 	this->enrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
// 	this->previousEnrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
// 	this->previousPreviousEnrichedDisplacements.resize(s->getEnrichmentFunctions().size()*ndof) ;
}

Matrix Mu::makeStressOrStrainMatrix(const Vector & stressOrStrain) 
{
	if(stressOrStrain.size() == 3)
	{
		Matrix ret2(2,2) ;
		ret2[0][0] = stressOrStrain[0] ;
		ret2[1][1] = stressOrStrain[1] ;
		ret2[1][0] = ret2[0][1] = stressOrStrain[2] * .5 ;
		return ret2 ;
	}
	else if(stressOrStrain.size() == 6)
	{
		Matrix ret3(3,3) ;
		ret3[0][0] = stressOrStrain[0] ;
		ret3[1][1] = stressOrStrain[1] ;
		ret3[2][2] = stressOrStrain[2] ;
		ret3[2][0] = ret3[0][2] = stressOrStrain[3] * .5 ;
		ret3[2][1] = ret3[1][2] = stressOrStrain[4] * .5 ;
		ret3[1][0] = ret3[0][1] = stressOrStrain[5] * .5 ;
		return ret3 ;
	}
	return Matrix(2 + (stressOrStrain.size()==6), 2 + (stressOrStrain.size()==6)) ;
}

bool isStrainField(FieldType f)
{
	return f == STRAIN_FIELD || f == NON_ENRICHED_STRAIN_FIELD/* || f == PRINCIPAL_STRAIN_FIELD*/ ;
}

bool isStressField(FieldType f)
{
	return f == REAL_STRESS_FIELD || f == NON_ENRICHED_REAL_STRESS_FIELD /*|| f == PRINCIPAL_REAL_STRESS_FIELD */
		|| f == EFFECTIVE_STRESS_FIELD || f == NON_ENRICHED_EFFECTIVE_STRESS_FIELD /*|| f == PRINCIPAL_EFFECTIVE_STRESS_FIELD*/  ;
}

bool isRealStressField(FieldType f)
{
	return f == REAL_STRESS_FIELD || f == NON_ENRICHED_REAL_STRESS_FIELD || f == PRINCIPAL_REAL_STRESS_FIELD  ;
}

bool isEffectiveStressField(FieldType f)
{
	return f == EFFECTIVE_STRESS_FIELD || f == NON_ENRICHED_EFFECTIVE_STRESS_FIELD || f == PRINCIPAL_EFFECTIVE_STRESS_FIELD  ;
}

Vector Mu::toPrincipal(const Vector & stressOrStrain)
{
	Vector ret(0., 2+(stressOrStrain.size() == 6)) ;
	if(ret.size() == 2)
	{
		ret[0] = ( stressOrStrain[0] + stressOrStrain[1] ) * .5 + 
			    sqrt( 0.25 *( stressOrStrain[0] - stressOrStrain[1] ) * ( stressOrStrain[0] - stressOrStrain[1] ) + 
			    ( stressOrStrain[2] * stressOrStrain[2] ) ) ;
		ret[1] = ( stressOrStrain[0] + stressOrStrain[1] ) * .5 - 
			    sqrt( 0.25 *( stressOrStrain[0] - stressOrStrain[1] ) * ( stressOrStrain[0] - stressOrStrain[1] ) + 
			    ( stressOrStrain[2] * stressOrStrain[2] ) ) ;
	}
	else if(ret.size() == 3)
	{
		Matrix mat = Mu::makeStressOrStrainMatrix(stressOrStrain) ;
		Matrix I( 3, 3 ) ;
		I[0][0] = 1 ;
		I[1][1] = 1 ;
		I[2][2] = 1 ;
		double m = ( mat[0][0] + mat[1][1] + mat[2][2] ) / 3. ;
		Matrix Am = mat - I * m ;
		double q = det( Am ) / 2. ;
		double r = std::inner_product( &Am.array()[0], &Am.array()[9], &Am.array()[0],  double( 0. ) ) / 6. ;
		double phi = atan2( sqrt( r * r * r - q * q ), q ) / 3. ;

		if( r * r * r - q * q < 1e-12 )
			phi = atan( 0 ) / 3. ;

		if( phi < 0 )
			phi += M_PI ;
		
		ret[0] = m + 2.*sqrt( r ) * cos( phi ) ;
		ret[1] = m - sqrt( r ) * ( cos( phi ) + sqrt( 3. ) * sin( phi ) ) ;
		ret[2] = m - sqrt( r ) * ( cos( phi ) - sqrt( 3. ) * sin( phi ) ) ;
	}
	return ret ;
}

void ElementState::getExternalField( Vector & nodalValues, int externaldofs, const Point & p, Vector & ret, bool local, VirtualMachine * vm) const 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	Point p_ = p ;
	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;

	ret = 0. ;
	
	for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
	{
		double f = vm->eval( parent->getShapeFunction( j ), p_) ;
		for(size_t k = 0 ; k < externaldofs ; k++)
			ret[k] += f * nodalValues[ j*externaldofs + k ] ;
	}
	if ( cleanup ) delete vm ;
}

void ElementState::getExternalFieldAtGaussPoints( Vector & nodalValues, int externaldofs, std::vector<Vector> & ret, VirtualMachine * vm) const 
{
	for(size_t p = 0 ; p < parent->getGaussPoints().gaussPoints.size() ; p++)
	{
		getExternalField( nodalValues, externaldofs, parent->getGaussPoints().gaussPoints[p].first, ret[p], true, vm) ;
	}
}

void ElementState::getField( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm, int )  const 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	int n = 0 ;
	Point p_ = p ;
	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;
	
	switch(f)
	{
		case DISPLACEMENT_FIELD:
			n =  parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f =  vm ->eval( parent->getShapeFunction( j ) , p_) ;
				for(size_t k = 0 ; k < n ; k++)
					ret[k] += f * displacements[j*n+k] ;
			}
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm ->eval( parent->getEnrichmentFunction( j ) , p_) ;
				for(size_t k = 0 ; k < n ; k++)
					ret[k] += f * enrichedDisplacements[j*n+k] ;
			}
			if ( cleanup ) delete vm ;
			return ;
		case ENRICHED_DISPLACEMENT_FIELD:
			n =  parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm ->eval( parent->getEnrichmentFunction( j ) , p_) ;
				for(size_t k = 0 ; k < n ; k++)
					ret[k] += f * enrichedDisplacements[j*n+k] ;
			}
			if ( cleanup ) delete vm ;
			return ;
		case STRAIN_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
			{
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;
				
				Vector dx(0., parent->getBehaviour()->getNumberOfDegreesOfFreedom()) ;
				Vector dy(0., parent->getBehaviour()->getNumberOfDegreesOfFreedom()) ;
				if(displacements.size()/2 == parent->getBoundingPoints().size())
				{
					for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
					{
						double f_xi = vm ->deval( parent->getShapeFunction( j ), XI, p_ ) ;
						double f_eta = vm ->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
						x_xi += f_xi * displacements[j * 2] ;
						x_eta += f_eta * displacements[j * 2] ;
						y_xi += f_xi * displacements[j * 2 + 1] ;
						y_eta += f_eta * displacements[j * 2 + 1] ;
					}
				}
				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm ->deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;

					x_xi += f_xi * enrichedDisplacements[j * 2] ;
					x_eta += f_eta * enrichedDisplacements[j * 2] ;
					y_xi += f_xi * enrichedDisplacements[j * 2 + 1] ;
					y_eta += f_eta * enrichedDisplacements[j * 2 + 1] ;
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
				ret[2] = ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3 )
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
					double f_xi = vm ->deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm ->deval( parent->getShapeFunction( j ), ZETA, p_ ) ;
					double x = displacements[j * 3] ;
					double y = displacements[j * 3 + 1] ;
					double z = displacements[j * 3 + 2] ;

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
					double f_xi = vm ->deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					double f_zeta = vm ->deval( parent->getEnrichmentFunction( j ), ZETA, p_ ) ;
					double x = enrichedDisplacements[j * 3] ;
					double y = enrichedDisplacements[j * 3 + 1] ;
					double z = enrichedDisplacements[j * 3 + 2] ;

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

				ret[3] = ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] =  ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] =  ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
			}
			if ( cleanup ) delete vm ;
			return ;
		case PRINCIPAL_STRAIN_FIELD:
		{
			Vector strains(0.,3) ;
			if(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
				strains.resize(6) ;
			this->getField(STRAIN_FIELD, p_, strains, true,vm) ;

			ret = toPrincipal(strains) ;
			if ( cleanup ) delete vm ;
			return ;
		}
		case NON_ENRICHED_STRAIN_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2 )
			{
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;

				if(displacements.size()/2 == parent->getBoundingPoints().size())
				{
					for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
					{
						double f_xi = vm ->deval( parent->getShapeFunction( j ), XI, p_ ) ;
						double f_eta = vm ->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
						x_xi += f_xi * displacements[j * 2] ;
						x_eta += f_eta * displacements[j * 2] ;
						y_xi += f_xi * displacements[j * 2 + 1] ;
						y_eta += f_eta * displacements[j * 2 + 1] ;
					}
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
				ret[2] =  ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
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
					double f_xi = vm ->deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm ->deval( parent->getShapeFunction( j ), ZETA, p_ ) ;
					double x = displacements[j * 3] ;
					double y = displacements[j * 3 + 1] ;
					double z = displacements[j * 3 + 2] ;

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


				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1]  + ( x_zeta ) * Jinv[0][2];
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
				ret[2] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];

				ret[3] =  ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] =  ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] =  ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
			}
			if ( cleanup ) delete vm ;
			return ;
		case VON_MISES_STRAIN_FIELD:
		{
			Vector eps(0., (size_t) parent->spaceDimensions()) ;
			getField( PRINCIPAL_STRAIN_FIELD, p_, eps, true,vm ) ;
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				ret[0] = ( 2. / 3. * ( eps[0] * eps[0] + eps[1] * eps[1] ) ) ;
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL )
			{
				ret[0] = sqrt( 2. / 3. * ( eps[0] * eps[0] + eps[1] * eps[1] + eps[2] * eps[2] ) ) ;
			}
			if ( cleanup ) delete vm ;
			return ;
		}
		case REAL_STRESS_FIELD:
			
			getField(STRAIN_FIELD, p_, ret, true, vm) ;
			if(parent->getBehaviour()->getTensor(p_, parent).numCols() != ret.size())
			{
				ret = 0 ;
				if ( cleanup ) delete vm ;
				return ;
			}
			ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			if ( cleanup ) delete vm ;
			return ;
		case PRINCIPAL_REAL_STRESS_FIELD:
		{
			Vector stress(0.,3+3*(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
			getField(REAL_STRESS_FIELD, p_, stress, true,vm) ;
			ret = toPrincipal(stress) ;
			if ( cleanup ) delete vm ;
			return ;
		}
		case NON_ENRICHED_REAL_STRESS_FIELD:
			getField(NON_ENRICHED_STRAIN_FIELD, p_, ret, true,vm) ;
			if(parent->getBehaviour()->getTensor(p_, parent).numCols() != ret.size())
			{
				ret = 0 ;
				if ( cleanup ) delete vm ;
				return ;
			}
			ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			if ( cleanup ) delete vm ;
			return ;
		case VON_MISES_REAL_STRESS_FIELD:
			if( parent->getOrder() == LINEAR )
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector sigma(0., 2) ;
					Point c(1./3., 1./3.) ;
					getField(PRINCIPAL_REAL_STRESS_FIELD, c, sigma, true, vm) ;
					ret[0] = sqrt( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
					if ( cleanup ) delete vm ;
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
					getField( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false, vm ) ;
					sigma[0] = sqrt( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
					sigma[1] = sqrt( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
					sigma[2] = sqrt( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
					sigma[3] = sqrt( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

					ret[0] = std::max( std::max( std::max( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
					if ( cleanup ) delete vm ;
					return ;
				}
				return ;
			}
			else
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector principalStresses(0., parent->getBoundingPoints().size()*2) ;
					getField(PRINCIPAL_REAL_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false, vm) ;
					double maxS = 0 ;

					for( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
					{
						maxS = std::max( maxS,
								sqrt( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
					}

					ret[0] = maxS ;
					if ( cleanup ) delete vm ;
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
					this->getField( PRINCIPAL_REAL_STRESS_FIELD, pts, sigma, false, vm) ;
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
			getField(STRAIN_FIELD, p_, ret, true,vm, 0) ;
			ret = (Vector) (parent->getBehaviour()->param * ret) - getParent()->getBehaviour()->getImposedStrain(p_, parent)*parent->getBehaviour()->param ;
			if ( cleanup ) delete vm ;
			return ;
		case PRINCIPAL_EFFECTIVE_STRESS_FIELD:
		{
			Vector stress(0.,3+3*(parent->spaceDimensions()== SPACE_THREE_DIMENSIONAL)) ;
			getField(EFFECTIVE_STRESS_FIELD, p_, stress, true,vm) ;
			ret = toPrincipal(stress) ; 
			if ( cleanup ) delete vm ;
			return ;
		}
		case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
			getField(NON_ENRICHED_STRAIN_FIELD, p_, ret, true,vm) ;
			ret = (Vector) (parent->getBehaviour()->param * ret) - getParent()->getBehaviour()->getImposedStrain(p_, parent)*parent->getBehaviour()->param ;
			if ( cleanup ) delete vm ;
			return ;
		case VON_MISES_EFFECTIVE_STRESS_FIELD:
			if( parent->getOrder() == LINEAR )
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector sigma(0., 2) ;
					Point c(1./3., 1./3.) ;
					getField(PRINCIPAL_EFFECTIVE_STRESS_FIELD, c, sigma, true,vm) ;
					ret[0] = sqrt( ( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + sigma[0] * sigma[0] + sigma[1] * sigma[1] ) / 2. ) ;
					if ( cleanup ) delete vm ;
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
					getField( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false, vm) ;
					sigma[0] = sqrt( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6. ;
					sigma[1] = sqrt( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6. ;
					sigma[2] = sqrt( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6. ;
					sigma[3] = sqrt( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6. ;

					ret[0] = std::max( std::max( std::max( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
					if ( cleanup ) delete vm ;
					return ;
				}
				if ( cleanup ) delete vm ;
				return ;
			}
			else
			{
				if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				{
					Vector principalStresses(0., parent->getBoundingPoints().size()*2) ;
					getField(PRINCIPAL_EFFECTIVE_STRESS_FIELD, parent->getBoundingPoints(), principalStresses, false, vm) ;
					double maxS = 0 ;

					for( size_t i = 0 ; i < principalStresses.size() / 2 ; i++ )
					{
						maxS = std::max( maxS,
								sqrt( ( ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) * ( principalStresses[i * 2 + 0] - principalStresses[i * 2 + 1] ) + principalStresses[i * 2 + 0] * principalStresses[i * 2 + 0] + principalStresses[i * 2 + 1] * principalStresses[i * 2 + 1] ) / 2. ) ) ;
					}

					ret[0] = maxS ;
					if ( cleanup ) delete vm ;
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
					getField( PRINCIPAL_EFFECTIVE_STRESS_FIELD, pts, sigma, false, vm) ;
					sigma[0] = sqrt( ( sigma[0] - sigma[1] ) * ( sigma[0] - sigma[1] ) + ( sigma[0] - sigma[2] ) * ( sigma[0] - sigma[2] ) + ( sigma[1] - sigma[2] ) * ( sigma[1] - sigma[2] ) ) / 6 ;
					sigma[1] = sqrt( ( sigma[6] - sigma[7] ) * ( sigma[6] - sigma[7] ) + ( sigma[6] - sigma[8] ) * ( sigma[6] - sigma[8] ) + ( sigma[7] - sigma[8] ) * ( sigma[7] - sigma[8] ) ) / 6 ;
					sigma[2] = sqrt( ( sigma[12] - sigma[13] ) * ( sigma[12] - sigma[13] ) + ( sigma[12] - sigma[14] ) * ( sigma[12] - sigma[14] ) + ( sigma[13] - sigma[14] ) * ( sigma[13] - sigma[14] ) ) / 6 ;
					sigma[3] = sqrt( ( sigma[18] - sigma[19] ) * ( sigma[18] - sigma[19] ) + ( sigma[18] - sigma[20] ) * ( sigma[18] - sigma[20] ) + ( sigma[19] - sigma[20] ) * ( sigma[19] - sigma[20] ) ) / 6 ;

					ret[0] = std::max( std::max( std::max( sigma[0], sigma[1] ), sigma[2] ), sigma[3] ) ;
					if ( cleanup ) delete vm ;
					return ;
				}
			}
			if ( cleanup ) delete vm ;
			return ;
		case PRINCIPAL_ANGLE_FIELD:
		{
			Vector strains(0., 3+3*(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
			getField(STRAIN_FIELD,  p_, strains, true) ;
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
				ret[0] =  0.5 * atan2( strains[2], strains[0] - strains[1] ) ;
			else
			{
				ret[0] =  0.5 * atan2(strains[3] , strains[0] - strains[1] ) ;
				ret[1] =  0.5 * atan2(strains[4] , strains[0] - strains[2] ) ;
				ret[2] =  0.5 * atan2(strains[5] , strains[1] - strains[2] ) ;
			}
			if ( cleanup ) delete vm ;
			return ;
		}
		case GRADIENT_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 1 )
			{
				double x_xi = 0;
				double x_eta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm ->deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					x_xi += f_xi * displacements[j] ;
					x_eta += f_eta * displacements[j] ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size(); j++ )
				{
					double f_xi = vm ->deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
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
					double f_xi = vm ->deval( parent->getShapeFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getShapeFunction( j ), ETA, p_ ) ;
					double f_zeta = vm ->deval( parent->getShapeFunction( j ), ZETA, p_ ) ;
					double x = displacements[j] ;

					x_xi   += f_xi   * x ;
					x_eta  += f_eta  * x ;
					x_zeta += f_zeta * x ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
				{
					double f_xi = vm ->deval( parent->getEnrichmentFunction( j ), XI, p_ ) ;
					double f_eta = vm ->deval( parent->getEnrichmentFunction( j ), ETA, p_ ) ;
					double f_zeta = vm ->deval( parent->getEnrichmentFunction( j ), ZETA, p_ ) ;
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
			if ( cleanup ) delete vm ;
			return ;
		case FLUX_FIELD:
			getField(GRADIENT_FIELD, p_, ret, true, vm) ;
			ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret) ;
			if ( cleanup ) delete vm ;
			return ;
	}
}

void ElementState::getField( FieldType f, const PointArray & p, Vector & ret, bool local, VirtualMachine * vm, int )  const 
{
	bool cleanup = !vm ;
	if( !vm ) vm = new VirtualMachine() ;
	Vector buffer(0., ret.size()/p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		getField(f, *p[i], buffer, local, vm) ;
		for(size_t j = buffer.size()*i ; j < buffer.size()*(i+1) ; j++)
			ret[j] = buffer[j - buffer.size()*i] ;
	}
	if(cleanup) delete vm ;
}

void ElementState::getField( FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, VirtualMachine * vm, int )  const 
{
	bool cleanup = !vm ;
	if( !vm ) vm = new VirtualMachine() ;
	Vector buffer(0., ret.size()/p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		getField(f, p[i].first, buffer, local, vm) ;
		for(size_t j = buffer.size()*i ; j < buffer.size()*(i+1) ; j++)
			ret[j] = buffer[j - buffer.size()*i] ;
	} 
	if(cleanup) delete vm ;
}

void ElementState::getFieldAtGaussPoint( FieldType f, size_t p, Vector & ret, VirtualMachine * vm, int i) 
{
	bool cleanup = !vm ;
	if( !vm ) vm = new VirtualMachine() ;
	Point p_ = parent->getGaussPoints().gaussPoints[p].first ;
	getField(f, p_, ret, true, vm, i) ;
	if(cleanup) delete vm ;
}

void ElementState::getAverageField( FieldType f, Vector & ret, VirtualMachine * vm, int dummy, double t) 
{
	bool cleanup = !vm ;
	if( !vm ) vm = new VirtualMachine() ;
	GaussPointArray gp = parent->getGaussPoints() ;
	ret = 0 ;
	double total = 0 ;
	
	switch(f)
	{
		case STRAIN_FIELD :
			if( strainAtGaussPoints.size() == 0)
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					strainAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				else
					strainAtGaussPoints.resize( 6*gp.gaussPoints.size() ) ;

				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(strainAtGaussPoints.size()/gp.gaussPoints.size()) ;
					getField(f, gp.gaussPoints[i].first, tmp, true,vm, i) ;
					for(size_t j = 0 ; j < strainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						strainAtGaussPoints[i*strainAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			else
			{
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(strainAtGaussPoints.size()/gp.gaussPoints.size()) ;
					for(size_t j = 0 ; j < strainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						tmp[j] = strainAtGaussPoints[i*strainAtGaussPoints.size()/gp.gaussPoints.size()+j];
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
		case PRINCIPAL_STRAIN_FIELD :
			if( pstrainAtGaussPoints.size() == 0)
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					pstrainAtGaussPoints.resize( 2*gp.gaussPoints.size() ) ;
				else
					pstrainAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(pstrainAtGaussPoints.size()/gp.gaussPoints.size()) ;
					getField(f, gp.gaussPoints[i].first, tmp, true, vm, i) ;
					for(size_t j = 0 ; j < pstrainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						pstrainAtGaussPoints[i*pstrainAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
					
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			else
			{
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(pstrainAtGaussPoints.size()/gp.gaussPoints.size()) ;
					for(size_t j = 0 ; j < pstrainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						tmp[j] = pstrainAtGaussPoints[i*pstrainAtGaussPoints.size()/gp.gaussPoints.size()+j];
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
		case REAL_STRESS_FIELD:
			if( stressAtGaussPoints.size() == 0)
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					stressAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				else
					stressAtGaussPoints.resize( 6*gp.gaussPoints.size() ) ;

				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(stressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					getField(f, gp.gaussPoints[i].first, tmp, true, vm, i) ;
					for(size_t j = 0 ; j < stressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						stressAtGaussPoints[i*stressAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			else
			{
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(stressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					for(size_t j = 0 ; j < stressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						tmp[j] = stressAtGaussPoints[i*stressAtGaussPoints.size()/gp.gaussPoints.size()+j];
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
		case PRINCIPAL_REAL_STRESS_FIELD :
			if( pstressAtGaussPoints.size() == 0)
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					pstressAtGaussPoints.resize( 2*gp.gaussPoints.size() ) ;
				else
					pstressAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(pstressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					getField(f, gp.gaussPoints[i].first, tmp, true, vm,i) ;
					for(size_t j = 0 ; j < pstressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						pstressAtGaussPoints[i*pstressAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			else
			{
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Vector tmp(pstressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					for(size_t j = 0 ; j < pstressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						tmp[j] = pstressAtGaussPoints[i*pstressAtGaussPoints.size()/gp.gaussPoints.size()+j];
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
		case EFFECTIVE_STRESS_FIELD:
			if( stressAtGaussPoints.size() == 0)
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					stressAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				else
					stressAtGaussPoints.resize( 6*gp.gaussPoints.size() ) ;

				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Point p_ = gp.gaussPoints[i].first ;
					Vector tmp(stressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					getField(f, p_, tmp, true, vm, i) ;
					for(size_t j = 0 ; j < stressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						stressAtGaussPoints[i*stressAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			else
			{
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Point p_ = gp.gaussPoints[i].first ;
					Vector tmp(stressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					for(size_t j = 0 ; j < stressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						tmp[j] = stressAtGaussPoints[i*stressAtGaussPoints.size()/gp.gaussPoints.size()+j];
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
		case PRINCIPAL_EFFECTIVE_STRESS_FIELD :
			if( pstressAtGaussPoints.size() == 0)
			{
				if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
					pstressAtGaussPoints.resize( 2*gp.gaussPoints.size() ) ;
				else
					pstressAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Point p_ = gp.gaussPoints[i].first ;
					Vector tmp(pstressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					getField(f, p_, tmp, true, vm, i) ;
					for(size_t j = 0 ; j < pstressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						pstressAtGaussPoints[i*pstressAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			else
			{
				for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
				{
					Point p_ = gp.gaussPoints[i].first ;
					Vector tmp(pstressAtGaussPoints.size()/gp.gaussPoints.size()) ;
					for(size_t j = 0 ; j < pstressAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
					{
						tmp[j] = pstressAtGaussPoints[i*pstressAtGaussPoints.size()/gp.gaussPoints.size()+j];
					}
					ret += tmp*gp.gaussPoints[i].second ;
					total += gp.gaussPoints[i].second ;
				}
				if(total > POINT_TOLERANCE_2D)
					ret /= total ;
				else
					ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
		default :
			for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			{
				Point p_ = gp.gaussPoints[i].first ;
				Vector tmp = ret ;
				getField(f, p_, tmp,true,  vm, dummy) ;
				ret += tmp*gp.gaussPoints[i].second ;
				total += gp.gaussPoints[i].second ;
			}
			if(total > POINT_TOLERANCE_2D)
				ret /= total ;
			else
				ret = 0 ;
			if(cleanup) delete vm ;
			return ;
	}
}

void ElementState::getAverageField( FieldType f, FieldType f_, Vector & ret, Vector & ret_, VirtualMachine * vm,  int dummy, double t) 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ; 
	GaussPointArray gp = parent->getGaussPoints() ;
	ret = 0 ;
	ret_ = 0 ;
	double total = 0 ;
	
	if(f == STRAIN_FIELD && (f_ == EFFECTIVE_STRESS_FIELD || f_ == REAL_STRESS_FIELD))
	{
		if( strainAtGaussPoints.size() == 0)
		{
			if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{
				strainAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				stressAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
			}
			else
			{
				strainAtGaussPoints.resize( 6*gp.gaussPoints.size() ) ;
				stressAtGaussPoints.resize( 6*gp.gaussPoints.size() ) ;
			}
			
			for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			{
				Vector tmp(strainAtGaussPoints.size()/gp.gaussPoints.size()) ;
				Vector tmp_(stressAtGaussPoints.size()/gp.gaussPoints.size()) ;
				getField(f,f_, gp.gaussPoints[i].first, tmp, tmp_, true,vm, i) ;
				
				for(size_t j = 0 ; j < strainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
				{
					strainAtGaussPoints[i*strainAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					stressAtGaussPoints[i*strainAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp_[j] ;
				}

				ret += tmp*gp.gaussPoints[i].second ;
				ret_ += tmp_*gp.gaussPoints[i].second ;
				total += gp.gaussPoints[i].second ;
			}
			if(total > POINT_TOLERANCE_2D)
			{
				ret_ /= total ;
				ret /= total ;
			}
			else
			{
				ret = 0 ;
				ret_ = 0 ;
			}
			if(cleanup) delete vm ;
			return ;
		}
		else
		{
			for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			{
				Vector tmp(strainAtGaussPoints.size()/gp.gaussPoints.size()) ;
				Vector tmp_(stressAtGaussPoints.size()/gp.gaussPoints.size()) ;
				for(size_t j = 0 ; j < strainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
				{
					tmp[j] = strainAtGaussPoints[i*strainAtGaussPoints.size()/gp.gaussPoints.size()+j];
					tmp_[j] = stressAtGaussPoints[i*stressAtGaussPoints.size()/gp.gaussPoints.size()+j];
				}
				ret += tmp*gp.gaussPoints[i].second ;
				ret_ += tmp_*gp.gaussPoints[i].second ;
				total += gp.gaussPoints[i].second ;
			}
			if(total > POINT_TOLERANCE_2D)
			{
				ret_ /= total ;
				ret /= total ;
			}
			else
			{
				ret = 0 ;
				ret_ = 0 ;
			}
			if(cleanup) delete vm ;
			return ;
		}
	}
	if(f == PRINCIPAL_STRAIN_FIELD && (f_ == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f == PRINCIPAL_REAL_STRESS_FIELD))
	{
		if( pstrainAtGaussPoints.size() == 0)
		{
			if(parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
			{
				pstrainAtGaussPoints.resize( 2*gp.gaussPoints.size() ) ;
				pstressAtGaussPoints.resize( 2*gp.gaussPoints.size() ) ;
			}
			else
			{
				pstrainAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
				pstressAtGaussPoints.resize( 3*gp.gaussPoints.size() ) ;
			}
			for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			{
				Vector tmp(pstrainAtGaussPoints.size()/gp.gaussPoints.size()) ;
				Vector tmp_(pstressAtGaussPoints.size()/gp.gaussPoints.size()) ;
				getField(f,f_, gp.gaussPoints[i].first, tmp,tmp_, true, vm, i) ;
				for(size_t j = 0 ; j < pstrainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
				{
					pstrainAtGaussPoints[i*pstrainAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp[j] ;
					pstressAtGaussPoints[i*pstressAtGaussPoints.size()/gp.gaussPoints.size()+j] = tmp_[j] ;
				}
				ret += tmp*gp.gaussPoints[i].second ;
				ret_ += tmp_*gp.gaussPoints[i].second ;
				total += gp.gaussPoints[i].second ;
				
			}
			if(total > POINT_TOLERANCE_2D)
			{
				ret_ /= total ;
				ret /= total ;
			}
			else
			{
				ret = 0 ;
				ret_ = 0 ;
			}
			if(cleanup) delete vm ;
			return ;
		}
		else
		{
			for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
			{
				Vector tmp(pstrainAtGaussPoints.size()/gp.gaussPoints.size()) ;
				Vector tmp_(pstressAtGaussPoints.size()/gp.gaussPoints.size()) ;
				for(size_t j = 0 ; j < pstrainAtGaussPoints.size()/gp.gaussPoints.size() ; j++)
				{
					tmp[j] = pstrainAtGaussPoints[i*pstrainAtGaussPoints.size()/gp.gaussPoints.size()+j];
					tmp_[j] = pstressAtGaussPoints[i*pstressAtGaussPoints.size()/gp.gaussPoints.size()+j];
				}
				ret += tmp*gp.gaussPoints[i].second ;
				ret_ += tmp_*gp.gaussPoints[i].second ;
				total += gp.gaussPoints[i].second ;
			}
			if(total > POINT_TOLERANCE_2D)
			{
				ret_ /= total ;
				ret /= total ;
			}
			else
			{
				ret = 0 ;
				ret_ = 0 ;
			}
			return ;
		}
	}

	getAverageField(f, ret,vm, dummy) ;
	getAverageField(f_, ret_,vm, dummy);
	if(cleanup) delete vm ;
	
}


void ElementState::getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int)  const 
{
	bool cleanup = !vm ;
	if (!vm) vm = new VirtualMachine() ;
	Point p_ = p ;
	if(!local)
		p_ = parent->inLocalCoordinates(p) ;
	if(isStrainField(f1) && isStressField(f2))
	{
		this->getField(f1, p, ret1, local, vm) ;
		if(parent->getBehaviour()->getTensor(p_, parent).numCols() != ret2.size())
		{
			ret2 = 0 ;
			return ;
		}
		if(isRealStressField(f2))
			ret2 = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret1) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		else
			ret2 = (Vector) (parent->getBehaviour()->param * ret1) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		return ;
	}
	if(f1 == PRINCIPAL_STRAIN_FIELD && (f2 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f2 == PRINCIPAL_REAL_STRESS_FIELD))
	{
		Vector v1(0., 3+3*(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
		Vector v2(0., v1.size()) ;
		if(isRealStressField(f2))
			getField(STRAIN_FIELD, REAL_STRESS_FIELD, p, v1, v2, local, vm) ;
		else
			getField(STRAIN_FIELD, EFFECTIVE_STRESS_FIELD, p, v1, v2, local, vm) ;
		ret1 = toPrincipal(v1) ;
		ret2 = toPrincipal(v2) ;
		if (cleanup) delete vm ;
		return ;
	}
	if(isStrainField(f2) && isStressField(f1))
	{
		getField(f2, p, ret2, local, vm) ;
		if(parent->getBehaviour()->getTensor(p_, parent).numCols() != ret1.size())
		{
			ret1 = 0 ;
			if (cleanup) delete vm ;
			return ;
		}
		if(isRealStressField(f1))
			ret1 = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret2) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		else
			ret1 = (Vector) (parent->getBehaviour()->param * ret2) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		if (cleanup) delete vm ;
		return ;
	}
	if(f2 == PRINCIPAL_STRAIN_FIELD && (f1 == PRINCIPAL_EFFECTIVE_STRESS_FIELD || f1 == PRINCIPAL_REAL_STRESS_FIELD))
	{
		Vector v1(0., 3+3*(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)) ;
		Vector v2(0., v1.size()) ;
		if(isRealStressField(f2))
			getField(REAL_STRESS_FIELD, STRAIN_FIELD, p, v1, v2, local,vm) ;
		else
			getField(EFFECTIVE_STRESS_FIELD, STRAIN_FIELD, p, v1, v2, local,vm) ;
		ret1 = toPrincipal(v1) ;
		ret2 = toPrincipal(v2) ;
		if (cleanup) delete vm ;
		return ;
	}
	if(f1 == GRADIENT_FIELD && f2 == FLUX_FIELD)
	{
		getField(f1, p, ret1, local, vm) ;
		ret2 = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret1) ;
	}
	if(f1 == FLUX_FIELD && f2 == GRADIENT_FIELD)
	{
		getField(f2, p, ret2, local, vm) ;
		ret1 = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret2) ;
	}
  if (cleanup) delete vm ;
}

void ElementState::getField( FieldType f1, FieldType f2, const PointArray & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int)  const 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	Vector b1(0., ret1.size()/p.size()) ;
	Vector b2(0., ret2.size()/p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		getField(f1, f2, *p[i], b1, b2, local, vm) ;
		for(size_t j = b1.size()*i ; j < b1.size()*(i+1) ; j++)
			ret1[j] = b1[j - b1.size()*i] ;
		for(size_t j = b2.size()*i ; j < b2.size()*(i+1) ; j++)
			ret2[j] = b2[j - b2.size()*i] ;
	}
	if (cleanup) delete vm ;
}

void ElementState::getField( FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int , int)  const 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	
	Vector b1(0., ret1.size()/p.size()) ;
	Vector b2(0., ret2.size()/p.size()) ;
	for(size_t i = 0 ; i < p.size() ; i++)
	{
		getField(f1, f2, p[i].first, b1, b2, local, vm) ;
		for(size_t j = b1.size()*i ; j < b1.size()*(i+1) ; j++)
			ret1[j] = b1[j - b1.size()*i] ;
		for(size_t j = b2.size()*i ; j < b2.size()*(i+1) ; j++)
			ret2[j] = b2[j - b2.size()*i] ;
	}  
	if (cleanup) delete vm ;
}

void ElementState::getFieldAtGaussPoint( FieldType f1, FieldType f2, size_t p, Vector & ret1, Vector & ret2, VirtualMachine * vm, int i, int j) 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	Point p_ = parent->getGaussPoints().gaussPoints[p].first ;
	getField(f1, f2, p_, ret1, ret2, true, vm, i, j) ;
	if (cleanup) delete vm ;
}

FunctionMatrix ElementState::getStressFunction( const Matrix &Jinv , StressCalculationMethod m) const
{
	if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		FunctionMatrix ret( 2, 2 ) ;
		FunctionMatrix temp( 3, 1 ) ;

		Function x_xi ;
		Function x_eta ;
		Function y_xi ;
		Function y_eta ;

		for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
		{
			x_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 2] ;
			x_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 2] ;
			y_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 2 + 1] ;
			y_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 2 + 1] ;
		}

		temp[0][0] = x_xi * Jinv[0][0] + x_eta * Jinv[0][1] ;
		temp[2][0] = ( x_xi * Jinv[1][0] + x_eta * Jinv[1][1] + y_xi * Jinv[0][0] + y_eta * Jinv[0][1] ) * 0.5 ;
		temp[1][0] = y_xi * Jinv[1][0] + y_eta * Jinv[1][1] ;

		Matrix cg( parent->getBehaviour()->getTensor( getParent()->getCenter(), parent ) ) ;

		if (m == EFFECTIVE_STRESS)
			temp *=  parent->getBehaviour()->param;
		else
			temp *= cg;

		ret[0][0] = temp[0][0];
		ret[0][1] = temp[2][0];
		ret[1][0] = temp[2][0];
		ret[1][1] = temp[1][0];

		return ret ;
	}
	else     //if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		FunctionMatrix ret( 3, 3 ) ;
		FunctionMatrix temp( 6, 1 ) ;

		Function x_xi ;
		Function x_eta ;
		Function x_zeta ;
		Function y_xi ;
		Function y_eta ;
		Function y_zeta ;
		Function z_xi ;
		Function z_eta ;
		Function z_zeta ;

		for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
		{
			x_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 3] ;
			x_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 3] ;
			x_zeta += parent->getShapeFunction( j ).d( ZETA ) * displacements[j * 3] ;
			y_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 3 + 1] ;
			y_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 3 + 1] ;
			y_zeta += parent->getShapeFunction( j ).d( ZETA ) * displacements[j * 3 + 1] ;
			z_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 3 + 2] ;
			z_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 3 + 2] ;
			z_zeta += parent->getShapeFunction( j ).d( ZETA ) * displacements[j * 3 + 2] ;
		}

		temp[0][0] = x_xi * Jinv[0][0] + x_eta * Jinv[0][1] + x_zeta * Jinv[0][2];
		temp[1][0] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
		temp[2][0] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];
		temp[3][0] = ( ( y_xi ) * Jinv[2][0] +
		               ( y_eta ) * Jinv[2][1] +
		               ( y_zeta ) * Jinv[2][2] +
		               ( z_xi ) * Jinv[1][0] +
		               ( z_eta ) * Jinv[1][1] +
		               ( z_zeta ) * Jinv[1][2] ) * .5;
		temp[4][0] = ( ( x_xi ) * Jinv[2][0] +
		               ( x_eta ) * Jinv[2][1] +
		               ( x_zeta ) * Jinv[2][2] +
		               ( z_xi ) * Jinv[0][0] +
		               ( z_eta ) * Jinv[0][1] +
		               ( z_zeta ) * Jinv[0][2] ) * .5;
		temp[5][0] = ( ( y_xi )  * Jinv[0][0] +
		               ( y_eta )  * Jinv[0][1] +
		               ( y_zeta ) * Jinv[0][2] +
		               ( x_xi )   * Jinv[1][0] +
		               ( x_eta )  * Jinv[1][1] +
		               ( x_zeta ) * Jinv[1][2] ) * .5;

		Matrix cg( parent->getBehaviour()->getTensor( getParent()->getCenter(), parent ) ) ;

		if (m == EFFECTIVE_STRESS)
			temp *=parent->getBehaviour()->param ;
		else
			temp *=cg ;

		ret[0][0] = temp[0][0] ;
		ret[0][1] = temp[5][0] ;
		ret[0][2] = temp[4][0] ;
		ret[1][0] = temp[5][0] ;
		ret[1][1] = temp[1][0] ;
		ret[1][2] = temp[3][0] ;
		ret[2][0] = temp[4][0] ;
		ret[2][1] = temp[3][0] ;
		ret[2][2] = temp[2][0] ;

		return ret ;
	}

}

FunctionMatrix ElementState::getStrainFunction( const Matrix &Jinv) const
{
	if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		FunctionMatrix ret( 2, 2 ) ;

		Function x_xi ;
		Function x_eta ;
		Function y_xi ;
		Function y_eta ;

		for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
		{
			x_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 2] ;
			x_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 2] ;
			y_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 2 + 1] ;
			y_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 2 + 1] ;
		}

		ret[0][0] = x_xi * Jinv[0][0] + x_eta * Jinv[0][1] ;
		ret[0][1] = ( x_xi * Jinv[1][0] + x_eta * Jinv[1][1] + y_xi * Jinv[0][0] + y_eta * Jinv[0][1] ) * 0.5;
		ret[1][0] = ( x_xi * Jinv[1][0] + x_eta * Jinv[1][1] + y_xi * Jinv[0][0] + y_eta * Jinv[0][1] ) * 0.5;
		ret[1][1] = y_xi * Jinv[1][0] + y_eta * Jinv[1][1] ;

		return ret ;
	}
	else     //if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		FunctionMatrix ret( 3, 3 ) ;
		FunctionMatrix temp( 6, 1 ) ;

		Function x_xi ;
		Function x_eta ;
		Function x_zeta ;
		Function y_xi ;
		Function y_eta ;
		Function y_zeta ;
		Function z_xi ;
		Function z_eta ;
		Function z_zeta ;

		for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
		{
			x_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 3] ;
			x_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 3] ;
			x_zeta += parent->getShapeFunction( j ).d( ZETA ) * displacements[j * 3] ;
			y_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 3 + 1] ;
			y_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 3 + 1] ;
			y_zeta += parent->getShapeFunction( j ).d( ZETA ) * displacements[j * 3 + 1] ;
			z_xi += parent->getShapeFunction( j ).d( XI ) * displacements[j * 3 + 2] ;
			z_eta += parent->getShapeFunction( j ).d( ETA ) * displacements[j * 3 + 2] ;
			z_zeta += parent->getShapeFunction( j ).d( ZETA ) * displacements[j * 3 + 2] ;
		}

		temp[0][0] = x_xi * Jinv[0][0] + x_eta * Jinv[0][1] + x_zeta * Jinv[0][2];
		temp[1][0] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1]  + ( y_zeta ) * Jinv[1][2];
		temp[2][0] = ( z_xi ) * Jinv[2][0] + ( z_eta ) * Jinv[2][1]  + ( z_zeta ) * Jinv[2][2];
		temp[3][0] = ( ( y_xi ) * Jinv[2][0] +
		               ( y_eta ) * Jinv[2][1] +
		               ( y_zeta ) * Jinv[2][2] +
		               ( z_xi ) * Jinv[1][0] +
		               ( z_eta ) * Jinv[1][1] +
		               ( z_zeta ) * Jinv[1][2] ) * .5;
		temp[4][0] = ( ( x_xi ) * Jinv[2][0] +
		               ( x_eta ) * Jinv[2][1] +
		               ( x_zeta ) * Jinv[2][2] +
		               ( z_xi ) * Jinv[0][0] +
		               ( z_eta ) * Jinv[0][1] +
		               ( z_zeta ) * Jinv[0][2] ) * .5;
		temp[5][0] = ( ( y_xi )  * Jinv[0][0] +
		               ( y_eta )  * Jinv[0][1] +
		               ( y_zeta ) * Jinv[0][2] +
		               ( x_xi )   * Jinv[1][0] +
		               ( x_eta )  * Jinv[1][1] +
		               ( x_zeta ) * Jinv[1][2] ) * .5;


		ret[0][0] = temp[0][0] ;
		ret[0][1] = temp[5][0] ;
		ret[0][2] = temp[4][0] ;
		ret[1][0] = temp[5][0] ;
		ret[1][1] = temp[1][0] ;
		ret[1][2] = temp[3][0] ;
		ret[2][0] = temp[4][0] ;
		ret[2][1] = temp[3][0] ;
		ret[2][2] = temp[2][0] ;

		return ret ;
	}

}

FunctionMatrix ElementState::getDisplacementFunction() const
{
	if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL )
	{
		FunctionMatrix ret( 2, 1 ) ;

		for( size_t j = 0 ; j < parent->getBoundingPoints().size() &&  j < displacements.size() * 2; j++ )
		{
			ret[0][0] += parent->getShapeFunction( j ) * displacements[j * 2] ;
			ret[1][0] += parent->getShapeFunction( j ) * displacements[j * 2 + 1] ;
		}

		for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
		{
			ret[0][0] += parent->getEnrichmentFunction( j ) * enrichedDisplacements[j * 2] ;
			ret[1][0] += parent->getEnrichmentFunction( j ) * enrichedDisplacements[j * 2 + 1] ;
		}

		return ret ;
	}
	else     //if (parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL)
	{
		FunctionMatrix ret( 3, 1 ) ;

		for( size_t j = 0 ; j < parent->getBoundingPoints().size() &&  j < displacements.size() * 3; j++ )
		{
			ret[0][0] += parent->getShapeFunction( j ) * displacements[j * 3] ;
			ret[1][0] += parent->getShapeFunction( j ) * displacements[j * 3 + 1] ;
			ret[2][0] += parent->getShapeFunction( j ) * displacements[j * 3 + 2] ;
		}

		for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 3; j++ )
		{
			ret[0][0] += parent->getEnrichmentFunction( j ) * enrichedDisplacements[j * 3] ;
			ret[1][0] += parent->getEnrichmentFunction( j ) * enrichedDisplacements[j * 3 + 1] ;
			ret[2][0] += parent->getEnrichmentFunction( j ) * enrichedDisplacements[j * 3 + 2] ;
		}

		return ret ;
	}
}

std::vector<double> ElementState::getEnrichedInterpolatingFactors( const Point &p, bool local ) const
{

	std::vector<double> ret( parent->getEnrichmentFunctions().size() ) ;
	VirtualMachine vm ;
	Point p_ = p ;

	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;

	for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
	{
		ret[j] = vm.eval( parent->getEnrichmentFunction( j ), p_ ) ;
	}

	return ret;
}

std::vector<double> ElementState::getNonEnrichedInterpolatingFactors( const Point &p, bool local ) const
{

	std::vector<double> ret( parent->getBoundingPoints().size() ) ;
	VirtualMachine vm ;
	Point p_( p ) ;

	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;

	for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
	{
		ret[j] = vm.eval( parent->getShapeFunction( j ), p_ ) ;
	}

	return ret;
}

std::vector<double> ElementState::getInterpolatingFactors( const Point &p, bool local ) const
{

	std::vector<double> ret( parent->getBoundingPoints().size() + parent->getEnrichmentFunctions().size() ) ;
	VirtualMachine vm ;
	Point p_ = p ;

	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;

	for( size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++ )
	{
		ret[j] = vm.eval( parent->getShapeFunction( j ), p_ ) ;
	}

	for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++ )
	{
		ret[j + parent->getBoundingPoints().size()] = vm.eval( parent->getEnrichmentFunction( j ), p_ ) ;
	}

	return ret;
}

void ElementState::initialize( bool initializeFractureCache)
{
	size_t ndofs = 0 ;
	if(parent->getBehaviour()) 
		ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	displacements.resize( parent->getBoundingPoints().size()*ndofs ) ;
	displacements = 0 ;
	previousDisplacements.resize( displacements.size() ) ;
	previousDisplacements = 0 ;

	buffer.resize( displacements.size() ) ;
	buffer = 0 ;

	if( std::abs( timePos - previousTimePos ) < POINT_TOLERANCE_3D && std::abs( timePos ) < POINT_TOLERANCE_3D )
	{
		timePos = -0.1 ;
		previousTimePos = -0.2 ;
	}

	if( initializeFractureCache && parent->getBehaviour()->getFractureCriterion() )
	{
		parent->getBehaviour()->getFractureCriterion()->initialiseCache( *this ) ;
	}

}

const Vector &ElementState::getBuffer() const
{
	return buffer ;
}

Vector &ElementState::getBuffer()
{
	return buffer ;
}

void ElementState::step( double dt, const Vector *d )
{
	
	strainAtGaussPoints.resize( 0 );
	stressAtGaussPoints.resize( 0 );	
	
	pstrainAtGaussPoints.resize( 0 );
	pstressAtGaussPoints.resize( 0 );
	
	effectivePStressAtGaussPoints.resize(0);

	if( parent->getBehaviour()&& parent->getBehaviour()->type != VOID_BEHAVIOUR )
	{
		previousTimePos = timePos ;
		timePos += dt ;
		size_t ndofs = parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;

		std::vector< size_t > ids = parent->getDofIds() ;

		if( ids.empty() )
			return ;

		if( buffer.size() != parent->getBoundingPoints().size()*ndofs + parent->getEnrichmentFunctions().size()*ndofs )
			buffer.resize( parent->getBoundingPoints().size()*ndofs + parent->getEnrichmentFunctions().size()*ndofs ) ;

		for( size_t i = 0 ; i < parent->getShapeFunctions().size() ; i++ )
		{
			for( size_t j = 0 ; j < ndofs ; j++ )
			{
				if(ids[i] * ndofs + j < d->size())
					buffer[i * ndofs + j] = ( *d )[ids[i] * ndofs + j] ;
				else
					buffer[i * ndofs + j] = 0 ;
			}
		}

		int nbp = parent->getBoundingPoints().size() ;

		for( size_t i = 0 ; i < parent->getEnrichmentFunctions().size() ; i++ )
		{
			for( size_t j = 0 ; j < ndofs ; j++ )
			{
				if(ids[i + nbp] * ndofs + j < d->size())
					buffer[i * ndofs + nbp * ndofs + j] = ( *d )[ids[i + nbp] * ndofs + j] ;
				else
					buffer[i * ndofs + nbp * ndofs + j] = 0 ;
			}
		}
	}
}


double ElementState::getTime() const
{
	return timePos ;
}

double ElementState::getDeltaTime() const
{
	return timePos - previousTimePos ;
}

double ElementState::getNodalDeltaTime() const
{
	return parent->getBoundingPoint(parent->getBoundingPoints().size() -1).t - parent->getBoundingPoint(0).t ;
}

FieldType stressFieldType(StressCalculationMethod m)
{
	return m == REAL_STRESS ? REAL_STRESS_FIELD : EFFECTIVE_STRESS_FIELD ;
}

FieldType principalStressFieldType(StressCalculationMethod m)
{
	return m == REAL_STRESS ? PRINCIPAL_REAL_STRESS_FIELD : PRINCIPAL_EFFECTIVE_STRESS_FIELD ;
}

std::vector<Point> ElementState::getPrincipalFrame( const Point &p, bool local, VirtualMachine * vm, StressCalculationMethod m ) 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	
	if(getParent()->spaceDimensions() == SPACE_TWO_DIMENSIONAL)
	{
		Vector principalStresses(0., 2) ;
		getAverageField(principalStressFieldType(m), principalStresses, vm) ;
		double principalAngle = 0.5*atan2(principalStresses[0]-principalStresses[1], 2.*principalStresses[2]) ;
		std::vector<Point> ret ;
		ret.push_back(Point(cos(principalAngle), -sin(principalAngle)));
		ret.push_back(Point(sin(principalAngle), cos(principalAngle)));
		ret.push_back(Point(0,0, 1));
		if(cleanup) delete vm ;
		return ret ;
	}
	
	Vector stress(0., 6) ;
	Point p_ = p ;
	if(!local)
		p_ = parent->inLocalCoordinates(p_) ;
	this->getField(principalStressFieldType(m), p_, stress, true) ;
	Matrix stressMatrix = makeStressOrStrainMatrix(stress) ;
	
	//then, we get the principal stresses
	Vector principalStresses(0., 3) ;
	this->getField(principalStressFieldType(m), p_, principalStresses, true) ;
	std::vector<Vector> principalVectors ;
	for(int i = 0 ; i <  principalStresses.size() ; ++i)
	{
		principalVectors.push_back(Vector(0.,stressMatrix.numCols())) ;
		Matrix m = (stressMatrix-identity(stressMatrix.numCols())*principalStresses[i]) ;
		//Assume that the first coefficient is 1 and get a reduced system.
		Matrix m_reduced (m) ;
		Vector v_reduced(0., stressMatrix.numCols()) ;
		v_reduced[0] = 1 ;
		m_reduced[0][0] = 1 ;
		for(size_t j = 1 ; j < stressMatrix.numCols() ; j++)
		{
			m_reduced[0][j] = 0 ;
		}
		for(size_t j = 1 ; j < stressMatrix.numRows() ; j++)
		{
			v_reduced[j] -= m_reduced[j][0] ;
			m_reduced[j][0] = 0 ;
		}
		
		int k = 1 ; 
		while(det(m_reduced) < 1e-6 && k < m_reduced.numCols())
		{
			m_reduced = m ;
			v_reduced.resize( stressMatrix.numCols(),0.) ;
			v_reduced[k] = 1 ;
			m_reduced[k][k] = 1 ;
			for(size_t j = 0 ; j < stressMatrix.numCols() ; j++)
			{
				if(j == k)
					continue ;
				m_reduced[k][j] = 0 ;
			}
			for(size_t j = 0 ; j < stressMatrix.numRows() ; j++)
			{
				if(j == k)
					continue ;
				v_reduced[j] -= m_reduced[j][k] ;
				m_reduced[j][k] = 0 ;
			}
		}
		solveSystem( m_reduced, v_reduced, principalVectors.back());
	}
	std::vector<Point> ret ;
	

	if(stressMatrix.numCols() == 2)
	{
		for(int i = 0 ; i <  principalVectors.size() ; ++i)
		{
			ret.push_back(Point(principalVectors[i][0], principalVectors[i][1]));
			ret.back() /= ret.back().norm() ;
		}
		ret.push_back(Point(0, 0, 1));
	}
	else
	{
		for(int i = 0 ; i <  principalVectors.size() ; ++i)
		{
			ret.push_back(Point(principalVectors[i][0], principalVectors[i][1], principalVectors[i][2]));
			ret.back() /= ret.back().norm() ;
		}
	}
	if(cleanup) delete vm ;
	return ret ;
}



double ElementState::elasticEnergy()
{
	Vector strain( 0., parent->getGaussPoints().gaussPoints.size()*(3+3*(parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL))) ;
	Vector stress( 0., strain.size()) ;
	this->getField( STRAIN_FIELD, REAL_STRESS_FIELD, parent->getGaussPoints().gaussPoints , strain, stress, true ) ;
	size_t stresscomponents = stress.size() / parent->getGaussPoints().gaussPoints.size() ;
	double e = 0 ;

	for( size_t i = 0 ; i < parent->getGaussPoints().gaussPoints.size() ; i++ )
	{
		double le = 0 ;

		for( size_t j = 0 ; j < stresscomponents ; j++ )
		{
			le += strain[i * stresscomponents + j] * stress[i * stresscomponents + j] ;
		}

		e += le * parent->getGaussPoints().gaussPoints[i].second ;
	}

	return .5 * e ;

}

ElementStateWithInternalVariables::ElementStateWithInternalVariables(IntegrableEntity * e, int n_, int p_) : ElementState(e), p(p_), n(n_)
{
	
}

ElementStateWithInternalVariables::ElementStateWithInternalVariables(const ElementStateWithInternalVariables &s)  : ElementState(dynamic_cast<const ElementState &>(s))
{
	n = s.numberOfInternalVariables() ;
	p = s.sizeOfInternalVariable() ;
  
	internalVariablesAtGaussPoints.resize( parent->getGaussPoints().gaussPoints.size() ) ;
	
	for(size_t g = 0 ; g < internalVariablesAtGaussPoints.size() ; g++)
	{
		internalVariablesAtGaussPoints[g].resize( n ) ;
		for(size_t k = 0 ; k < n ; k++)
			internalVariablesAtGaussPoints[g][k].resize(p) ;
	}
  
}
						
ElementStateWithInternalVariables & ElementStateWithInternalVariables::operator =(const ElementStateWithInternalVariables & s) 
{
	strainAtGaussPoints.resize(0);
	stressAtGaussPoints.resize(0);
	effectivePStressAtGaussPoints.resize(0);
	displacements.resize( s.getDisplacements().size() ) ;
	displacements = s.getDisplacements() ;
	enrichedDisplacements.resize( s.getEnrichedDisplacements().size() ) ;
	enrichedDisplacements  = s.getEnrichedDisplacements();

	previousDisplacements.resize( s.getPreviousDisplacements().size() ) ;
	previousDisplacements = s.getPreviousDisplacements() ;
	previousEnrichedDisplacements.resize( s.getPreviousEnrichedDisplacements().size() ) ;
	previousEnrichedDisplacements = s. getPreviousEnrichedDisplacements();

	buffer.resize( s.getBuffer().size() ) ;
	buffer = s.getBuffer() ;

	timePos = s.getTime();
	previousTimePos = s.getTime() - s.getDeltaTime();

	parent = s.getParent();

	n = s.numberOfInternalVariables() ;
	p = s.sizeOfInternalVariable() ;
	internalVariablesAtGaussPoints.resize( parent->getGaussPoints().gaussPoints.size() ) ;
	
	for(size_t g = 0 ; g < internalVariablesAtGaussPoints.size() ; g++)
	{
		internalVariablesAtGaussPoints[g].resize( n ) ;
		for(size_t k = 0 ; k < n ; k++)
			internalVariablesAtGaussPoints[g][k].resize(p) ;
	}
	
	return *this ;  
}

void ElementStateWithInternalVariables::getField(FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine * vm, int i) const
{
	if(f != INTERNAL_VARIABLE_FIELD)
		ElementState::getField(f, p, ret, local, vm, i) ;
}

void ElementStateWithInternalVariables::getField(FieldType f, const PointArray & p, Vector & ret, bool local, VirtualMachine * vm, int i) const
{
	if(f != INTERNAL_VARIABLE_FIELD)
		ElementState::getField(f, p, ret, local, vm, i) ;
}

void ElementStateWithInternalVariables::getField(FieldType f, const std::valarray<std::pair<Point, double> > & p, Vector & ret, bool local, VirtualMachine * vm, int i) const
{
	if(f != INTERNAL_VARIABLE_FIELD)
		ElementState::getField(f, p, ret, local, vm, i) ;
}

// void ElementStateWithInternalVariables::getFieldAtNodes( FieldType f, Vector & ret, int i) 
// {
// 	if(f != INTERNAL_VARIABLE_FIELD)
// 		ElementState::getFieldAtNodes(f, ret, i) ; 
// }

void ElementStateWithInternalVariables::getFieldAtGaussPoint( FieldType f, size_t g, Vector & ret, VirtualMachine * vm, int i) 
{
	if(f == INTERNAL_VARIABLE_FIELD)
	{
		ret = internalVariablesAtGaussPoints[g][i] ;
		return ;
	}
	else
	{
		ElementState::getFieldAtGaussPoint(f, g, ret, vm, i) ;
// 		Point p_ = parent->getGaussPoints().gaussPoints[g].first ;
// 		this->getField(f, p_, ret, true, i) ;
	}
}

void ElementStateWithInternalVariables::getField(FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine *vm, int i, int j) const
{
	if(f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD)
		ElementState::getField(f1, f2, p, ret1, ret2, local, vm, i, j) ;
}

void ElementStateWithInternalVariables::getField(FieldType f1, FieldType f2, const PointArray & p,  Vector & ret1, Vector & ret2, bool local, VirtualMachine *vm, int i, int j) const
{
	if(f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD)
		ElementState::getField(f1, f2, p, ret1, ret2, local, vm, i, j) ;
}

void ElementStateWithInternalVariables::getField(FieldType f1, FieldType f2, const std::valarray<std::pair<Point, double> > & p,  Vector & ret1, Vector & ret2, bool local, VirtualMachine *vm, int i, int j) const
{
	if(f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD)
		ElementState::getField(f1, f2, p, ret1, ret2, local, vm, i, j) ;
}

// void ElementStateWithInternalVariables::getFieldAtNodes(FieldType f1, FieldType f2,  Vector & ret1, Vector & ret2, int i, int j) 
// {
// 	if(f1 != INTERNAL_VARIABLE_FIELD && f2 != INTERNAL_VARIABLE_FIELD)
// 		ElementState::getFieldAtNodes(f1, f2, ret1, ret2, i, j) ;
// }

void ElementStateWithInternalVariables::getFieldAtGaussPoint(FieldType f1, FieldType f2, size_t g,  Vector & ret1, Vector & ret2, VirtualMachine * vm, int i, int j) 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	bool done1 = false ;
	bool done2 = false ;
	if(f1 == INTERNAL_VARIABLE_FIELD)
	{
		ret1 = internalVariablesAtGaussPoints[g][i] ;
		done1 = true ;
	}
	if(f2 == INTERNAL_VARIABLE_FIELD)
	{
		ret2 = internalVariablesAtGaussPoints[g][j] ;
		done2 = true ;
	}
	if(done1 && done2)
		return ;
	Point p_ = parent->getGaussPoints().gaussPoints[g].first ;
	if(done1)
	{
		ElementState::getField(f2, p_, ret2, true, vm, j) ;
		if(cleanup) delete vm ;
		return ;
	}
	if(done2)
	{
		ElementState::getField(f1, p_, ret1, true, vm, i) ;
		if(cleanup) delete vm ;
		return ;
	}
	ElementState::getField(f1, f2, p_, ret1, ret2, true, vm, i, j) ;
	if(cleanup) delete vm ;
}

void ElementStateWithInternalVariables::initialize( bool initializeFractureCache)
{
	ElementState::initialize(initializeFractureCache) ;
	
	size_t ngp = parent->getGaussPoints().gaussPoints.size() ;
	std::vector<Vector> dummy(0) ;
	Vector vdummy(0., p) ;

	for(size_t g = 0 ; g < ngp ; g++)
	{
		internalVariablesAtGaussPoints.push_back(dummy) ;
		for(size_t k = 0 ; k < n ; k++)
		{
			internalVariablesAtGaussPoints[g].push_back(vdummy) ;
		}
	}
	
}

void ElementStateWithInternalVariables::setNumberOfGaussPoints(int g) 
{
	while(internalVariablesAtGaussPoints.size() < g)
	{
		std::vector<Vector> dummy(0) ;
		for(size_t i = 0 ; i < n ; i++)
		{
			dummy.push_back(internalVariablesAtGaussPoints[0][i]) ;
		}
		internalVariablesAtGaussPoints.push_back(dummy) ;
	}
	while(internalVariablesAtGaussPoints.size() > g)
	{
		internalVariablesAtGaussPoints.pop_back() ;
	}
}


void ElementStateWithInternalVariables::setInternalVariableAtGaussPoint(Vector & v, size_t g, int i) 
{
	internalVariablesAtGaussPoints[g][i] = v ;
}

int Mu::isGaussPoint(const Point & p, IntegrableEntity * e) 
{
	if(e)
	{
		GaussPointArray gauss = e->getGaussPoints() ; 
		for(size_t i = 0 ; i < gauss.gaussPoints.size() ; i++)
		{
			if(p == gauss.gaussPoints[i].first)
				return i ;
		}
	}
	return -1 ;
}

ParallelElementState::ParallelElementState(IntegrableEntity * e, std::vector<ElementState *> s) : ElementState(e)
{
	for(size_t i = 0 ; i < s.size() ; i++)
		states.push_back(s[i]) ;
}

ParallelElementState::ParallelElementState(const ParallelElementState &s) : ElementState(s)
{
	for(size_t i = 0 ; i < s.getNumberOfStates() ; i++)
		states.push_back( new ElementState( s.getState(i) ) ) ;
}


ParallelElementState & ParallelElementState::operator =(const ParallelElementState & s) 
{
	ElementState::operator =(s) ;
	for(size_t i = 0 ; i < s.getNumberOfStates() ; i++)
		this->getState(i) = s.getState(i) ;
	return *this ;
}

ElementState & ParallelElementState::getState(size_t i) 
{
	return *states[i] ;
}

const ElementState & ParallelElementState::getState(size_t i) const
{
	return *states[i] ;
}

void ParallelElementState::initialize(bool initializeFractureCache) 
{
	ElementState::initialize( initializeFractureCache) ;
	for(size_t i = 0 ; i < states.size() ; i++)
		states[i]->initialize( initializeFractureCache) ;
}

void ParallelElementState::step(double dt, const Vector* d ) 
{
	ElementState::step(dt, d) ;
	for(size_t i = 0 ; i < states.size() ; i++)
		states[i]->step(dt,d) ;
}

SerialElementState::SerialElementState(IntegrableEntity * e, std::vector<ElementState *> s) : ElementState(e)
{
	for(size_t i = 0 ; i < s.size() ; i++)
		states.push_back(s[i]) ;
}

SerialElementState::SerialElementState(const SerialElementState &s) : ElementState(s)
{
	for(size_t i = 0 ; i < s.getNumberOfStates() ; i++)
		states.push_back( new ElementState( s.getState(i) ) ) ;
}


SerialElementState & SerialElementState::operator =(const SerialElementState & s) 
{
	ElementState::operator =(s) ;
	for(size_t i = 0 ; i < s.getNumberOfStates() ; i++)
		this->getState(i) = s.getState(i) ;
	
	return *this ;
}

ElementState & SerialElementState::getState(size_t i) 
{
	return *states[i] ;
}

const ElementState & SerialElementState::getState(size_t i) const
{
	return *states[i] ;
}

void SerialElementState::initialize(bool initializeFractureCache) 
{
	ElementState::initialize( initializeFractureCache) ;
	for(size_t i = 0 ; i < states.size() ; i++)
		states[i]->initialize( initializeFractureCache) ;
}

void SerialElementState::step(double dt, const Vector* d ) 
{
	ElementState::step(dt, d) ;
	for(size_t i = 0 ; i < states.size() ; i++)
		states[i]->step(dt,d) ;
}

KelvinVoightSpaceTimeElementState::KelvinVoightSpaceTimeElementState(IntegrableEntity * e) : ElementState(e)
{

  
}

KelvinVoightSpaceTimeElementState::KelvinVoightSpaceTimeElementState(const KelvinVoightSpaceTimeElementState &s) : ElementState(s)
{

}


KelvinVoightSpaceTimeElementState & KelvinVoightSpaceTimeElementState::operator =(const KelvinVoightSpaceTimeElementState & s) 
{
	ElementState::operator =(s) ;
	
	return *this ;
}

void KelvinVoightSpaceTimeElementState::getField( FieldType f, const Point & p, Vector & ret, bool local, VirtualMachine *vm, int )  const 
{  
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	int n = 0 ;
	Point p_ = p ;
	if( !local )
		p_ = parent->inLocalCoordinates( p ) ;

	switch(f)
	{
		case SPEED_FIELD:
			n =  parent->getBehaviour()->getNumberOfDegreesOfFreedom() ;
			for(size_t j = 0 ; j < parent->getBoundingPoints().size() ; j++)
			{
				double f =  vm ->deval( parent->getShapeFunction( j ) , TIME_VARIABLE, p_) ;
				for(size_t k = 0 ; k < n ; k++)
					ret[k] += f * displacements[j*n+k] ;
			}
			for(size_t j = 0 ; j < parent->getEnrichmentFunctions().size() ; j++)
			{
				double f =  vm ->deval( parent->getShapeFunction( j ) , TIME_VARIABLE, p_) ;
				for(size_t k = 0 ; k < n ; k++)
					ret[k] += f * enrichedDisplacements[j*n+k] ;
			}
			if(cleanup) delete vm ;
			return ;
		case STRAIN_RATE_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2 )
			{
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm ->ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ , 1e-5) ;
					double f_eta = vm ->ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ , 1e-5) ;
					x_xi += f_xi * displacements[j * 2] ;
					x_eta += f_eta * displacements[j * 2] ;
					y_xi += f_xi * displacements[j * 2 + 1] ;
					y_eta += f_eta * displacements[j * 2 + 1] ;
				}

				for( size_t j = 0 ; j < parent->getEnrichmentFunctions().size() && j < enrichedDisplacements.size() * 2; j++ )
				{
					double f_xi = vm ->ddeval( parent->getEnrichmentFunction( j ), XI, TIME_VARIABLE, p_ ) ;
					double f_eta = vm ->ddeval( parent->getEnrichmentFunction( j ), ETA, TIME_VARIABLE, p_ ) ;

					x_xi += f_xi * enrichedDisplacements[j * 2] ;
					x_eta += f_eta * enrichedDisplacements[j * 2] ;
					y_xi += f_xi * enrichedDisplacements[j * 2 + 1] ;
					y_eta += f_eta * enrichedDisplacements[j * 2 + 1] ;
				}

				Matrix Jinv ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
				ret[2] =  ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
				ret *= Jinv[2][2] ;
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3 )
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
					double f_xi = vm ->ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ ) ;
					double f_eta = vm ->ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ ) ;
					double f_zeta = vm ->ddeval( parent->getShapeFunction( j ), ZETA, TIME_VARIABLE, p_ ) ;
					double x = displacements[j * 3] ;
					double y = displacements[j * 3 + 1] ;
					double z = displacements[j * 3 + 2] ;

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
					double f_xi = vm ->ddeval( parent->getEnrichmentFunction( j ), XI, TIME_VARIABLE, p_ ) ;
					double f_eta = vm ->ddeval( parent->getEnrichmentFunction( j ), ETA, TIME_VARIABLE, p_ ) ;
					double f_zeta = vm ->ddeval( parent->getEnrichmentFunction( j ), ZETA, TIME_VARIABLE, p_ ) ;
					double x = enrichedDisplacements[j * 3] ;
					double y = enrichedDisplacements[j * 3 + 1] ;
					double z = enrichedDisplacements[j * 3 + 2] ;

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

				ret[3] =  ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] =  ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] =  ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
			}
			if(cleanup) delete vm ;
			return ;
		case NON_ENRICHED_STRAIN_RATE_FIELD:
			if( parent->spaceDimensions() == SPACE_TWO_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 2 )
			{
				double x_xi = 0;
				double x_eta = 0;
				double y_xi = 0;
				double y_eta = 0;

				for( size_t j = 0 ; j < parent->getBoundingPoints().size(); j++ )
				{
					double f_xi = vm ->ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ ,1e-4) ;
					double f_eta = vm ->ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ ,1e-4) ;
					x_xi += f_xi * displacements[j * 2] ;
					x_eta += f_eta * displacements[j * 2] ;
					y_xi += f_xi * displacements[j * 2 + 1] ;
					y_eta += f_eta * displacements[j * 2 + 1] ;
				}

				Matrix Jinv( 3, 3 ) ;
				parent->getInverseJacobianMatrix( p_, Jinv ) ;
				ret[0] = ( x_xi ) * Jinv[0][0] + ( x_eta ) * Jinv[0][1] ;
				ret[1] = ( y_xi ) * Jinv[1][0] + ( y_eta ) * Jinv[1][1] ;
				ret[2] =  ( ( x_xi ) * Jinv[1][0] + ( x_eta ) * Jinv[1][1]  + ( y_xi ) * Jinv[0][0] + ( y_eta ) * Jinv[0][1] );
			}
			else if( parent->spaceDimensions() == SPACE_THREE_DIMENSIONAL && parent->getBehaviour()->getNumberOfDegreesOfFreedom() == 3 )
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
					double f_xi = vm ->ddeval( parent->getShapeFunction( j ), XI, TIME_VARIABLE, p_ ) ;
					double f_eta = vm ->ddeval( parent->getShapeFunction( j ), ETA, TIME_VARIABLE, p_ ) ;
					double f_zeta = vm ->ddeval( parent->getShapeFunction( j ), ZETA, TIME_VARIABLE, p_ ) ;
					double x = displacements[j * 3] ;
					double y = displacements[j * 3 + 1] ;
					double z = displacements[j * 3 + 2] ;

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

				ret[3] =  ( ( y_xi ) * Jinv[2][0] +
						    ( y_eta ) * Jinv[2][1] +
						    ( y_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[1][0] +
						    ( z_eta ) * Jinv[1][1] +
						    ( z_zeta ) * Jinv[1][2] );

				ret[4] =  ( ( x_xi ) * Jinv[2][0] +
						    ( x_eta ) * Jinv[2][1] +
						    ( x_zeta ) * Jinv[2][2] +
						    ( z_xi ) * Jinv[0][0] +
						    ( z_eta ) * Jinv[0][1] +
						    ( z_zeta ) * Jinv[0][2] );

				ret[5] =  ( ( y_xi )   * Jinv[0][0] +
						    ( y_eta )  * Jinv[0][1] +
						    ( y_zeta ) * Jinv[0][2] +
						    ( x_xi )   * Jinv[1][0] +
						    ( x_eta )  * Jinv[1][1] +
						    ( x_zeta ) * Jinv[1][2] );
			}
			if(cleanup) delete vm ;
			return ;
		case REAL_STRESS_FIELD:
		{
			Vector strain(0., ret.size()) ;
			Vector rate(0., ret.size()) ;
			this->getField(STRAIN_FIELD, STRAIN_RATE_FIELD, p_, strain, rate, true, vm) ;
			if(parent->getBehaviour()->getTensor(p_, parent).numCols() != ret.size())
			{
				ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * strain) ;
			ret += (Vector) (dynamic_cast<KelvinVoight *>(parent->getBehaviour())->eta * rate ) ;
			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			if(cleanup) delete vm ;
			return ;
		}
		case NON_ENRICHED_REAL_STRESS_FIELD:
		{
			Vector strain(0., ret.size()) ;
			Vector rate(0., ret.size()) ;
			this->getField(NON_ENRICHED_STRAIN_FIELD, NON_ENRICHED_STRAIN_RATE_FIELD, p_, strain, rate, true, vm) ;
			if(parent->getBehaviour()->getTensor(p_, parent).numCols() != ret.size())
			{
				ret = 0 ;
				if(cleanup) delete vm ;
				return ;
			}
			ret = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * strain) ;
			ret += (Vector) (dynamic_cast<KelvinVoight *>(parent->getBehaviour())->eta * rate ) ;
			ret -= getParent()->getBehaviour()->getImposedStress(p_, parent) ;
			if(cleanup) delete vm ;
			return ;
		}
		case EFFECTIVE_STRESS_FIELD:
		{
			Vector strain(0., ret.size()) ;
			Vector rate(0., ret.size()) ;
			this->getField(STRAIN_FIELD, STRAIN_RATE_FIELD, p_, strain, rate, true, vm) ;
			
			ret = (Vector) (parent->getBehaviour()->param * strain) ;
			ret += (Vector) (dynamic_cast<KelvinVoight *>(parent->getBehaviour())->eta * rate ) ;
			ret -= getParent()->getBehaviour()->getImposedStrain(p_, parent)*parent->getBehaviour()->param ;
			if(cleanup) delete vm ;
			return ;
		}
		case NON_ENRICHED_EFFECTIVE_STRESS_FIELD:
		{
			Vector strain(0., ret.size()) ;
			Vector rate(0., ret.size()) ;
			this->getField(NON_ENRICHED_STRAIN_FIELD, NON_ENRICHED_STRAIN_RATE_FIELD, p_, strain, rate, true, vm) ;
			ret = (Vector) (parent->getBehaviour()->param * strain) ;
			ret += (Vector) (dynamic_cast<KelvinVoight *>(parent->getBehaviour())->eta * rate ) ;
			ret -= getParent()->getBehaviour()->getImposedStrain(p_, parent)*parent->getBehaviour()->param ;
			if(cleanup) delete vm ;
			return ;
		}
	}
	
	ElementState::getField( f, p, ret, local, vm) ;
	if(cleanup) delete vm ;
}

void KelvinVoightSpaceTimeElementState::getField( FieldType f1, FieldType f2, const Point & p, Vector & ret1, Vector & ret2, bool local, VirtualMachine * vm, int i, int j)  const 
{
	bool cleanup = !vm ;
	if(!vm) vm = new VirtualMachine() ;
	Point p_ = p ;
	if(!local)
		p_ = parent->inLocalCoordinates(p) ;
	if(isStrainField(f1) && isStressField(f2))
	{
		this->getField(f1, p, ret1, local) ;
		if(isRealStressField(f2))
			ret2 = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret1) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		else
			ret2 = (Vector) (parent->getBehaviour()->param * ret1) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		Vector rate(0., ret1.size()) ;
		if(f1 == STRAIN_FIELD)
			this->getField( STRAIN_RATE_FIELD, p, rate, local, vm) ;
		else
			this->getField( NON_ENRICHED_STRAIN_RATE_FIELD, p, rate, local, vm) ;
		ret2 += (Vector) (dynamic_cast<KelvinVoight *>(parent->getBehaviour())->eta * rate) ;
		if(cleanup) delete vm ;
		return ;
	}
	if(isStrainField(f2) && isStressField(f1))
	{
		this->getField(f2, p, ret2, local, vm) ;
		if(isRealStressField(f1))
			ret1 = (Vector) (parent->getBehaviour()->getTensor(p_, parent) * ret2) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		else
			ret1 = (Vector) (parent->getBehaviour()->param * ret2) - getParent()->getBehaviour()->getImposedStress(p_, parent) ;
		Vector rate(0., ret2.size()) ;
		if(f2 == STRAIN_FIELD)
			this->getField( STRAIN_RATE_FIELD, p, rate, local, vm) ;
		else
			this->getField( NON_ENRICHED_STRAIN_RATE_FIELD, p, rate, local, vm) ;
		ret1 += (Vector) (dynamic_cast<KelvinVoight *>(parent->getBehaviour())->eta * rate) ;
		if(cleanup) delete vm ;
		return ;
	}
	
	
	this->getField(f1, p, ret1, local, vm, i) ;
	this->getField(f2, p, ret2, local, vm, j) ;
	if(cleanup) delete vm ;
 
}

// void KelvinVoightSpaceTimeElementState::getFieldAtNodes( FieldType f1, FieldType f2, Vector & ret1, Vector & ret2, int i, int j) 
// {	
// 	ElementState::getFieldAtNodes(f1, ret1, i) ;  
// 	ElementState::getFieldAtNodes(f2, ret2, j) ;  
// }

Vector Form::getForcesFromAppliedStress( const Vector & data, Function & shape, const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic, const Vector & normal )
{
	if(isVolumic)
		return VirtualMachine().ieval(Gradient(shape) * data, gp, Jinv, v)  ;
	Vector ret(0., normal.size()) ;
		
	if(normal.size() == 2)
	{
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			ret[0] += (data[0]*normal[0]+data[2]*normal[1])*gp.gaussPoints[i].second ;
			ret[1] += (data[2]*normal[0]+data[1]*normal[1])*gp.gaussPoints[i].second ;
		}
	}
	else
	{
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			ret[0] += (data[0]*normal[0]+data[3]*normal[1]+data[4]*normal[2])*gp.gaussPoints[i].second ;
			ret[1] += (data[1]*normal[1]+data[3]*normal[0]+data[5]*normal[2])*gp.gaussPoints[i].second ;
			ret[2] += (data[2]*normal[2]+data[4]*normal[0]+data[5]*normal[1])*gp.gaussPoints[i].second ;
		}
	}
	return ret * VirtualMachine().ieval(shape, gp) ;
}

Vector Form::getForcesFromAppliedStress( const Function & data, size_t index, size_t externaldofs,  Function & shape, IntegrableEntity * e,const GaussPointArray & gp, const std::valarray<Matrix> & Jinv, std::vector<Variable> & v, bool isVolumic , const Vector& normal) 
{

	VirtualMachine vm ;
	
	size_t n = e->getBoundingPoints().size() ;
	Vector field(0., n*externaldofs) ;
	for(size_t i = 0 ; i < n ; i++)
		field[ i*externaldofs + index ] = vm.eval( data, e->getBoundingPoint(i) ) ;
	
	std::vector<Vector> g(gp.gaussPoints.size(), Vector(0., externaldofs)) ;
	e->getState().getExternalFieldAtGaussPoints( field, externaldofs, g) ;

	if(isVolumic)
		return vm.ieval( Gradient( shape ) * g, gp, Jinv, v) ;
	
	Vector ret(0., normal.size()) ;
	if(normal.size() == 2)
	{
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			ret[0] += (g[0][i]*normal[0]+g[2][i]*normal[1])*gp.gaussPoints[i].second ;
			ret[1] += (g[1][i]*normal[1]+g[2][i]*normal[1])*gp.gaussPoints[i].second ;
		}
	}
	else
	{
		for(size_t i = 0 ; i < gp.gaussPoints.size() ; i++)
		{
			ret[0] += (g[0][i]*normal[0]+g[3][i]*normal[1]+g[4][i]*normal[2])*gp.gaussPoints[i].second ;
			ret[1] += (g[1][i]*normal[1]+g[3][i]*normal[0]+g[5][i]*normal[2])*gp.gaussPoints[i].second ;
			ret[2] += (g[2][i]*normal[2]+g[4][i]*normal[0]+g[5][i]*normal[1])*gp.gaussPoints[i].second ;
		}
	}
	return ret * VirtualMachine().ieval(shape, gp) ;
}

void IntegrableEntity::setCachedGaussPoints(GaussPointArray * gp)
{
	if(cachedGps && cachedGps != gp)
		delete cachedGps ;
	cachedGps = gp ;
}

