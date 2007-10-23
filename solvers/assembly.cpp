//
// C++ Implementation: assembly
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "assembly.h"
#include "gausseidell.h"
#include "polakribiereconjugategradient.h"
#include "biconjugategradientstabilized.h"
using namespace Mu ;



LagrangeMultiplier::LagrangeMultiplier(std::valarray<unsigned int> i, Vector c, double b, int my_id ) : ids(i), coefs(c), id(my_id), value(b), type(GENERAL)
{
} ;

CoordinateIndexedIncompleteSparseMatrix LagrangeMultiplier::getMatrix() const
{
	return CoordinateIndexedIncompleteSparseMatrix( ids, coefs, id) ;
}
void LagrangeMultiplier::setId(const int i)
{
	id = i ;
}


int LagrangeMultiplier::getId() const 
{
	return id ;
}

double LagrangeMultiplier::getValue() const 
{
	return value ;
}


const std::valarray<unsigned int> LagrangeMultiplier::getDofIds() const 
{
	return ids ;
}

std::vector<std::pair<unsigned int, double> > LagrangeMultiplier::getHints() const
{
	return hints ;
}

void LagrangeMultiplier::setHint(std::pair<unsigned int, double> h)
{
	hints.push_back(h) ;
}




Assembly::Assembly() 
{
	this->coordinateIndexedMatrix = NULL ;
	this->nonLinearPartialMatrix = NULL;
	multiplier_offset = 0 ;
	this->displacements.resize(0) ;
	this->externalForces.resize(0) ;
	this->boundaryMatrix = NULL ;
// 	multiplier_offset = 2 ;//bookmark...chk if =3
}

Assembly::~Assembly()
{
	delete this->coordinateIndexedMatrix ;
	delete this->nonLinearPartialMatrix ;
}



Vector & Assembly::getForces()
{
	if(this->coordinateIndexedMatrix == NULL)
		make_final() ;
	
	return this->externalForces ;
}

Vector & Assembly::getNonLinearForces()
{
	if(this->coordinateIndexedMatrix == NULL)
		make_final() ;
	
	return this->nonLinearExternalForces ;
}

void Assembly::operator +=( ElementarySurface * e)
{
	this->add(e) ;
}

void Assembly::operator +=( ElementaryVolume * e)
{
	this->add(e) ;
}


std::vector<ElementarySurface *> Assembly::getElements2d() const
{
	return element2d ;
}

std::vector<ElementarySurface *> & Assembly::getElements2d()
{
	return element2d ;
}

ElementarySurface * Assembly::getElement2d(const size_t i) const
{
	return element2d[i];

}

ElementarySurface * Assembly::getElement2d(const size_t i)
{
	return element2d[i] ;
}

std::vector<ElementaryVolume *> Assembly::getElements3d() const
{
	return element3d ;
}

std::vector<ElementaryVolume *> & Assembly::getElements3d()
{
	return element3d ;
}

ElementaryVolume * Assembly::getElement3d(const size_t i) const
{
	return element3d[i] ;
}

ElementaryVolume * Assembly::getElement3d(const size_t i)
{
	return element3d[i] ;
}


void Assembly::add(ElementarySurface * e)
{
	this->has3Dims = false ;
	std::vector<size_t> ids  = e->getDofIds() ;
	std::sort(ids.begin(), ids.end()) ;
	size_t ndof = e->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	multiplier_offset = std::max(multiplier_offset-ndof, *ids.rbegin()*ndof) + ndof;
	element2d.push_back(e) ;
}
void Mu::Assembly::add(Mu::ElementaryVolume * e)
{
	this->has3Dims = true ;
	std::vector<size_t> ids  = e->getDofIds() ;
	std::sort(ids.begin(), ids.end()) ;
	size_t ndof = e->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	multiplier_offset = std::max(multiplier_offset-ndof, *ids.rbegin()*ndof) + ndof;
	element3d.push_back(e) ;
}

bool Assembly::nonLinearStep()
{

// 	if(displacements.size() == 0)
// 		return ;

	bool nl = false ;
	std::set<std::pair<size_t, size_t> > nonlinmap ;
	
	for(size_t i = 0 ; i < element2d.size() ; i++)
	{
// 		if(i%100 == 0)
// 			std::cerr << "\r computing sparsness pattern... triangle " << i+1 << "/" << element2d.size() << std::flush ;
		if(element2d[i]->getNonLinearBehaviour() != NULL)
		{
			element2d[i]->nonLinearStep(0, &this->displacements) ;
			nl = true ;
		}
		
		if(element2d[i]->getNonLinearBehaviour() != NULL && element2d[i]->getNonLinearBehaviour()->hasInducedMatrix() )
		{
			
			if(element2d[i]->getNonLinearBehaviour()->isActive())
			{
				
				std::vector<size_t> ids = element2d[i]->getDofIds() ;
				for(size_t j = 0 ; j< ids.size() ;j++)
				{
					for(size_t k = 0 ; k< ids.size() ;k++)
					{
						nonlinmap.insert(std::make_pair(ids[j]*2, ids[k]*2)) ;
						nonlinmap.insert(std::make_pair(ids[j]*2+1, ids[k]*2)) ;
						nonlinmap.insert(std::make_pair(ids[j]*2, ids[k]*2+1)) ;
						nonlinmap.insert(std::make_pair(ids[j]*2+1, ids[k]*2+1)) ;
					}
				}
			}
		}
	}
	
	std::valarray<unsigned int> nonlin_column_index((unsigned int)0, nonlinmap.size()) ;
	std::valarray<unsigned int> nonlin_row_index((unsigned int)0, nonlinmap.size()) ;
	size_t current = 0 ;
	for(std::set<std::pair<size_t, size_t> >::const_iterator i = nonlinmap.begin() ; i != nonlinmap.end() ; ++i)
	{
		nonlin_row_index[current] = i->first ;
		nonlin_column_index[current] = i->second ;
		current++ ;
	}
	
	
	delete this->nonLinearPartialMatrix ;
	this->nonLinearPartialMatrix = new CoordinateIndexedIncompleteSparseMatrix(nonlin_row_index, nonlin_column_index) ;
	if(this->nonLinearExternalForces.size() != this->externalForces.size())
	{
		this->nonLinearExternalForces.resize(this->externalForces.size()) ;
		this->nonLinearExternalForces = 0 ;
	}
	
	this->nonLinearExternalForces = 0 ;
	
	for(size_t i = 0 ; i < element2d.size() ; i++)
	{
		if(element2d[i]->getNonLinearBehaviour() != NULL && element2d[i]->getNonLinearBehaviour()->hasInducedMatrix())
		{
// 			if(i%100 == 0)
// 				std::cerr << "\r computing stiffness matrix... triangle " << i+1 << "/" << element2d.size() << std::flush ;
			if(element2d[i]->getNonLinearBehaviour()->isActive())
			{
				std::vector<size_t> ids = element2d[i]->getDofIds() ;
				std::vector<std::vector<Matrix > > mother  = element2d[i]->getNonLinearElementaryMatrix();
				for(size_t j = 0 ; j < ids.size() ;j++)
				{
					getNonLinearMatrix()[ids[j]*2][ids[j]*2] += mother[j][j][0][0] ;
					getNonLinearMatrix()[ids[j]*2][ids[j]*2+1] += mother[j][j][0][1] ;
					getNonLinearMatrix()[ids[j]*2+1][ids[j]*2] += mother[j][j][1][0] ;
					getNonLinearMatrix()[ids[j]*2+1][ids[j]*2+1] += mother[j][j][1][1] ;
					for(size_t k = j+1 ; k < ids.size() ;k++)
					{
						getNonLinearMatrix()[ids[j]*2][ids[k]*2] += mother[j][k][0][0] ;
						getNonLinearMatrix()[ids[j]*2][ids[k]*2+1] += mother[j][k][0][1] ;
						getNonLinearMatrix()[ids[j]*2+1][ids[k]*2] += mother[j][k][1][0] ;
						getNonLinearMatrix()[ids[j]*2+1][ids[k]*2+1] += mother[j][k][1][1] ;
						
						getNonLinearMatrix()[ids[k]*2][ids[j]*2] += mother[k][j][0][0] ;
						getNonLinearMatrix()[ids[k]*2][ids[j]*2+1] += mother[k][j][0][1] ;
						getNonLinearMatrix()[ids[k]*2+1][ids[j]*2] += mother[k][j][1][0] ;
						getNonLinearMatrix()[ids[k]*2+1][ids[j]*2+1] += mother[k][j][1][1] ;
					}
				}
			}
		}
		
		if(element2d[i]->getNonLinearBehaviour() != NULL && element2d[i]->getNonLinearBehaviour()->hasInducedForces())
		{
			if(element2d[i]->getNonLinearBehaviour()->isActive())
			{
				std::vector<size_t> ids = element2d[i]->getDofIds() ;

				Vector forces = element2d[i]->getNonLinearForces() ;
				nl = true ;
				for(size_t j = 0 ; j < ids.size() ;j++)
				{
// 					std::cerr << forces[j*2] << ", " << forces[j*2+1] << std::endl ;
					
					bool hasBC = false;
					for(size_t k = 0 ; k < multipliers.size() ; k++)
					{
						if( (size_t)multipliers[k].getId() == ids[j]*2 ||
						    (size_t)multipliers[k].getId() == ids[j]*2 + 1)
						{
							hasBC = true ;
							break ;
						}
					}	
					
					if(!hasBC)
					{
						this->nonLinearExternalForces[ids[j]*2] += /*this->nonLinearExternalForces[ids[j]*2]*0.2+*/forces[j*2]/**0.8*/ ;
						this->nonLinearExternalForces[ids[j]*2+1] += /*this->nonLinearExternalForces[ids[j]*2+1]*0.2+*/forces[j*2+1]/**0.8*/ ;
					}
				}
			}
		}
	}
	
	return nl ;
// 	std::cerr << " ...done" << std::endl ;
	
}


void Assembly::setBoundaryConditions()
{

	this->externalForces.resize(coordinateIndexedMatrix->row_size.size()) ;
	this->externalForces = 0 ;
	this->nonLinearExternalForces.resize(coordinateIndexedMatrix->row_size.size()) ;
	this->nonLinearExternalForces = 0 ;
	
	size_t ndofs = 0 ;
	if(!element2d.empty())
		ndofs = element2d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	
	std::cerr << " setting BCs... forces element " << 0 << "/" << std::max(element2d.size(),element3d.size()) << std::flush ;
	for(size_t i = 0 ; i < element2d.size() ; i++)
	{
		if(i%1000 == 0)
			std::cerr << "\r setting BCs... forces element " << i << "/" << element2d.size() << std::flush ;
		std::vector<size_t> ids = element2d[i]->getDofIds() ;
// 		std::vector<std::vector<Matrix > > mother  = element2d[i]->getElementaryMatrix();
		
		if(element2d[i]->getBehaviour()->hasInducedForces())
		{
			Vector f = element2d[i]->getForces() ;
			
			for(size_t j = 0 ; j < ids.size() ; j++)
			{
				for(size_t l = 0 ; l < ndofs  ; l++)
				{
					this->externalForces[ndofs*ids[j]+l] += f[ndofs*j+l] ;
				}
			}
		}
	}
	
	
	
	if(!element3d.empty())
		ndofs = element3d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;
	
	for(size_t i = 0 ; i < element3d.size() ; i++)
	{
		if(i%1000 == 0)
			std::cerr << "\r setting BCs... forces element " << i << "/" << element3d.size() << std::flush ;
		
		std::vector<size_t> ids = element3d[i]->getDofIds() ;
// 		std::vector<std::vector<Matrix > > mother  = element3d[i]->getElementaryMatrix();
		
		if(element3d[i]->getBehaviour()->hasInducedForces())
		{
			Vector f = element3d[i]->getForces() ;
			
			for(size_t j = 0 ; j < ids.size() ;j++)
			{
				for(size_t l = 0 ; l < ndofs  ; l++)
				{
					this->externalForces[ndofs*ids[j]+l] += f[ndofs*j+l] ;
				}
			}
		}
	}
	
	std::cerr << " ...done" << std::endl ;
	
	for(size_t i = 0 ; i < multipliers.size() ; i++)
	{
		if(multipliers[i].type == GENERAL)
		{
			getMatrix()+=multipliers[i].getMatrix() ;
		}
		else if(multipliers[i].type == SET_FORCE_XI 
				|| multipliers[i].type == SET_FORCE_ETA
				|| multipliers[i].type == SET_FORCE_ZETA)
		{
			this->externalForces[multipliers[i].getId()] += multipliers[i].getValue() ; 
		}
	}

	
	
	//we might have voids in the numbering...
	for(size_t i = 0 ; i < coordinateIndexedMatrix->row_size.size() ; i++)
	{
		if(getMatrix()[i][i] == 0)
			getMatrix()[i][i] = 1 ;
	}
	
	std::sort(multipliers.begin(), multipliers.end()) ;
	
	std::cerr << " setting BCs... displacement dof " << 0 << "/" << coordinateIndexedMatrix->row_size.size() << std::flush ;
	
	for(size_t k = 0 ; k < coordinateIndexedMatrix->row_size.size() ; k++)
	{
		if(k% 1000 == 0)
			std::cerr << "\r setting BCs... displacement dof " << k << "/" << coordinateIndexedMatrix->row_size.size() << std::flush ;

		int array_index = coordinateIndexedMatrix->accumulated_row_size[k] ;

		for(int m = 0 ;  m < (int)multipliers.size() ; m++)
		{
			int id = multipliers[m].getId() ;
			double val = getMatrix()[k][id] ;
			
			if(id == (int)k && multipliers[m].type != SET_FORCE_XI 
				&&  multipliers[m].type != SET_FORCE_ETA 
				&&  multipliers[m].type != SET_FORCE_ZETA 
				&& multipliers[m].type != GENERAL
			  )
			{
				this->externalForces[id] = multipliers[m].getValue() ;
				

				for(size_t l = 0 ; l < coordinateIndexedMatrix->row_size[k] ; l++)
				{
					if((int)coordinateIndexedMatrix->column_index[array_index+l] == id)
					{
						coordinateIndexedMatrix->array[array_index+l] = 1 ;
					}
					else
					{
						coordinateIndexedMatrix->array[array_index+l] = 0 ;
					}
				}

				
				break ;
				
			}
			else if( multipliers[m].type != SET_FORCE_XI 
						&&  multipliers[m].type != SET_FORCE_ETA 
						&&  multipliers[m].type != SET_FORCE_ZETA 
						&& multipliers[m].type != GENERAL
			       )
			{
				this->externalForces[k] -= multipliers[m].getValue()*val ;
				getMatrix()[k][id] = 0 ;
			}
		}
	}
	
	std::cerr << " ...done" << std::endl ;
	
	multipliers.clear() ;
	
// 	std::cerr << " ...done." << std::endl ;

}

void Assembly::make_final()
{
	if (has3Dims == false)
	{
		if( element2d.empty())
		{
			std::cerr << "no elements in mesh !" << std::endl ;
			exit(0) ;
		}
		size_t ndof = element2d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ; 
		
		std::cerr << "ndof = " << ndof << std::endl ;
		
		delete this->coordinateIndexedMatrix ;

		std::set<std::pair<unsigned int, unsigned int> > * map  = new std::set<std::pair<unsigned int, unsigned int> >();
		
		for(size_t i = 0 ; i < element2d.size() ; i++)
		{
			if(i%1000 == 0)
				std::cerr << "\r computing sparsness pattern... triangle " << i+1 << "/" << element2d.size() << std::flush ;
			std::vector<size_t> ids = element2d[i]->getDofIds() ;
			
			for(size_t j = 0 ; j< ids.size() ;j++)
			{
				for(size_t n = 0 ; n < ndof ; n++)
				{
					for(size_t m = 0 ; m < ndof ; m++)
					{
						map->insert(std::make_pair(ids[j]*ndof+n, ids[j]*ndof+m)) ;
					}
				}

				for(size_t k = j+1 ; k< ids.size() ;k++)
				{
					for(size_t n = 0 ; n < ndof ; n++)
					{
						for(size_t m = 0 ; m < ndof ; m++)
						{
							map->insert(std::make_pair(ids[j]*ndof+n, ids[k]*ndof+m)) ;
							map->insert(std::make_pair(ids[k]*ndof+n, ids[j]*ndof+m)) ;
						}
					}
				}
			}
		}
		
		size_t max = map->rbegin()->first +1;
		size_t realDofs = max ;
		for(size_t i = 0 ; i < multipliers.size() ; i++)
		{
			if(multipliers[i].type == GENERAL)
			{
				std::valarray<unsigned int> ids = multipliers[i].getDofIds() ;
				for(size_t j = 0 ; j< ids.size() ;j++)
				{
					map->insert(std::make_pair( ids[j], multipliers[i].getId())) ;
					map->insert(std::make_pair( multipliers[i].getId(), ids[j])) ;
				}
				max++ ;
			}
		}
		
		//filling eventual gaps
		for(size_t i = 0 ; i < realDofs ; i++)
		{
			map->insert(std::make_pair(i,i)) ;
		}
		
		std::cerr << " ...done" << std::endl ;
		size_t total_num_dof = map->rbegin()->first+1 ;
		nonLinearExternalForces.resize(total_num_dof, 0.) ;
		
		std::valarray<unsigned int> row_length((unsigned int)0, total_num_dof) ;
		std::valarray<unsigned int> column_index((unsigned int)0, map->size()) ;
		
		size_t current = 0 ;
		for(std::set<std::pair<unsigned int, unsigned int> >::const_iterator i = map->begin() ; i != map->end() ; ++i)
		{
			row_length[i->first]++ ;
			column_index[current] = i->second ;
			current++ ;
		}
		
		delete map ;
		this->coordinateIndexedMatrix = new CoordinateIndexedSparseMatrix(row_length, column_index) ;

		if(this->displacements.size() != max)
		{
			this->displacements.resize(max) ;
			this->displacements = 0 ;
		}
		
		for(size_t i = 0 ; i < element2d.size() ; i++)
		{
			if(i%100 == 0)
				std::cerr << "\r computing stiffness matrix... triangle " << i+1 << "/" << element2d.size() << std::flush ;
			
			std::vector<size_t> ids = element2d[i]->getDofIds() ;
			std::vector<std::vector<Matrix > > mother  = element2d[i]->getElementaryMatrix();
			
			for(size_t j = 0 ; j < ids.size() ;j++)
			{
				for(size_t n = 0 ; n < ndof ; n++)
				{
					for(size_t m = 0 ; m < ndof ; m++)
					{
						getMatrix()[ids[j]*ndof+n][ids[j]*ndof+m] += mother[j][j][n][m] ;
					}
				}
	
				for(size_t k = j+1 ; k < ids.size() ;k++)
				{
					for(size_t n = 0 ; n < ndof ; n++)
					{
						for(size_t m = 0 ; m < ndof ; m++)
						{
							getMatrix()[ids[j]*ndof+n][ids[k]*ndof+m] += mother[j][k][n][m] ;
							getMatrix()[ids[k]*ndof+n][ids[j]*ndof+m] += mother[k][j][n][m] ;
						}
					}
				}
			}
		}
		
		std::cerr << " ...done" << std::endl ;
		
		setBoundaryConditions() ;
		
	}
	else 
	{
		std::set<std::pair<unsigned int, unsigned int> > * map  = new std::set<std::pair<unsigned int, unsigned int> >();
		size_t ndof = 0 ;
		if(!element3d.empty())
			ndof = element3d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ; 
			
		for(size_t i = 0 ; i < element3d.size() ; i++)
		{
			if(i%1000 == 0)
				std::cerr << "\r computing sparsness pattern... tetrahedron " << i+1 << "/" << element3d.size() << std::flush ;
			std::vector<size_t> ids = element3d[i]->getDofIds() ;
			
			for(size_t j = 0 ; j< ids.size() ;j++)
			{
				for(size_t n = 0 ; n < ndof ; n++)
				{
					for(size_t m = 0 ; m < ndof ; m++)
					{
						map->insert(std::make_pair(ids[j]*ndof+n, ids[j]*ndof+m)) ;
					}
				}
				
				
				for(size_t k = j+1 ; k< ids.size() ;k++)
				{
					
					for(size_t n = 0 ; n < ndof ; n++)
					{
						for(size_t m = 0 ; m < ndof ; m++)
						{
							map->insert(std::make_pair(ids[j]*ndof+n, ids[k]*ndof+m)) ;
							map->insert(std::make_pair(ids[k]*ndof+n, ids[j]*ndof+m)) ;
						}
					}
				}
			}
		}
			
			size_t max = map->rbegin()->first +1;
			size_t realDofs = max ;
			for(size_t i = 0 ; i < multipliers.size() ; i++)
			{
				if(multipliers[i].type == GENERAL)
				{
					std::valarray<unsigned int> ids = multipliers[i].getDofIds() ;
					for(size_t j = 0 ; j< ids.size() ;j++)
					{
						map->insert(std::make_pair( ids[j], multipliers[i].getId())) ;
						map->insert(std::make_pair( multipliers[i].getId(), ids[j])) ;
					}
					max++ ;
				}
			}
			
	//filling eventual gaps
			for(size_t i = 0 ; i < realDofs ; i++)
			{
				map->insert(std::make_pair(i,i)) ;
			}
			
			size_t total_num_dof = map->rbegin()->first+1 ;
			nonLinearExternalForces.resize(total_num_dof, 0.) ;
			std::valarray<unsigned int> row_length((unsigned int)0, total_num_dof) ;
			std::valarray<unsigned int> column_index((unsigned int)0, map->size()) ;
			size_t current = 0 ;
			for(std::set<std::pair<unsigned int, unsigned int> >::const_iterator i = map->begin() ; i != map->end() ; ++i)
			{
		
				
				row_length[i->first]++ ;
				column_index[current] = i->second ;
				current++ ;
			}
			delete map ;
			this->coordinateIndexedMatrix = new CoordinateIndexedSparseMatrix(row_length, column_index) ;
			this->displacements.resize(max, 0.) ;
			std::cerr << " ...done" << std::endl ;
			for(size_t i = 0 ; i < element3d.size() ; i++)
			{
				if(i%1000 == 0)
					std::cerr << "\r computing stiffness matrix... tetrahedron " << i+1 << "/" << element3d.size() << std::flush ;
				
				std::vector<size_t> ids = element3d[i]->getDofIds() ;
				std::vector<std::vector<Matrix > > mother  = element3d[i]->getElementaryMatrix();
				for(size_t j = 0 ; j < ids.size() ;j++)
				{
					ids[j] *= ndof ;
				}
				
				for(size_t j = 0 ; j < ids.size() ;j++)
				{

					for(size_t l = 0 ; l < ndof  ; l++)
					{
						for(size_t m = 0 ; m < ndof  ; m++)
						{
							getMatrix()[ids[j]+l][ids[j]+m] += mother[j][j][l][m] ;
						}
					}
					
					for(size_t k = j+1 ; k < ids.size() ;k++)
					{
						for(size_t l = 0 ; l < ndof  ; l++)
						{
							for(size_t m = 0 ; m < ndof  ; m++)
							{
								getMatrix()[ids[j]+l][ids[k]+m] += mother[j][k][l][m] ;
								getMatrix()[ids[k]+l][ids[j]+m] += mother[k][j][l][m] ;
							}
						}
					}
				}
			}
			
			std::cerr << " ...done" << std::endl ;
			
		setBoundaryConditions() ;
	
}
}

CoordinateIndexedSparseMatrix & Assembly::getMatrix()
{
	return *this->coordinateIndexedMatrix ;
}

const CoordinateIndexedSparseMatrix & Assembly::getMatrix() const
{
	return *this->coordinateIndexedMatrix ;
}

CoordinateIndexedIncompleteSparseMatrix & Assembly::getNonLinearMatrix() 
{
	return *this->nonLinearPartialMatrix ;
}


void Assembly::print() 
{
	if(this->coordinateIndexedMatrix == NULL)
		make_final() ;
	
	for(size_t i = 0 ; i < getMatrix().row_size.size() ; i++)
	{
		for(size_t j = 0 ; j < getMatrix().row_size.size() ; j++)
			std::cerr << getMatrix()[i][j] << "   " << std::flush ;
		
		std::cerr << std::endl ;
	}
	std::cerr << std::endl ;
	for(size_t i = 0 ; i < externalForces.size() ; i++)
	{
		std::cerr << externalForces[i] << std::endl ;
	}
}

double Assembly::froebeniusNorm()
{
	return sqrt(std::inner_product(&coordinateIndexedMatrix->array[0], 
	                               &coordinateIndexedMatrix->array[coordinateIndexedMatrix->array.size()], 
	                               &coordinateIndexedMatrix->array[0], (double)(0.))) ;
}

void  Assembly::fixPoint(size_t id)
{
	setPoint(0,0, id) ;
}


void Assembly::setPoint(double ex, size_t id)
{		
	Vector c(2) ;
	std::valarray<unsigned int> i(2) ;
	if(multipliers.empty() || std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id)) == multipliers.end() )
	{
		multipliers.push_back(LagrangeMultiplier(i,c, ex, id)) ;
		multipliers.rbegin()->type = SET_ALONG_XI ;
	}
	
	return ;
}

void Assembly::setPoint(double ex, double ey, size_t id)
{		
	Vector c(2) ;
	std::valarray<unsigned int> i(2) ;
	if(multipliers.empty() || std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2)) == multipliers.end() )
	{
		multipliers.push_back(LagrangeMultiplier(i,c, ex, id*2)) ;
		multipliers.rbegin()->type = SET_ALONG_XI ;
	}
	
	if(multipliers.empty() || std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2+1)) == multipliers.end())
	{
		multipliers.push_back(LagrangeMultiplier(i,c, ey, id*2+1)) ;
		multipliers.rbegin()->type = SET_ALONG_ETA ;

	}


	return ;
}

void Assembly::setPoint(double ex, double ey, double ez, size_t id)
{
	Vector c(2) ;

	std::valarray<unsigned int> i(2) ;
	
	if(std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*3)) == multipliers.end())
	{
		multipliers.push_back(LagrangeMultiplier(i,c, ex, id*3)) ;
		multipliers.rbegin()->type = SET_ALONG_XI ;
	}
	if(std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*3+1)) == multipliers.end())
	{
		multipliers.push_back(LagrangeMultiplier(i,c, ey, id*3+1)) ;
		multipliers.rbegin()->type = SET_ALONG_ETA ;
	}
	if(std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*3+2)) == multipliers.end())
	{
		multipliers.push_back(LagrangeMultiplier(i,c, ez, id*3+2)) ;
		multipliers.rbegin()->type = SET_ALONG_ZETA ;
	}
	return ;
}

void Assembly::setPeriodicPoint(size_t id0, size_t id1)
{
	if( has3Dims == false)
	{
		Vector c(2) ;
		c[0] = 1 ;
		c[1] = -1 ;
		std::valarray<unsigned int> i(2) ;
		i[0] = id0*2 ;
		i[1] = id1*2 ;
		multipliers.push_back(LagrangeMultiplier(i,c, 0.)) ;
		i[0] = id0*2 + 1;
		i[1] = id1*2 + 1;
		multipliers.push_back(LagrangeMultiplier(i,c, 0.)) ;
	}
	else
	{
		Vector c(2) ;
		c[0] = 1 ;
		c[1] = -1 ;
		std::valarray<unsigned int> i(2) ;
		i[0] = id0*2 ;
		i[1] = id1*2 ;
		multipliers.push_back(LagrangeMultiplier(i,c, 0.)) ;
		i[0] = id0*2 + 1;
		i[1] = id1*2 + 1;
		multipliers.push_back(LagrangeMultiplier(i,c, 0.)) ;
		i[0] = id0*3 + 1;
		i[1] = id1*3 + 1;
		multipliers.push_back(LagrangeMultiplier(i,c, 0.)) ;
	}
	//change ..push bac again
}


void Assembly::set3D()
{
	has3Dims = true ;
}

void Assembly::set2D()
{
	has3Dims = false ;
}

void Assembly::clear()
{
	element2d.clear() ;
	element3d.clear() ;
	
	delete coordinateIndexedMatrix ;
	coordinateIndexedMatrix = NULL ;
	delete nonLinearPartialMatrix ;
	nonLinearPartialMatrix = NULL ;
	delete boundaryMatrix ;
	boundaryMatrix = NULL ;
}

void Assembly::setForceOn(double val, size_t id)
{
	std::valarray<unsigned int> i(2) ;
	Vector c(2) ;
	multipliers.push_back(LagrangeMultiplier(i,c,val, id)) ;
	multipliers.rbegin()->type = SET_FORCE_XI ;
}

void Assembly::setForceOn(Variable v, double val, size_t id)
{
	std::valarray<unsigned int> i(2) ;
	Vector c(2) ;
	if(!has3Dims)
	{
		switch(v)
		{
		case XI:
			{
// 				if(std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2)) == multipliers.end())
// 				{
					multipliers.push_back(LagrangeMultiplier(i,c,val, id*2)) ;
					multipliers.rbegin()->type = SET_FORCE_XI ;
// 				}
				break ;
			}
		case ETA:
			{
// 				if(std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2+1)) == multipliers.end())
// 				{
					multipliers.push_back(LagrangeMultiplier(i,c,val, id*2+1)) ;
					multipliers.rbegin()->type = SET_FORCE_ETA ;
// 				}
				break ;
			}
		default:
			{
				break ;
			}
		}
	}
	else
	{
		switch(v)
		{
		case XI:
			{
				multipliers.push_back(LagrangeMultiplier(i,c,val, id*3)) ;
				multipliers.rbegin()->type = SET_FORCE_XI ;
				break ;
			}
		case ETA:
			{
				multipliers.push_back(LagrangeMultiplier(i,c,val, id*3+1)) ;
				multipliers.rbegin()->type = SET_FORCE_ETA ;
				break ;
			}
		case ZETA:
			{
				multipliers.push_back(LagrangeMultiplier(i,c,val, id*3+2)) ;
				multipliers.rbegin()->type = SET_FORCE_ZETA ;
				break ;
			}
		default:
			{
				break ;
			}
		}
	}
	return ;
}

void Assembly::setPointAlong(Variable v, double val, size_t id) 
{
	std::valarray<unsigned int> i(2) ;
	Vector c(2) ;
	if(has3Dims == false)
	{
		switch(v)
		{
			case XI:
			{
				if(multipliers.empty() ||std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2)) == multipliers.end())
				{
					multipliers.push_back(LagrangeMultiplier(i,c,val, id*2)) ;
					multipliers.rbegin()->type = SET_ALONG_XI ;
				}
				break ;
			}
			case ETA:
			{
				if(multipliers.empty() ||std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2+1)) == multipliers.end())
				{
					multipliers.push_back(LagrangeMultiplier(i,c,val, id*2+1)) ;
					multipliers.rbegin()->type = SET_ALONG_ETA ;
				}
				break ;
			}
		default:
			{
				break ;
			}
		}
	}
	else
	{
		switch(v)
		{
		case XI:
			{
				multipliers.push_back(LagrangeMultiplier(i,c,val, id*3)) ;
				multipliers.rbegin()->type = SET_ALONG_XI ;
				break ;
			}
		case ETA:
			{
				multipliers.push_back(LagrangeMultiplier(i,c,val, id*3+1)) ;
				multipliers.rbegin()->type = SET_ALONG_ETA ;
				break ;
			}
		case ZETA:
			{
				multipliers.push_back(LagrangeMultiplier(i,c,val, id*3+2)) ;
				multipliers.rbegin()->type = SET_ALONG_ETA ;
				break ;
			}
		default:
			{
				break ;
			}
		}
	}
	return ;
	
}


void Assembly::fixPoint(size_t id, Mu::Variable v)
{
	setPointAlong(v, 0, id) ;
}

Vector & Assembly::solve(Vector x0, size_t maxit, const bool verbose)
{
	if(this->coordinateIndexedMatrix == NULL)
		make_final() ;
	
	return GaussSeidel(getMatrix(), externalForces).solve(x0, NULL) ;
}

Vector & Assembly::cgsolve(Vector x0, size_t maxit) 
{
	if(this->coordinateIndexedMatrix == NULL)
		make_final() ;
	
// 	print() ;
	
	if(x0.size() == 0)
	{
// 		displacements = ConjugateGradient(getMatrix(), externalForces).solve(displacements, NULL,1e-12, 16000, true) ;
		displacements = BiConjugateGradientStabilized(getMatrix(), externalForces).solve(displacements, NULL,1e-12, 16000, true) ;
// 		ConjugateGradientWithSecant(this).solve() ;
	}
	else
	{
		displacements = BiConjugateGradientStabilized(getMatrix(), externalForces).solve(displacements, NULL,1e-12, 16000, true) ;
// 		displacements = ConjugateGradient(getMatrix(), externalForces).solve(x0, NULL,1e-12, 16000, true) ;
// 		ConjugateGradientWithSecant(this).solve();
	}
	return displacements ;
// 	return ConjugateGradientWithSecant(this).solve() ;
	
}


void Assembly::fix()
{
	make_final() ;
}

Vector Assembly::cgnpsolve(const Vector b, size_t maxit) 
{
	NullPreconditionner np ;
	return ConjugateGradient(getMatrix(), externalForces).solve(b, &np, 1e-15, -1) ;
}

void Assembly::printDiag() const
{
	for(size_t i = 0 ; i < getMatrix().row_size.size() ; i++)
		std::cerr << getMatrix()[i][i] << std::endl ;
}

Vector operator *(const std::map< std::pair<size_t, size_t>, Matrix > A , const Vector x)
{
	Vector ret(0., x.size()) ;
	
	size_t ddl = A.begin()->second.numRows() ;
	for(std::map< std::pair<size_t, size_t>, Matrix >::const_iterator ij = A.begin() ; ij != A.end() ; ++ij)
	{
		for(size_t i = 0 ;  i < ddl ; i++)
		{
			for(size_t j = 0 ;  j < ddl ; j++)
			{
				ret[ij->first.first*ddl+i] += x[ij->first.second*ddl+j] * ij->second[i][j] ;
			}
		}
	}
	
	return ret ;
}

Vector operator *(const std::map< std::pair<size_t, size_t>, double > A , const Vector x)
{
	Vector ret(0., x.size()) ;
	
	for(std::map< std::pair<size_t, size_t>, double >::const_iterator ij = A.begin() ; ij != A.end() ; ++ij)
	{

		ret[ij->first.first] += x[ij->first.second] * ij->second ;

	}
	
	return ret ;
}


std::map< std::pair<size_t, size_t>, Matrix > operator *(const std::map< std::pair<size_t, size_t>, Matrix > A , const std::map< std::pair<size_t, size_t>, Matrix > B)
{
	std::map< std::pair<size_t, size_t>, Matrix > ret ;
	for(std::map< std::pair<size_t, size_t>, Matrix >::const_iterator ik = A.begin() ; ik != A.end() ; ++ik)
	{
		for(std::map< std::pair<size_t, size_t>, Matrix >::const_iterator kj = B.begin() ; kj != B.end() ; ++kj)
		{
			if(ret.find(std::pair<size_t, size_t>(ik->first.first, kj->first.second)) == ret.end())
			{
				ret[std::pair<size_t, size_t>(ik->first.first, kj->first.second)] = ik->second*kj->second ;
			}
			else
			{
				ret[std::pair<size_t, size_t>(ik->first.first, kj->first.second)] += ik->second*kj->second ;
			}
		}
	}
	
	return ret ;
}


std::map<std::pair<size_t, size_t>, Matrix> incompleteCholeskyDecomposition(std::map<std::pair<size_t, size_t>, Matrix>  & morseMatrix)
{
	
	std::map<std::pair<size_t, size_t>, double> cholesky ;
	std::map<std::pair<size_t, size_t>, Matrix> cholesky_final ;
	
	size_t ddl = morseMatrix.begin()->second.numRows() ;
	
	//first, let us get the empty cholesky structure
	
	for(std::map<std::pair<size_t, size_t>, Matrix>::iterator ij =  morseMatrix.begin() ; ij != morseMatrix.end() ; ++ij)
	{
		if(ij->first.first <= ij->first.second) //if line < column
		{
			cholesky_final[ij->first] = Matrix(ddl,ddl) ;
			
			for(size_t i = 0 ; i < ddl ; i++)
			{
				for(size_t j = 0 ; j < ddl ; j++)
				{
					std::pair<size_t, size_t> kl(ij->first.first*ddl+i, ij->first.second*ddl+j) ;
					cholesky[kl] = 0;
				}
			}
		}
	}
	
	size_t matrix_size  = cholesky.rbegin()->first.first ;

	cholesky[std::pair<size_t, size_t>(0,0)] = sqrt(morseMatrix.begin()->second[0][0]) ;
	
	
	for(size_t k = 0 ; k < matrix_size ; k++)
	{
		for(size_t s = k+1 ; s < matrix_size ; s++)
		{
			std::pair<size_t, size_t> sk(s,k) ;
			std::pair<size_t, size_t> kk(k,k) ;
			if(cholesky.find(sk) != cholesky.end())
			{
				cholesky[sk] = morseMatrix[std::pair<size_t, size_t>(s/ddl, k/ddl)][s%ddl][k%ddl]/cholesky[kk] ;
			}
		}
		
		
		for(size_t j = k+1 ; j < matrix_size ; j++)
		{
			for(size_t i = j ; i < matrix_size ; i++)
			{
				std::pair<size_t, size_t> ik(i,k) ;
				std::pair<size_t, size_t> jk(j,k) ;
				std::pair<size_t, size_t> ij(i,j) ;
				if(cholesky.find(ij) != cholesky.end() && cholesky.find(ik) != cholesky.end() && cholesky.find(jk) != cholesky.end())
					cholesky[ij] +=  (morseMatrix[std::pair<size_t, size_t>(i/ddl, k/ddl)][i%ddl][k%ddl] - cholesky[ik]* cholesky[jk]) ;
			}
		}
		
		std::pair<size_t, size_t> kk(k+1,k+1) ;
		cholesky[kk] = sqrt(morseMatrix[std::pair<size_t, size_t>((k+1)/ddl, (k+1)/ddl)][(k+1)%ddl][(k+1)%ddl] );
	}

	for( std::map<std::pair<size_t, size_t>, double>::iterator ij = cholesky.begin() ;  ij != cholesky.end() ; ++ij)
	{
		size_t i = ij->first.first ;
		size_t j = ij->first.second ;
		cholesky_final[std::pair<size_t, size_t>(i/ddl, j/ddl)][i%ddl][j%ddl]  = ij->second ;
	}
	
	return cholesky_final ;
}

Vector solveCholeskyDecomposedSystem(std::map<std::pair<size_t, size_t>, Matrix> & choleskyDecomposition, const Vector &b)
{
	return Vector() ; //shut up the compiler
}


