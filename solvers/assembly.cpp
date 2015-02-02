//
// C++ Implementation: assembly
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2013
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "assembly.h"
#include "gausseidell.h"
#include "polakribiereconjugategradient.h"
#include "biconjugategradientstabilized.h"
#include "../physics/dual_behaviour.h"
#include "eigenvalues.h"
#include <valarray>
#include <sys/time.h>

using namespace Amie ;

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
    rowstart = 0 ;
    colstart = 0 ;
    ndofmax = 0 ;
    this->coordinateIndexedMatrix = nullptr ;
    this->nonLinearPartialMatrix = nullptr;
    multiplier_offset = 0 ;
    this->displacements.resize(0) ;
    this->externalForces.resize(0) ;
    this->prevDisplacements.resize(0) ;
    this->naturalBoundaryConditionForces.resize(0) ;
    this->boundaryMatrix = nullptr ;
    ndof = 1 ;
    dim = SPACE_THREE_DIMENSIONAL ;
    epsilon = 1e-16 ;
// 	multiplier_offset = 2 ;//bookmark...chk if =3
}

Assembly::~Assembly()
{
    delete this->coordinateIndexedMatrix ;
    delete this->nonLinearPartialMatrix ;
}

Vector & Assembly::getForces()
{
    if(this->coordinateIndexedMatrix == nullptr)
        make_final() ;

    return this->externalForces ;
}

Vector & Assembly::getNaturalBoundaryConditionForces()
{
    if(this->coordinateIndexedMatrix == nullptr)
        make_final() ;

    return this->naturalBoundaryConditionForces ;
}

Vector & Assembly::getNonLinearForces()
{
    if(this->coordinateIndexedMatrix == nullptr)
        make_final() ;

    return this->nonLinearExternalForces ;
}

void Assembly::add(ElementarySurface * e, double scale)
{
    dim = SPACE_TWO_DIMENSIONAL ;
    ndof = e->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    multiplier_offset =  ndof;
    element2d.push_back(e) ;
    scales.push_back(scale);
}

void Assembly::add(ElementaryVolume * e, double scale)
{
    dim = SPACE_THREE_DIMENSIONAL ;
    ndof = e->getBehaviour()->getNumberOfDegreesOfFreedom() ;
    multiplier_offset =  ndof;
    element3d.push_back(e) ;
    scales.push_back(scale);
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
        if(element2d[i]->getNonLinearBehaviour() != nullptr)
        {
            element2d[i]->nonLinearStep(0, &this->displacements) ;
            nl = true ;
        }

        if(element2d[i]->getNonLinearBehaviour() != nullptr && element2d[i]->getNonLinearBehaviour()->hasInducedMatrix() )
        {

            if(element2d[i]->getNonLinearBehaviour()->isActive())
            {

                std::vector<size_t> ids = element2d[i]->getDofIds() ;
                for(size_t j = 0 ; j< ids.size() ; j++)
                {
                    for(size_t k = 0 ; k< ids.size() ; k++)
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
    for(auto i = nonlinmap.begin() ; i != nonlinmap.end() ; ++i)
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
        if(element2d[i]->getNonLinearBehaviour() != nullptr && element2d[i]->getNonLinearBehaviour()->hasInducedMatrix())
        {
// 			if(i%100 == 0)
// 				std::cerr << "\r computing stiffness matrix... triangle " << i+1 << "/" << element2d.size() << std::flush ;
            if(element2d[i]->getNonLinearBehaviour()->isActive())
            {
                std::vector<size_t> ids = element2d[i]->getDofIds() ;
                std::valarray<std::valarray<Matrix > > mother  = element2d[i]->getNonLinearElementaryMatrix();
                for(size_t j = 0 ; j < ids.size() ; j++)
                {
                    getNonLinearMatrix()[ids[j]*2][ids[j]*2] += mother[j][j][0][0] ;
                    getNonLinearMatrix()[ids[j]*2][ids[j]*2+1] += mother[j][j][0][1] ;
                    getNonLinearMatrix()[ids[j]*2+1][ids[j]*2] += mother[j][j][1][0] ;
                    getNonLinearMatrix()[ids[j]*2+1][ids[j]*2+1] += mother[j][j][1][1] ;
                    for(size_t k = j+1 ; k < ids.size() ; k++)
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

        if(element2d[i]->getNonLinearBehaviour() != nullptr )
        {
            if(element2d[i]->getNonLinearBehaviour()->isActive())
            {
                std::vector<size_t> ids = element2d[i]->getDofIds() ;

                Vector forces = element2d[i]->getNonLinearForces() ;
                nl = true ;
                for(size_t j = 0 ; j < ids.size() ; j++)
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
    this->externalForces.resize(coordinateIndexedMatrix->row_size.size()*coordinateIndexedMatrix->stride, 0.) ;
    this->externalForces = 0 ;
    this->naturalBoundaryConditionForces.resize(coordinateIndexedMatrix->row_size.size()*coordinateIndexedMatrix->stride) ;
    this->naturalBoundaryConditionForces = 0 ;
    this->nonLinearExternalForces.resize(coordinateIndexedMatrix->row_size.size()*coordinateIndexedMatrix->stride) ;
    this->nonLinearExternalForces = 0 ;

    std::sort(multipliers.begin(), multipliers.end()) ;
    std::valarray<int> multiplierIds(multipliers.size()) ;
    for(size_t i = 0 ; i < multiplierIds.size() ; i++)
    {
        multiplierIds[i] = multipliers[i].getId() ;
    }



    int stride = coordinateIndexedMatrix->stride ;
    std::cerr << " setting BCs... displacement dof " << 0 << "/" << coordinateIndexedMatrix->row_size.size() << std::flush ;
    for(size_t k = 0 ; k < coordinateIndexedMatrix->row_size.size() ; k++)
    {
        if(k% 1000 == 0)
            std::cerr << "\r setting BCs... displacement dof " << k*stride << "/" << coordinateIndexedMatrix->row_size.size()*stride << std::flush ;
        int lineBlockIndex = k ;
        int * start_multiplier = std::lower_bound(&multiplierIds[0], &multiplierIds[multiplierIds.size()], coordinateIndexedMatrix->column_index[coordinateIndexedMatrix->accumulated_row_size[k]]*stride) ;
        int * end_multiplier = std::upper_bound(start_multiplier,  &multiplierIds[multiplierIds.size()],coordinateIndexedMatrix->column_index[coordinateIndexedMatrix->accumulated_row_size[k]+coordinateIndexedMatrix->row_size[k]-1]*stride+stride-1 ) ;
        int * start_multiplier_in_line = std::lower_bound(start_multiplier, end_multiplier,k*stride) ;
        int * end_multiplier_in_line = std::upper_bound(start_multiplier,  end_multiplier,k*stride+stride-1) ;

// 		int start_multiplier_index = start_multiplier-&multiplierIds[0] ;
// 		int end_multiplier_index =  end_multiplier-&multiplierIds[0] ;
        int start_multiplier_in_line_index = start_multiplier_in_line-&multiplierIds[0] ;
        int end_multiplier_in_line_index =  end_multiplier_in_line-&multiplierIds[0] ;

        for(size_t l = 0 ; l <  coordinateIndexedMatrix->row_size[k] ; l++)
        {
            int columnBlockIndex = coordinateIndexedMatrix->column_index[coordinateIndexedMatrix->accumulated_row_size[k]+l] ;

            int * start_multiplier_in_block = std::lower_bound(start_multiplier, end_multiplier,columnBlockIndex*stride) ;
            int * end_multiplier_in_block = std::upper_bound(start_multiplier,  end_multiplier,columnBlockIndex*stride+stride-1) ;

            int start_multiplier_in_block_index = start_multiplier_in_block-&multiplierIds[0] ;
            int end_multiplier_in_block_index =  end_multiplier_in_block-&multiplierIds[0] ;

            for(int p = start_multiplier_in_line_index ; p != end_multiplier_in_line_index ; p++)
            {
                if(multipliers[p].type != SET_FORCE_XI
                        &&  multipliers[p].type != SET_FORCE_ETA
                        &&  multipliers[p].type != SET_FORCE_ZETA
                        &&  multipliers[p].type != SET_FORCE_INDEXED_AXIS
                        &&  multipliers[p].type != GENERAL)
                {
                    int id = multipliers[p].getId() ;
                    for(int m = 0 ; m < stride ; m++)
                    {
                        if( id != (lineBlockIndex*stride+m))
                        {
                            for(int n = 0 ; n < stride ; n++)
                            {
                                if(id == (columnBlockIndex*stride+n))
                                {
                                    double & val = getMatrix()[lineBlockIndex*stride+m][columnBlockIndex*stride+n] ;
                                    this->externalForces[lineBlockIndex*stride+m] -= multipliers[p].getValue()*val ;
                                    this->naturalBoundaryConditionForces[lineBlockIndex*stride+m] -= multipliers[p].getValue()*val ;
                                    val = 0 ;
                                }
                            }
                        }
                        else
                        {
                            for(int n = 0 ; n < stride ; n++)
                            {
                                if((columnBlockIndex*stride+n) == id)
                                {
                                    getMatrix()[lineBlockIndex*stride+m][columnBlockIndex*stride+n] = 1 ;
//									count++ ;
                                }
                                else
                                {
                                    getMatrix()[lineBlockIndex*stride+m][columnBlockIndex*stride+n] = 0 ;
                                }
                            }
                        }
                    }
                }
            }


            for(int p = start_multiplier_in_block_index ; p != end_multiplier_in_block_index ; p++)
            {
                if(multipliers[p].type != SET_FORCE_XI
                        &&  multipliers[p].type != SET_FORCE_ETA
                        &&  multipliers[p].type != SET_FORCE_ZETA
                        &&  multipliers[p].type != SET_FORCE_INDEXED_AXIS
                        &&  multipliers[p].type != GENERAL)
                {
                    int id = multipliers[p].getId() ;
                    for(int m = 0 ; m < stride ; m++)
                    {
                        if( id != (lineBlockIndex*stride+m))
                        {
                            for(int n = 0 ; n < stride ; n++)
                            {
                                if(id == (columnBlockIndex*stride+n))
                                {
                                    double & val = getMatrix()[lineBlockIndex*stride+m][columnBlockIndex*stride+n] ;
                                    this->externalForces[lineBlockIndex*stride+m] -= multipliers[p].getValue()*val ;
                                    this->naturalBoundaryConditionForces[lineBlockIndex*stride+m]  -= multipliers[p].getValue()*val ;
                                    val = 0 ;
                                }
                            }
                        }
                        else
                        {
                            externalForces[id] = multipliers[p].getValue() ;
                            for(int n = 0 ; n < stride ; n++)
                            {
                                if((columnBlockIndex*stride+n) == id)
                                {
                                    getMatrix()[lineBlockIndex*stride+m][columnBlockIndex*stride+n] = 1 ;
                                }
                                else
                                {
                                    getMatrix()[lineBlockIndex*stride+m][columnBlockIndex*stride+n] = 0 ;
                                }
                            }
                        }
                    }
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
                || multipliers[i].type == SET_FORCE_ZETA
                || multipliers[i].type == SET_FORCE_INDEXED_AXIS)
        {
            this->externalForces[multipliers[i].getId()] += multipliers[i].getValue() ;
        }
        else if(multipliers[i].type == SET_GLOBAL_FORCE_VECTOR)
        {
            this->externalForces += multipliers[i].coefs ;
        }

    }

    if( dim == SPACE_TWO_DIMENSIONAL &&  element2d[0]->getOrder() >= CONSTANT_TIME_LINEAR)
    {

        size_t totaldofs = getMatrix().row_size.size()*getMatrix().stride ;
        rowstart = totaldofs/2 ;
        colstart = totaldofs/2 ;
        if( element2d[0]->getOrder() == CONSTANT_TIME_QUADRATIC || element2d[0]->getOrder() == LINEAR_TIME_QUADRATIC || element2d[0]->getOrder() == QUADRATIC_TIME_QUADRATIC )
        {
            rowstart = 2*totaldofs/3 ;
            colstart = 2*totaldofs/3 ;
        }
        if(totaldofs < 128)
        {
            colstart = 0 ;
            rowstart = 0 ;
        }
     }
     if( dim == SPACE_THREE_DIMENSIONAL &&  element3d[0]->getOrder() >= CONSTANT_TIME_LINEAR)
     {

        size_t totaldofs = getMatrix().row_size.size()*getMatrix().stride ;
        rowstart = totaldofs/2 ;
        colstart = totaldofs/2 ;
        if( element3d[0]->getOrder() == CONSTANT_TIME_QUADRATIC || element3d[0]->getOrder() == LINEAR_TIME_QUADRATIC || element3d[0]->getOrder() == QUADRATIC_TIME_QUADRATIC )
        {
            rowstart = 2*totaldofs/3 ;
            colstart = 2*totaldofs/3 ;
        }
        if(totaldofs < 128)
        {
            colstart = 0 ;
            rowstart = 0 ;
        }
    }


    prevDisplacements.resize(rowstart) ;
    prevDisplacements = 0. ;

    for(size_t p = 0 ; p < multipliers.size() ; p++)
    {
        size_t id = multipliers[p].getId() ;
        if(id < rowstart)
            prevDisplacements[id] = multipliers[p].value ;
    }

    multipliers.clear() ;


    element2d.clear() ;
    element3d.clear() ;

// 	std::cerr << " ...done." << std::endl ;

}

void Assembly::initialiseElementaryMatrices()
{
    timeval time0, time1 ;
    gettimeofday(&time0, nullptr);

    if(dim == SPACE_TWO_DIMENSIONAL)
    {
// 			#pragma omp parallel for
        for(size_t i = 0 ; i < element2d.size() ; i++)
        {
            std::cerr << "\rGenerating elementary matrices... triangle " << i+1 << "/" << element2d.size() << std::flush ;
            if(element2d[i]->getBehaviour())
            {
                element2d[i]->getElementaryMatrix() ;
                if(element2d[i]->getBehaviour()->isViscous())
                    element2d[i]->getViscousElementaryMatrix() ;
            }
        }
    }
    else if(dim == SPACE_THREE_DIMENSIONAL)
    {
// 			#pragma omp parallel for
        for(size_t i = 0 ; i < element3d.size() ; i++)
        {
            if(i%1000 == 0)
                std::cerr << "\rGenerating elementary matrices... tetrahedron " << i << "/" << element3d.size() << std::flush ;
            if(element3d[i]->getBehaviour())
            {
                element3d[i]->getElementaryMatrix() ;
                if(element3d[i]->getBehaviour()->isViscous())
                    element3d[i]->getViscousElementaryMatrix() ;
            }
        }
    }

    gettimeofday(&time1, nullptr);
    double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
    std::cerr << " ...done. Time to generate (s) " << delta/1e6 << std::endl ;
}

void Assembly::initialiseElementaryMatrices(TetrahedralElement * father)
{
    timeval time0, time1 ;
    gettimeofday(&time0, nullptr);

// 	#pragma omp parallel for
    for(size_t i = 0 ; i < element3d.size() ; i++)
    {
        if(i%1000 == 0)
            std::cerr << "\rGenerating elementary matrices... tetrahedron " << i << "/" << element3d.size() << std::flush ;
        if(element3d[i]->getBehaviour())
        {
            element3d[i]->getElementaryMatrix() ;

            if(element3d[i]->getBehaviour()->isViscous())
                element3d[i]->getViscousElementaryMatrix() ;
        }
    }

    gettimeofday(&time1, nullptr);
    double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
    std::cerr << " ...done. Time to generate (s) " << delta/1e6 << std::endl ;
}

void Assembly::initialiseElementaryMatrices(TriElement * father)
{
    timeval time0, time1 ;
    gettimeofday(&time0, nullptr);
//	std::cerr << "Generating elementary matrices..." << std::flush ;
// 	#pragma omp parallel for
    for(size_t i = 0 ; i < element2d.size() ; i++)
    {

        std::cerr << "\rGenerating elementary matrices... triangle " << i+1 << "/" << element2d.size() << std::flush ;
        if(element2d[i]->getBehaviour())
        {
            element2d[i]->getElementaryMatrix() ;
            if(element2d[i]->getBehaviour()->isViscous())
            {
                element2d[i]->getViscousElementaryMatrix() ;
            }
        }

    }

    gettimeofday(&time1, nullptr);
    double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
    std::cerr << " ...done. Time to generate (s) " << delta/1e6 << std::endl ;
}

void Assembly::checkZeroLines()
{
    std::cerr << "removing 0-only lines..." << std::flush ;
    bool zeros = true ;
    for(size_t i = 0 ; i < ndof ; i++)
    {
        zeros = true ;
        size_t j = 0 ;
        while(zeros && (j < ndof))
        {
            zeros = (getMatrix()[i][j] < POINT_TOLERANCE_2D) ;
            j++ ;
        }
        if(zeros && (j==ndof))
        {
            getMatrix()[i][i] = 1 ;
            externalForces[i] = 0. ;
        }
    }
    std::cerr << "done. " << std::endl ;
}

bool Assembly::make_final()
{
    bool symmetric = true ;

//	size_t ndof = 2 ;

    if(dim == SPACE_TWO_DIMENSIONAL)
    {

        if( element2d.empty() && coordinateIndexedMatrix != nullptr)
        {
            std::cerr << "no elements in mesh (2D) !" << std::endl ;
            return false ;
        }

//		ndof = element2d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        size_t max ;
        if(coordinateIndexedMatrix == nullptr)
        {

            std::set<std::pair<unsigned int, unsigned int> > * map  = new std::set<std::pair<unsigned int, unsigned int> >();
            size_t instants = element2d[0]->timePlanes() ;
            size_t dofsperplane = element2d[0]->getBoundingPoints().size() / instants ;

            for(size_t i = 0 ; i < element2d.size() ; i++)
            {
                if(i%1000 == 0)
                    std::cerr << "\r computing sparsness pattern... triangle " << i+1 << "/" << element2d.size() << std::flush ;
                std::vector<size_t> ids = element2d[i]->getDofIds() ;
                size_t additionalDofPerPlane = ids.size()/instants - dofsperplane ;

                for(size_t j = 0 ; j < dofsperplane * instants ; j++)
                {
                    map->insert(std::make_pair(ids[j], ids[j])) ;
                    if( j >=  dofsperplane * instants - dofsperplane )
                    {
                        for( size_t k = 0 ; k < j ; k++)
                            map->insert( std::make_pair(ids[j], ids[k])) ;
                        for( size_t k = j+1 ; k < ids.size() ; k++)
                            map->insert( std::make_pair(ids[j], ids[k])) ;
                    }
                }

                for(size_t j = dofsperplane * instants ; j < ids.size() ; j++)
                {
                    map->insert(std::make_pair(ids[j], ids[j])) ;

                    if( j >= ids.size() - additionalDofPerPlane )
                    {
                        for( size_t k = 0 ; k < j ; k++)
                            map->insert( std::make_pair(ids[j], ids[k])) ;
                        for( size_t k = j+1 ; k < ids.size() ; k++)
                            map->insert( std::make_pair(ids[j], ids[k])) ;
                    }

                }

// 				for(size_t j = 0 ; j< ids.size() ;j++)
// 				{
//
// 					map->insert(std::make_pair(ids[j], ids[j])) ;
//
// 					for(size_t k = j+1 ; k< ids.size() ;k++)
// 					{
// 						if(k >= dofsperplane*(instants-1))
// 						{
// 							if(j >= dofsperplane*(instants-1))
// 								map->insert(std::make_pair(ids[j], ids[k])) ;
// 							map->insert(std::make_pair(ids[k], ids[j])) ;
// 						}
// 					}
// 				}
            }
            max = map->rbegin()->first +1;
            size_t realDofs = max ;
            for(size_t i = 0 ; i < multipliers.size() ; i++)
            {
                if(multipliers[i].type == GENERAL)
                {
                    std::valarray<unsigned int> ids = multipliers[i].getDofIds() ;
                    for(size_t j = 0 ; j< ids.size() ; j++)
                    {
                        multipliers[i].setId( max ) ;
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
            map->insert(std::make_pair(ndofmax-1,ndofmax-1)) ;

            std::cerr << " ...done" << std::endl ;
            size_t total_num_dof = map->rbegin()->first+1 ;

            std::valarray<unsigned int> row_length((unsigned int)0, total_num_dof) ;
            std::valarray<unsigned int> column_index((unsigned int)0, map->size()) ;

            size_t current = 0 ;
            for(auto i = map->begin() ; i != map->end() ; ++i)
            {
                row_length[i->first]++ ;
                column_index[current] = i->second ;
                current++ ;
            }

            delete map ;
            coordinateIndexedMatrix = new CoordinateIndexedSparseMatrix(row_length, column_index, ndof) ;
            if(max < ndofmax)
                max = ndofmax ;
            if(displacements.size() != max)
            {
                displacements.resize(max) ;
            }
            displacements = 0 ;
        }
        else
        {
            max = coordinateIndexedMatrix->accumulated_row_size.size() ;
        }

        coordinateIndexedMatrix->array = 0 ;
        double dmax = 0 ;
        double vmax = 0 ;
        
        for(size_t i = 0 ; i < element2d.size() ; i++)
        {
            if(!element2d[i]->getBehaviour())
                continue ;
            if(i%1000 == 0)
                std::cerr << "\r computing stiffness matrix... triangle " << i+1 << "/" << element2d.size() << std::flush ;
            std::vector<size_t> ids = element2d[i]->getDofIds() ;
            std::valarray<std::valarray<Matrix > > mother = element2d[i]->getElementaryMatrix();
            Matrix test(ids.size()*ndof, ids.size()*ndof) ;
            for(size_t j = 0 ; j < ids.size() ; j++)
            {
                for(size_t n = 0 ; n < ndof ; n++)
                {
                    for(size_t m = 0 ; m < ndof ; m++)
                    {
                        getMatrix()[ids[j]*ndof+n][ids[j]*ndof+m] += scales[i] * mother[j][j][n][m] ;
                    }
                }
                for(size_t k = j+1 ; k < ids.size() ; k++)
                {
                    for(size_t n = 0 ; n < ndof ; n++)
                    {
                        for(size_t m = 0 ; m < ndof ; m++)
                        {
                            getMatrix()[ids[j]*ndof+n][ids[k]*ndof+m] += scales[i] * mother[j][k][n][m] ;
                            getMatrix()[ids[k]*ndof+n][ids[j]*ndof+m] += scales[i] * mother[k][j][n][m] ;
                            test[j*ndof+n][k*ndof+m] = scales[i] * mother[j][k][n][m] ;
                            test[k*ndof+n][j*ndof+m] = scales[i] * mother[k][j][n][m] ;
                        }
                    }
                }
            }

            if(element2d[i]->getBehaviour()->isViscous())
            {
                std::valarray<std::valarray<Matrix > > vmother = element2d[i]->getViscousElementaryMatrix();
                for(size_t j = 0 ; j < ids.size() ; j++)
                {
                    for(size_t n = 0 ; n < ndof ; n++)
                    {
                        for(size_t m = 0 ; m < ndof ; m++)
                        {
                            getMatrix()[ids[j]*ndof+n][ids[j]*ndof+m] += scales[i] * vmother[j][j][n][m] ;
                        }
                    }
                    for(size_t k = j+1 ; k < ids.size() ; k++)
                    {
                        for(size_t n = 0 ; n < ndof ; n++)
                        {
                            for(size_t m = 0 ; m < ndof ; m++)
                            {
                                getMatrix()[ids[j]*ndof+n][ids[k]*ndof+m] += scales[i] * vmother[j][k][n][m] ;
                                getMatrix()[ids[k]*ndof+n][ids[j]*ndof+m] += scales[i] * vmother[k][j][n][m] ;
                                test[j*ndof+n][k*ndof+m] = scales[i] * vmother[j][k][n][m] ;
                                test[k*ndof+n][j*ndof+m] = scales[i] * vmother[k][j][n][m] ;
                            }
                        }
                    }
                }
            }

            dmax = std::abs(test.array()).max() ;
            if(dmax > POINT_TOLERANCE_2D)
            {
                for(size_t j = 0 ; j < test.numRows() ; j++)
                {
                    for(size_t k = j+1 ; k < test.numCols() ; k++)
                    {
                        vmax = std::abs(test[j][k]-test[k][j]) ;
                        symmetric = symmetric && (vmax/dmax < 1e-12) ;

                    }
                }
            }
//  			element2d[i]->clearElementaryMatrix() ;
        }

        std::cerr << " ...done" << std::endl ;
        getMatrix().stride =  element2d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;

        setBoundaryConditions() ;
        checkZeroLines() ;

    }
    if(dim == SPACE_THREE_DIMENSIONAL)
    {

        if(element3d.empty() && coordinateIndexedMatrix != nullptr)
        {
            std::cerr << "no elements in mesh (3D) !" << std::endl ;
            return false ;
        }

        size_t max ;
//		ndof = element3d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;

        if( coordinateIndexedMatrix == nullptr)
        {
            std::set<std::pair<unsigned int, unsigned int> > * map  = new std::set<std::pair<unsigned int, unsigned int> >();

            for(size_t i = 0 ; i < element3d.size() ; i++)
            {
                if(i%1000 == 0)
                    std::cerr << "\r computing sparsness pattern... tetrahedron " << i+1 << "/" << element3d.size() << std::flush ;
                std::vector<size_t> ids = element3d[i]->getDofIds() ;

                for(size_t j = 0 ; j< ids.size() ; j++)
                {
                    map->insert(std::make_pair(ids[j], ids[j])) ;
                }

                for(size_t j = 0 ; j< ids.size() ; j++)
                {
                    for(size_t k = j+1 ; k< ids.size() ; k++)
                    {
                        map->insert(std::make_pair(ids[j], ids[k])) ;
                        map->insert(std::make_pair(ids[k], ids[j])) ;
                    }
                }
            }

            max = map->rbegin()->first +1;
            size_t realDofs = max ;
            for(size_t i = 0 ; i < multipliers.size() ; i++)
            {
                if(multipliers[i].type == GENERAL)
                {
                    std::valarray<unsigned int> ids = multipliers[i].getDofIds() ;
                    for(size_t j = 0 ; j< ids.size() ; j++)
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
            this->coordinateIndexedMatrix = new CoordinateIndexedSparseMatrix(row_length, column_index, ndof) ;
            if(max < ndofmax)
                max = ndofmax ;
            if(this->displacements.size() != max)
            {
                this->displacements.resize(max) ;
                this->displacements = 0 ;
            }

            std::cerr << " ...done" << std::endl ;
        }
        else
        {
            max = this->coordinateIndexedMatrix->accumulated_row_size.size() ;
        }

        coordinateIndexedMatrix->array = 0 ;
        double dmax = 0 ;
        double vmax = 0 ;

        for(size_t i = 0 ; i < element3d.size() ; i++)
        {
            if(i%1000 == 0)
                std::cerr << "\r computing stiffness matrix... tetrahedron " << i+1 << "/" << element3d.size() << std::flush ;

            std::vector<size_t> ids = element3d[i]->getDofIds() ;
            std::valarray<std::valarray<Matrix > > mother = element3d[i]->getElementaryMatrix();
            for(size_t j = 0 ; j < ids.size() ; j++)
            {
                ids[j] *= ndof ;
            }
            Matrix test(ids.size()*ndof, ids.size()*ndof) ;
            for(size_t j = 0 ; j < ids.size() ; j++)
            {

                for(size_t l = 0 ; l < ndof  ; l++)
                {
                    for(size_t m = 0 ; m < ndof  ; m++)
                    {
                        getMatrix()[ids[j]+l][ids[j]+m] += scales[i] * mother[j][j][l][m] ;
                    }
                }

                for(size_t k = j+1 ; k < ids.size() ; k++)
                {
                    for(size_t l = 0 ; l < ndof  ; l++)
                    {
                        for(size_t m = 0 ; m < ndof  ; m++)
                        {
                            getMatrix()[ids[j]+l][ids[k]+m] += scales[i] * mother[j][k][l][m] ;
                            getMatrix()[ids[k]+l][ids[j]+m] += scales[i] * mother[k][j][l][m] ;
                            test[j*ndof+l][k*ndof+m] = scales[i] * mother[j][k][l][m] ;
                            test[k*ndof+l][j*ndof+m] = scales[i] * mother[k][j][l][m] ;
                        }
                    }
                }
            }

            if(element3d[i]->getBehaviour()->isViscous())
            {
                std::valarray<std::valarray<Matrix > > vmother = element3d[i]->getViscousElementaryMatrix();
                for(size_t j = 0 ; j < ids.size() ; j++)
                {

                    for(size_t l = 0 ; l < ndof  ; l++)
                    {
                        for(size_t m = 0 ; m < ndof  ; m++)
                        {
                            getMatrix()[ids[j]+l][ids[j]+m] += scales[i] * vmother[j][j][l][m] ;
                        }
                    }

                    for(size_t k = j+1 ; k < ids.size() ; k++)
                    {
                        for(size_t l = 0 ; l < ndof  ; l++)
                        {
                            for(size_t m = 0 ; m < ndof  ; m++)
                            {
                                getMatrix()[ids[j]+l][ids[k]+m] += scales[i] * vmother[j][k][l][m] ;
                                getMatrix()[ids[k]+l][ids[j]+m] += scales[i] * vmother[k][j][l][m] ;
                                test[j*ndof+l][k*ndof+m] = scales[i] * vmother[j][k][l][m] ;
                                test[k*ndof+l][j*ndof+m] = scales[i] * vmother[k][j][l][m] ;
                            }
                        }
                    }
                }
            }

            dmax = std::abs(test.array()).max() ;
            if(dmax > POINT_TOLERANCE_3D)
            {
                for(size_t j = 0 ; j < test.numRows() ; j++)
                {
                    for(size_t k = j+1 ; k < test.numCols() ; k++)
                    {
                        vmax = std::abs(test[j][k]-test[k][j]) ;
                        symmetric = symmetric && (vmax/dmax < 1e-12) ;

                    }
                }
            }
// 			element3d[i]->clearElementaryMatrix() ;
        }

//
        getMatrix().stride =   element3d[0]->getBehaviour()->getNumberOfDegreesOfFreedom() ;
        std::cerr << " ...done" << std::endl ;

        setBoundaryConditions() ;
        checkZeroLines() ;
    }
// 	std::cerr << smallestEigenValue(getMatrix()) << std::endl;
    return symmetric ;
}

bool Assembly::mgprepare()
{
    return make_final() ;
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
    if(this->coordinateIndexedMatrix == nullptr)
        make_final() ;

    for(size_t i = 0 ; i < externalForces.size() ; i++)
    {
        for(size_t j = 0 ; j < externalForces.size() ; j++)
            std::cerr << std::setprecision(16) <<getMatrix()[i][j] << "   " << std::flush ;

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

void Assembly::setDisplacementByDof(size_t dof, double val)
{
    displacements[dof] = val ;
}


void Assembly::setPoint(double ex, size_t id)
{
    Vector c(2) ;
    std::valarray<unsigned int> i(2) ;
    auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id)) ;
    if(!(multipliers.empty() || duplicate == multipliers.end()))
        multipliers.erase(duplicate) ;
    multipliers.push_back(LagrangeMultiplier(i,c, ex, id)) ;
    multipliers.back().type = SET_ALONG_XI ;

}

void Assembly::setPoint(double ex, double ey, size_t id)
{
    setSpaceDimension(SPACE_TWO_DIMENSIONAL) ;
    Vector c(2) ;
    std::valarray<unsigned int> i(2) ;

    auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2)) ;
    if(!(multipliers.empty() || duplicate == multipliers.end()))
        multipliers.erase(duplicate) ;

    multipliers.push_back(LagrangeMultiplier(i,c, ex, id*ndof)) ;
    multipliers.back().type = SET_ALONG_XI ;

    if(ndof > 1)
    {
        duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*2+1)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
            multipliers.erase(duplicate) ;

        multipliers.push_back(LagrangeMultiplier(i,c, ey, id*ndof+1)) ;
        multipliers.back().type = SET_ALONG_ETA ;
    }


    return ;
}

void Assembly::setPoint(double ex, double ey, double ez, size_t id)
{
    setSpaceDimension(SPACE_THREE_DIMENSIONAL) ;
    Vector c(2) ;
    std::valarray<unsigned int> i(2) ;

    auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*3)) ;
    if(!(multipliers.empty() || duplicate == multipliers.end()))
        multipliers.erase(duplicate) ;

    multipliers.push_back(LagrangeMultiplier(i,c, ex, id*ndof)) ;
    multipliers.back().type = SET_ALONG_XI ;

    if(ndof > 1)
    {
        duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*3+1)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
            multipliers.erase(duplicate) ;

        multipliers.push_back(LagrangeMultiplier(i,c, ey, id*ndof+1)) ;
        multipliers.back().type = SET_ALONG_ETA ;
    }

    if(ndof > 2)
    {
        duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*3+2)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
            multipliers.erase(duplicate) ;

        multipliers.push_back(LagrangeMultiplier(i,c, ez, id*ndof+2)) ;
        multipliers.back().type = SET_ALONG_ZETA ;
    }
}

void Assembly::setPeriodicPoint(size_t id0, size_t id1)
{
    Vector c(2) ;
    c[0] = 1 ;
    c[1] = -1 ;
    std::valarray<unsigned int> i(2) ;
    for(size_t n = 0 ; n < ndof ; n++)
    {
        i[0] = id0*ndof+n ;
        i[1] = id1*ndof+n ;
        multipliers.push_back(LagrangeMultiplier(i,c, 0.)) ;
    }
}


void Assembly::setNumberOfDregreesOfFreedom(int dof) {
    ndof = dof ;
}
int Assembly::getNumberOfDegreesOfFreedom() const {
    return ndof ;
}

void Assembly::setSpaceDimension(SpaceDimensionality d) {
    dim = d ;
}
SpaceDimensionality Assembly::getSpaceDimension() const {
    return dim ;
}

void Assembly::clear()
{
    multipliers.clear();
    delete coordinateIndexedMatrix ;
    coordinateIndexedMatrix = nullptr ;
    delete nonLinearPartialMatrix ;
    nonLinearPartialMatrix = nullptr ;
    delete boundaryMatrix ;
    boundaryMatrix = nullptr ;
}

void Assembly::clearElements()
{
    element2d.clear() ;
    element3d.clear() ;
    scales.clear();
}

void Assembly::setForceOn(double val, size_t id)
{
    std::valarray<unsigned int> i(2) ;
    Vector c(2) ;
    multipliers.push_back(LagrangeMultiplier(i,c,val, id)) ;
    multipliers.back().type = SET_FORCE_XI ;
}

void Assembly::setForceOn(Variable v, double val, size_t id)
{
    std::valarray<unsigned int> i(2) ;
    Vector c(2) ;
    switch(v)
    {
    case XI:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof)) ;

        if(duplicate != multipliers.end())
        {
            duplicate->value += val ;
        }
        else
        {
            multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof)) ;
            multipliers.back().type = SET_FORCE_XI ;
        }
        break ;
    }
    case ETA:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+1)) ;

        if(duplicate != multipliers.end())
        {
            duplicate->value += val ;
        }
        else
        {
            multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+1)) ;
            multipliers.back().type = SET_FORCE_ETA ;
        }
        break ;
    }
    case ZETA:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+2)) ;

        if(duplicate != multipliers.end())
        {
            duplicate->value += val ;
        }
        else if(!(multipliers.empty() || duplicate == multipliers.end()))
        {
            multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+2)) ;
            multipliers.back().type = SET_FORCE_ZETA ;
        }
        break ;
    }
    default:
    {
        break ;
    }
    }


    return ;
}

void Amie::Assembly::addForceVector(const Vector & v)
{
    std::valarray<unsigned int> i(2) ;
    multipliers.push_back( LagrangeMultiplier(i, v, 0., -1) ) ;
    multipliers.back().type = SET_GLOBAL_FORCE_VECTOR ;
}

void Assembly::addForceOn(Variable v, double val, size_t id)
{
    std::valarray<unsigned int> i(2) ;
    Vector c(2) ;

    switch(v)
    {
    case XI:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
        {
            if((*duplicate).type != SET_FORCE_XI)
                return ;
            (*duplicate).value += val;
            return ;
        }

        multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof)) ;
        multipliers.back().type = SET_FORCE_XI ;
        break ;
    }
    case ETA:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+1)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
        {
            if((*duplicate).type != SET_FORCE_ETA)
                return ;
            (*duplicate).value += val;
            return ;
        }
        multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+1)) ;
        multipliers.back().type = SET_FORCE_ETA ;
        break ;
    }
    case ZETA:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+2)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
        {
            if((*duplicate).type != SET_FORCE_ZETA)
                return ;
            (*duplicate).value += val;
            return ;
        }
        multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+2)) ;
        multipliers.back().type = SET_FORCE_ZETA ;
        break ;
    }
    default:
    {
        break ;
    }
    }

    return ;
}

void Assembly::addMultiplier(const LagrangeMultiplier & l)
{
    auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(l.getId())) ;
    if(!(multipliers.empty() || duplicate == multipliers.end() || duplicate->getId() == -1))
    {
        multipliers.erase(duplicate) ;
    }
    multipliers.push_back(l) ;
}

void Assembly::setPointAlongIndexedAxis(int axis, double val, size_t id, bool force)
{
    std::valarray<unsigned int> i(2) ;
    Vector c(2) ;
    if(!force)
    {
	    auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+axis)) ;
	    if(!(multipliers.empty() || duplicate == multipliers.end()))
	    {
		multipliers.erase(duplicate) ;
	    }
    }

    multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+axis)) ;
    multipliers.back().type = SET_ALONG_INDEXED_AXIS ;
    return ;
}

void Assembly::addForceOnIndexedAxis(int axis, double val, size_t id)
{
    std::valarray<unsigned int> i(2) ;
    Vector c(2) ;
    auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+axis)) ;
    if(!(multipliers.empty() || duplicate == multipliers.end()))
    {
        if((*duplicate).type != SET_FORCE_INDEXED_AXIS)
            return ;
        (*duplicate).value += val;
        return ;
    }

    multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+axis)) ;
    multipliers.back().type = SET_FORCE_INDEXED_AXIS ;

    return ;
}

void Assembly::setPointAlong(Variable v, double val, size_t id)
{
    std::valarray<unsigned int> i(2) ;
    Vector c(2) ;

    switch(v)
    {
    case XI:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
            multipliers.erase(duplicate) ;

        multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof)) ;
        multipliers.back().type = SET_ALONG_XI ;
        break ;
    }
    case ETA:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+1)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
            multipliers.erase(duplicate) ;

        multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+1)) ;
        multipliers.back().type = SET_ALONG_ETA ;
        break ;
    }
    case ZETA:
    {
        auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id*ndof+2)) ;
        if(!(multipliers.empty() || duplicate == multipliers.end()))
            multipliers.erase(duplicate) ;

        multipliers.push_back(LagrangeMultiplier(i,c,val, id*ndof+2)) ;
        multipliers.back().type = SET_ALONG_ZETA ;
        break ;
    }
    default:
    {
        break ;
    }
    }
    return ;

}

void Assembly::fixPoint(size_t id, Amie::Variable v)
{
    auto duplicate = std::find_if(multipliers.begin(), multipliers.end(), MultiplierHasId(id)) ;
    if(!(multipliers.empty() || duplicate == multipliers.end()))
        multipliers.erase(duplicate) ;
    setPointAlong(v, 0, id) ;
}

bool Assembly::solve(Vector x0, size_t maxit, const bool verbose)
{
    if(this->coordinateIndexedMatrix == nullptr)
        make_final() ;

    GaussSeidel gs(this) ;
    bool ret = gs.solve(x0, nullptr) ;
    displacements.resize(gs.x.size()) ;
    displacements = gs.x ;
    return ret ;
}

bool Assembly::cgsolve(Vector x0, int maxit, bool verbose)
{
    bool ret = true ;
// 	if(this->coordinateIndexedMatrix == nullptr)
// 	double lambda_min = smallestEigenValue(getMatrix()) ;
// 	double lambda_max = largestEigenValue(getMatrix()) ;
// // 		std::cout << "largest eigenvalue = " << lambda_max << std::endl ;
// // 		std::cout << "smallest eigenvalue = " << lambda_min << std::endl ;
// 	std::cout << "condition = " << (lambda_max)/(lambda_min) << std::endl ;
    timeval time0, time1 ;
    gettimeofday(&time0, nullptr);

    if(make_final() )
    {

        ConjugateGradientWithSecant cg(this) ;
        if(rowstart > 0 || colstart > 0)
        {

            cg.rowstart = rowstart;
            cg.colstart = colstart;
        }

//         for(size_t i = 0 ; i < 192 ; i++)
//         {
//             for(size_t j = 0 ; j < 192 ; j++)
//                 std::cout << getMatrix()[i][j] << "   " << std::flush ;
//             std::cout << std::endl ;
//         }
        for(size_t i = 0 ; i < 192 ; i++)
            std::cout << getForces()[i] << "   " << std::flush ;

        ret = cg.solve(x0, nullptr, epsilon, -1, verbose) ;
        gettimeofday(&time1, nullptr);
        double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
        std::cerr << "Time to solve (s) " << delta/1e6 << std::endl ;
        displacements.resize(cg.x.size()) ;
        displacements = cg.x ;

        if(rowstart > 0 || colstart > 0)
        {
            for(size_t p = 0 ; p < prevDisplacements.size() ; p++)
            {
                displacements[p] = prevDisplacements[p] ;
            }
        }

    }
    else
    {
        std::cerr << "non-symmetrical problem" << std::endl ;
        timeval time0, time1 ;
        gettimeofday(&time0, nullptr);

        BiConjugateGradientStabilized cg(this) ;
        ret = cg.solve(displacements, nullptr,1e-22, -1, true) ;

        gettimeofday(&time1, nullptr);
        double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
        std::cerr << "Time to solve (s) " << delta/1e6 << std::endl ;
        displacements.resize(cg.x.size()) ;
        displacements = cg.x ;
    }
    return ret ;

}

bool Assembly::mgsolve(LinearSolver * mg, Vector x0, Preconditionner *pg, int Maxit)
{

    std::cerr << "symmetrical problem" << std::endl ;
    timeval time0, time1 ;
    gettimeofday(&time0, nullptr);
    bool ret = mg->solve(x0, pg, 1e-8, -1, true) ;
    gettimeofday(&time1, nullptr);
    double delta = time1.tv_sec*1000000 - time0.tv_sec*1000000 + time1.tv_usec - time0.tv_usec ;
    std::cerr << "Time to solve (s) " << delta/1e6 << std::endl ;
    displacements.resize(mg->x.size()) ;
    displacements = mg->x ;

    return ret ;
}

void Assembly::fix()
{
    make_final() ;
}

bool Assembly::cgnpsolve(const Vector b, size_t maxit)
{
    NullPreconditionner np ;
    ConjugateGradient cg(this) ;
    bool ret = cg.solve(b, &np, 1e-15, -1) ;
    displacements.resize(cg.x.size()) ;
    displacements = cg.x ;
    return ret ;
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

size_t Assembly::getMaxDofID() const
{
    if(coordinateIndexedMatrix)
    {
        if(coordinateIndexedMatrix->row_size.size()*coordinateIndexedMatrix->stride > 0)
            return (coordinateIndexedMatrix->row_size.size()*coordinateIndexedMatrix->stride)/ndof ;
    }

    size_t max = 0 ;
    if(has2DElements())
    {
        for(size_t i = 0 ; i < element2d.size() ; i++)
        {
            for(size_t j = 0 ; j < element2d[i]->getBoundingPoints().size() ; j++)
            {
                max = (element2d[i]->getBoundingPoint(j).getId() > max ? element2d[i]->getBoundingPoint(j).getId() : max ) ;
            }
        }
        return max+1 ;
    }
    if(has3DElements())
    {
        for(size_t i = 0 ; i < element3d.size() ; i++)
        {
            for(size_t j = 0 ; j < element3d[i]->getBoundingPoints().size() ; j++)
            {
                max = (element3d[i]->getBoundingPoint(j).getId() > max ? element3d[i]->getBoundingPoint(j).getId() : max ) ;
            }
        }
        return max+1 ;
    }
    return 0 ;
}


ParallelAssembly::ParallelAssembly(const std::vector<Geometry *> & domains) : domains(domains)
{
    for(size_t i = 0 ; i < domains.size() ; i++)
        assembly.push_back(Assembly());
}

ParallelAssembly::~ParallelAssembly() { } ;

size_t ParallelAssembly::getMaxDofID() const
{
    size_t ret = assembly[0].getMaxDofID() ;
    for(size_t i = 1 ; i < domains.size() ; i++)
        ret = std::max(assembly[i].getMaxDofID(), ret) ;
    
    return ret ;
}

void ParallelAssembly::add(ElementarySurface * e, double scale )
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
    {
        if(domains[i]->in(e->getCenter()))
        {
            assembly[i].add(e, scale) ;
            break ;
        }
    }
}

void ParallelAssembly::add(ElementaryVolume * e, double scale )
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
    {
        if(domains[i]->in(e->getCenter()))
        {
            assembly[i].add(e, scale) ;
            break ;
        }
    }
}

void ParallelAssembly::setMaxDof(size_t n)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setMaxDof(n) ;
}

void ParallelAssembly::setRemoveZeroOnlyLines(bool r) 
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setRemoveZeroOnlyLines(r) ;
}

void ParallelAssembly::initialiseElementaryMatrices()
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].initialiseElementaryMatrices() ;
}
void ParallelAssembly::initialiseElementaryMatrices(TriElement * father)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].initialiseElementaryMatrices(father) ;
}
void ParallelAssembly::initialiseElementaryMatrices(TetrahedralElement * father)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].initialiseElementaryMatrices(father) ;
}

bool  ParallelAssembly::has2DElements() const 
{
    bool ret = false ;
    for( size_t i = 0 ; i < domains.size() ; i++ )
    {
        ret = ret || assembly[i].has2DElements() ;
    }
    return ret ;
}

bool  ParallelAssembly::has3DElements() const 
{
    bool ret = false ;
    for( size_t i = 0 ; i < domains.size() ; i++ )
    {
        ret = ret || assembly[i].has3DElements() ;
    }
    return ret ;
}

void ParallelAssembly::print()
{
    for( size_t i = 0 ; i < domains.size() ; i++ )
        assembly[i].print() ;
}

void ParallelAssembly::printDiag() const
{
    for( size_t i = 0 ; i < domains.size() ; i++ )
        assembly[i].printDiag() ;
}

void ParallelAssembly::setBoundaryConditions()
{
    for( size_t i = 0 ; i < domains.size() ; i++ )
        assembly[i].setBoundaryConditions() ;
}

bool ParallelAssembly::nonLinearStep()
{
    for( size_t i = 0 ; i < domains.size() ; i++ )
        assembly[i].nonLinearStep() ;
}

bool ParallelAssembly::solve(Vector x, size_t maxit, const bool verbose)
{
    std::cout << "solve in parallel !" << std::endl ;
    exit(0) ;
}

bool ParallelAssembly::cgsolve(Vector x0, int maxit, bool verbose )
{
    std::cout << "solve in parallel !" << std::endl ;
    exit(0) ;
}

void ParallelAssembly::setEpsilon(double e) 
{ 
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setEpsilon(e) ;
}

double ParallelAssembly::getEpsilon() const 
{
    return assembly[0].getEpsilon() ;
}

bool ParallelAssembly::mgprepare()
{
    bool ret = false ;
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        ret = assembly[i].mgprepare() || ret ;
    
}

bool ParallelAssembly::mgsolve(LinearSolver * mg, Vector x0, Preconditionner * pg, int maxit)
{
    std::cout << "solve in parallel !" << std::endl ;
    exit(0) ;
}

bool ParallelAssembly::cgnpsolve(Vector b, size_t maxit)
{
    std::cout << "solve in parallel !" << std::endl ;
    exit(0) ;
}

CoordinateIndexedSparseMatrix & ParallelAssembly::getMatrix(int i)
{
    return assembly[i].getMatrix() ;
}

const CoordinateIndexedSparseMatrix & ParallelAssembly::getMatrix(int i) const
{
    return assembly[i].getMatrix() ;
}

CoordinateIndexedIncompleteSparseMatrix & ParallelAssembly::getNonLinearMatrix(int i)
{
    return assembly[i].getNonLinearMatrix() ;
}

Vector & ParallelAssembly::getForces(int i)
{
    return assembly[i].getForces() ;
}

Vector & ParallelAssembly::getNaturalBoundaryConditionForces(int i)
{
    return assembly[i].getNaturalBoundaryConditionForces() ;
}

Vector & ParallelAssembly::getNonLinearForces(int i)
{
    return assembly[i].getNonLinearForces() ;
}

Vector & ParallelAssembly::getDisplacements(int i) 
{
    return assembly[i].getDisplacements() ;
}

void ParallelAssembly::fixPoint(size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].fixPoint(id) ;
}

void ParallelAssembly::setPoint(double ex, size_t id)    
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setPoint(ex,id) ;
}

void ParallelAssembly::setPoint(double ex, double ey, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setPoint(ex, ey,id) ;
}

void ParallelAssembly::setPoint(double ex, double ey, double ez, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setPoint(ex, ey, ez, id) ;
}

void ParallelAssembly::setPointAlong(Variable v, double val, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setPointAlong(v, val, id) ;
}

void ParallelAssembly::setPointAlongIndexedAxis(int index, double val, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setPointAlongIndexedAxis(index, val, id) ;
}

void ParallelAssembly::setForceOn(Variable var, double val, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setForceOn(var, val, id) ;    
}

void ParallelAssembly::addForceOn(Variable var, double val, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].addForceOn(var, val, id) ;      
}

void ParallelAssembly::addForceOnIndexedAxis(int index, double val, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].addForceOnIndexedAxis(index, val, id) ;   
}

void ParallelAssembly::setDisplacementByDof(size_t dof, double val)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setDisplacementByDof(dof, val) ;    
}

void ParallelAssembly::setForceOn(double val, size_t id)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setForceOn(val, id) ;  
}

void ParallelAssembly::addMultiplier(const LagrangeMultiplier & l)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].addMultiplier(l) ;  
}

void ParallelAssembly::addForceVector(const Vector & v)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].addForceVector(v) ;     
}

void ParallelAssembly::fixPoint(size_t id, Variable v)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].fixPoint(id, v) ;        
}

void ParallelAssembly::setPeriodicPoint(size_t id0, size_t id1)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setPeriodicPoint(id0, id1) ;      
}
    
double ParallelAssembly::froebeniusNorm()
{
    double n = 0 ;
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        n += assembly[i].froebeniusNorm() ;
}

void ParallelAssembly::setNumberOfDregreesOfFreedom(int dof)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setNumberOfDregreesOfFreedom(dof) ; 
}

int ParallelAssembly::getNumberOfDegreesOfFreedom() const
{
    return  assembly[0].getNumberOfDegreesOfFreedom() ;
}

void ParallelAssembly::setSpaceDimension(SpaceDimensionality d)
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].setSpaceDimension(d) ;  
}

SpaceDimensionality ParallelAssembly::getSpaceDimension() const
{
    return  assembly[0].getSpaceDimension() ;
}

void ParallelAssembly::fix()
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].fix() ; 
}

void ParallelAssembly::clear()
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].clear() ; 
}

void ParallelAssembly::clearElements()
{
    for( size_t i = 0 ;  i < domains.size() ; i++ )
        assembly[i].clearElements() ; 
}



