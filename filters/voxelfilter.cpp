//
// C++ Implementation: voxelfilter
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "voxelfilter.h"
#include <fstream>
#include <string>

#include "../polynomial/vm_base.h"

using namespace Amie ;

VoxelFilter::VoxelFilter()
{
}


VoxelFilter::~VoxelFilter()
{
}

bool VoxelFilter::existsPath(std::vector<std::vector<std::vector<unsigned char> > > & phase,
                             int isource, int jsource, int ksource,
                             int itarget, int jtarget, int ktarget,
                             int istart, int jstart, int kstart,
                             int iend, int jend, int kend
                            ) const
{

    std::vector<ConnectedNode *> local ;

    ConnectedNode * start = nullptr;
    ConnectedNode * end = nullptr;
    for(int i =istart ; i < iend+1 ; i++)
    {
        for(int j = jstart ; j < jend+1 ; j++)
        {
            for(int k = kstart ; k < kend+1 ; k++)
            {
                if(i >-1 && j>-1 && k >-1&&
                        i < (int)phase.size()&&
                        j < (int)phase[i].size() &&
                        k < (int)phase[i][j].size() &&
                        phase[i][j][k])
                {
                    local.push_back(new ConnectedNode(i, j, k)) ;

                    if(i == isource && j == jsource && k== ksource)
                        start = local[local.size()-1] ;

                    if(i == itarget && j == jtarget && k== ktarget)
                        end = local[local.size()-1] ;
                }
            }
        }
    }

    if(!end)
    {
        for (size_t i = 0 ; i < local.size() ; i++)
            delete local[i] ;

        return true ;
    }

    for(size_t i = 0 ; i < local.size() ; i++)
    {
        for(size_t j = i+1 ; j < local.size() ; j++)
        {
            if(local[i]->isNeighbour(local[j]))
            {
                local[i]->neighbour.push_back(local[j]) ;
                local[j]->neighbour.push_back(local[i]) ;
            }
        }
    }

    std::vector<ConnectedNode *> toCheck = start->neighbour ;
    start->visited = true ;

    while(!toCheck.empty())
    {
        std::vector<ConnectedNode *> temp ;

        for(size_t i = 0 ; i < toCheck.size() ; i++)
        {
            if(!toCheck[i]->visited)
            {
                toCheck[i]->visited = true ;
                if(toCheck[i] == end)
                {
                    for(size_t j = 0 ; j < local.size() ; j++)
                        delete local[j] ;

                    return true ;
                }

                for(size_t j = 0 ; j < toCheck[i]->neighbour.size() ; j++)
                {
                    temp.push_back(toCheck[i]->neighbour[j]) ;
                }
            }
        }

        toCheck = temp ;
    }

    for(size_t i = 0 ; i < local.size() ; i++)
        delete local[i] ;

    return false ;

}



void VoxelFilter::read(const char * filename, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * mesh )
{
    points.clear() ;
    elems.clear() ;
    if(behaviourMap.empty())
    {
        std::cerr << "no behaviours !" << std::endl ;
        return ;
    }

    std::fstream file(filename) ;
    int r ;
    int c ;
    int s ;


    std::string dummy ;

    if(file.is_open())
    {
        file >> dummy  ;
        file >> dummy  ;
        file >> r  ;
        file >> c  ;
        file >> s  ;
    }
    else
    {
        std::cerr << "could not open file !" << std::endl ;
        return ;
    }
    long position = file.tellp() ;
    file.close() ;
    file.open(filename, std::ios::binary|std::ios::in) ;
    file.seekp(position)  ;
    std::cout << "volume is " << r << " x " << c << " x " << s << std::endl ;
//     r = 3 ;
//     c = 3 ;
//     s = 3 ;
    int index = 0 ;
    for( int i = 0 ; i < r+1 ; i++)
    {
        for( int j = 0 ; j < c+1 ; j++)
        {
            for( int k = 0 ; k < s+1 ; k++)
            {
                points.push_back(new Point(100.*((double)i/r), 100.*(double)j/c  ,100.*(double)k/s)) ;
                points.back()->setId(index++) ;
            }
        }
    }
    std::cerr << "generated " << points.size()  << " points" << std::endl ;

    TetrahedralElement * father = new TetrahedralElement(LINEAR) ;

// 	father->compileAndPrecalculate() ;

    std::vector<std::vector<std::vector<unsigned char> > > phase ;
    for( int i = 0 ; i < r ; i++)
    {
        phase.push_back(std::vector<std::vector<unsigned char> >(0)) ;
        for( int j = 0 ; j < c ; j++)
        {
            phase[i].push_back(std::vector<unsigned char>(0)) ;
            for( int k = 0 ; k < s ; k++)
            {
                if(!file.eof())
                {
                    unsigned char behaviourKey ;
                    file >> behaviourKey ;
                    phase[i][j].push_back(behaviourKey) ;
                }
            }
        }
    }

    std::cout << "generated phases" << std::endl ;
    int eindex = 0 ;
    for( int i = 0 ; i < r ; i++)
    {
        for( int j = 0 ; j < c ; j++)
        {
            for( int k = 0 ; k < s ; k++)
            {

                std::vector<Point *> corner ;
                corner.push_back(points[i    *(s+1)*(c+1)    +(j+1)*(s+1) +k  ]) ;
                corner.push_back(points[(i+1)*(s+1)*(c+1)    +(j+1)*(s+1) +k  ]) ;
                corner.push_back(points[i    *(s+1)*(c+1)    +(j+1)*(s+1) +k+1]) ;
                corner.push_back(points[(i+1)*(s+1)*(c+1)    +(j+1)*(s+1) +k+1]) ;
                corner.push_back(points[i    *(s+1)*(c+1)    +j*(s+1)     +k  ]) ;
                corner.push_back(points[(i+1)*(s+1)*(c+1)    +j*(s+1)     +k  ]) ;
                corner.push_back(points[ i   *(s+1)*(c+1)    +j*(s+1)     +k+1] ) ;
                corner.push_back(points[(i+1)*(s+1)*(c+1)    +j*(s+1)     +k+1]) ;


                DelaunayTetrahedron * tet = new DelaunayTetrahedron(mesh, nullptr, corner[0], corner[1], corner[2], corner[4], nullptr) ;
                if(!mesh)
                    tet->index = eindex++ ;
                tet->refresh(father) ;
                elems.push_back(tet) ;
                elems.back()->setBehaviour(mesh,behaviourMap[phase[i][j][k]]->getCopy()) ;

                tet = new DelaunayTetrahedron( mesh,nullptr, corner[3], corner[5], corner[7], corner[6], nullptr) ;
                if(!mesh)
                    tet->index = eindex++ ;
                tet->refresh(father) ;
                elems.push_back(tet) ;
                elems.back()->setBehaviour(mesh,behaviourMap[phase[i][j][k]]->getCopy()) ;

                tet = new DelaunayTetrahedron( mesh,nullptr, corner[2], corner[3], corner[4], corner[6], nullptr) ;
                if(!mesh)
                    tet->index = eindex++ ;
                tet->refresh(father) ;
                elems.push_back(tet) ;
                elems.back()->setBehaviour(mesh,behaviourMap[phase[i][j][k]]->getCopy()) ;

                tet = new DelaunayTetrahedron( mesh,nullptr, corner[1], corner[2], corner[3], corner[4], nullptr) ;
                if(!mesh)
                    tet->index = eindex++ ;
                tet->refresh(father) ;
                elems.push_back(tet) ;
                elems.back()->setBehaviour(mesh,behaviourMap[phase[i][j][k]]->getCopy()) ;

                tet = new DelaunayTetrahedron( mesh,nullptr, corner[1], corner[3], corner[4], corner[5], nullptr) ;
                if(!mesh)
                    tet->index = eindex++ ;
                tet->refresh(father) ;
                elems.push_back(tet) ;
                elems.back()->setBehaviour(mesh,behaviourMap[phase[i][j][k]]->getCopy()) ;

                tet = new DelaunayTetrahedron( mesh,nullptr, corner[3], corner[4], corner[6], corner[5], nullptr) ;
                if(!mesh)
                    tet->index = eindex++ ;
                tet->refresh(father) ;
                elems.push_back(tet) ;
                elems.back()->setBehaviour(mesh,behaviourMap[phase[i][j][k]]->getCopy()) ;
            }
        }
    }

    for( int i = 0 ; i < r ; i++)
    {
        for( int j = 0 ; j < c ; j++)
        {
            for( int k = 0 ; k < s ; k++)
            {
                int elemIndex = (i*s*c+j*s+k)*6 ;
                for(int delta = 0 ; delta < 6 ; delta++ )
                {
                    for(int l = std::max(i-1, 0) ; l <= std::min(i+1, r-1) ; l++)
                    {
                        for(int m = std::max(j-1, 0) ; m <= std::min(j+1, c-1) ; m++)
                        {
                            for(int n = std::max(k-1, 0) ; n <= std::min(k+1, s-1) ; n++)
                            {
                                int auxElemIndex((l*s*c+m*s+n)*6) ;

                                for(int auxdelta = 0 ; auxdelta < 6 ; auxdelta++ )
                                {
                                    if(elemIndex+delta == auxElemIndex+auxdelta)
                                        continue ;
                                    if(elems[elemIndex+delta]->isNeighbour(elems[auxElemIndex+auxdelta]))
                                        elems[elemIndex+delta]->addNeighbour(elems[auxElemIndex+auxdelta]) ;

                                    if(elems[elemIndex+delta]->numberOfCommonVertices(elems[auxElemIndex+auxdelta]) > 0)
                                        elems[elemIndex+delta]->addNeighbourhood(elems[auxElemIndex+auxdelta]) ;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

std::vector<Point *> & VoxelFilter::getPoints()
{
    return points ;
}

std::vector<DelaunayTetrahedron *> & VoxelFilter::getElements()
{
    return elems ;
}

const std::vector<Point *> & VoxelFilter::getPoints() const
{
    return points ;
}

const std::vector<DelaunayTetrahedron *> & VoxelFilter::getElements() const
{
    return elems ;
}

