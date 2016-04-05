// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "parser.h"
#include "../font.h"
#include "../enumeration_translator.h"
#include "../../polynomial/vm_function_extra.h"
#include <sys/stat.h>

using namespace Amie ;

void InclusionParser::readData()
{
    if (!file.fail())
    {
        circs->clear();
        size_t numberOfCircles ;
        size_t index ;
        file >> numberOfCircles ;

        for (size_t i = 0 ; i < numberOfCircles ; i++)
        {
            double centerX, centerY, radius ;
            file >> index >> radius >> centerX >> centerY ;
#ifdef DEBUG
            std::cout << "index : " << index << ", r : " << radius << ", x : " << centerX << ", y : " << centerY << std::endl ;
#endif
            circs->push_back(Circle(radius, centerX, centerY )) ;
        }
    }
}

void SegmentsParser::readData()
{
    if(!file.fail())
    {
        segments->clear() ;
        size_t numberOfSegments ;
        size_t index ;
        file >> numberOfSegments ;

        for (size_t i = 0 ; i < numberOfSegments ; i++)
        {
            if( file.eof())
                return ;

            size_t numberOfHeads, numberOfPoints ;
            file >> index >> numberOfHeads >> numberOfPoints ;
            std::valarray<Point *> m(numberOfPoints) ;

            for ( size_t j = 0 ; j < numberOfPoints ; j++)
            {
                double x, y ;
                file >> x >> y ;
                m[j] = new Point(x,y) ;

            }
            segments->push_back(SegmentedLine(m)) ;
        }
    }
}

void PointParser::readData()
{
    if(!file.fail())
    {
        pointSets->clear() ;

        size_t numberOfPoints, numberOfPointSets ;
        size_t index ;
        file >> numberOfPointSets ;

        for (size_t i = 0 ; i < numberOfPointSets ; i++)
        {
            file >> numberOfPoints ;
            PointSet * m = new PointSet(numberOfPoints) ;
            for (size_t j = 0 ;  j < numberOfPoints ; j++)
            {
                if( file.eof())
                    return ;

                file >> index >> m->getPoint(j)->getX() >> m->getPoint(j)->getY() ;
#ifdef DEBUG
                std::cout << "index : " << index << ", x : "<< m->getPoint(j)->getX() << ", y : " << m->getPoint(j)->getY() << std::endl ;
#endif
            }
            pointSets->push_back(m) ;
        }
    }
}



void PointParser3D::readData()
{
    if(!file.fail())
    {
        pointSets->clear() ;

        size_t numberOfPoints;
        size_t index ;
        size_t indexPeriodic;
        file >> numberOfPoints ;
        for (size_t i = 0 ;  i < numberOfPoints ; i++)
        {
            if( file.eof())
                return ;

            double x,y,z ;
//  			file >>  index>> x >> y >> z ;
            file >> index >> indexPeriodic >> x >> y >> z;
            Point * p = new Point(x,y,z) ;
            p->getId() = index-1 ;
            pointSets->push_back(p) ;
            this->PeriodicIds.push_back(indexPeriodic-1) ;

#ifdef DEBUG
            std::cout << "index : " << index << " indexPeriodic: " << indexPeriodic << ", x : "<< p->getX() << ", y : " << p->getY() << ", z : " << p->getZ() << std::endl ;
#endif
        }
    }

}

void PeriodParser::readData()
{
    if(!file.fail())
    {
        periodSets->clear() ;

        size_t numberOfindex;

        file >> numberOfindex ;
        std::vector<Point *>  p;

        for (size_t i = 0 ;  i < numberOfindex ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> index_id(1) ;
            size_t index ;
            file >> index >> index_id[0] ;

            periodSets->push_back(std::pair<size_t,std::valarray<int> >(index,index_id )) ;
#ifdef DEBUG
            std::cout << "index : " << index << "id "<< index_id[0]<< std::endl ;
#endif
        }
    }

}
void TetrahedronParser::readData()
{
    if(!file.fail())
    {
        tets->clear() ;
        size_t numberOfTets ;
        file >> numberOfTets ;
        for (size_t i = 0 ;  i < numberOfTets ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(4) ;
            size_t mat_index ;
            size_t index ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3] >>  mat_index;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            vertex_id[3]-- ;

            tets->push_back(std::pair<std::valarray<int>, size_t>(vertex_id, mat_index)) ;
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  << ", point 3 :" << vertex_id[3]  << ", material : "<< mat_index<< std::endl ;
#endif
        }
    }
}
void HexahedronParser::readData()
{
    if(!file.fail())
    {
        hex->clear() ;
        size_t numberOfHex ;
        file >> numberOfHex ;
        for (size_t i = 0 ;  i < numberOfHex ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(8) ;
            size_t mat_index ;
            size_t index ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3] >> vertex_id[4] >> vertex_id[5] >> vertex_id[6]  >> vertex_id[7] >> mat_index;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            vertex_id[3]-- ;
            vertex_id[4]-- ;
            vertex_id[5]-- ;
            vertex_id[6]-- ;
            vertex_id[7]-- ;

            hex->push_back(std::pair<std::valarray<int>, size_t>(vertex_id, mat_index)) ;
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  << ", point 3 :" << vertex_id[3]  << ", point 4 :" << vertex_id[4]  << ", point 5 :" << vertex_id[5]  << ", point 6 :" << vertex_id[6]  << ", point 7 :" << vertex_id[7]  << ", material : "<< mat_index<< std::endl ;
#endif
        }
    }
}
void BoundaryParser::readData()
{
    if(!file.fail())
    {
        tri0->clear() ;
        tri1->clear() ;
        tri2->clear() ;
        tri3->clear() ;
        tri4->clear() ;
        tri5->clear() ;
        size_t numberOfTris ;
        file >> numberOfTris ;
        for (size_t i = 0 ;  i < numberOfTris ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(3) ;
            size_t face_index ;
            size_t index ;
            size_t to_discard ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >>  face_index >> to_discard ;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            if(face_index == 1)
            {
                tri0->push_back(vertex_id) ;
            }
            else if(face_index == 2)
            {
                tri1->push_back(vertex_id) ;
            }
            else if(face_index == 3)
            {
                tri2->push_back(vertex_id) ;
            }
            else if(face_index == 4)
            {
                tri3->push_back(vertex_id) ;
            }
            else if(face_index == 5)
            {
                tri4->push_back(vertex_id) ;
            }
            else if(face_index == 6)
            {
                tri5->push_back(vertex_id) ;
            }
            else
            {
                assert(false) ;
            }
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  <<  ", face : "<< face_index<< std::endl ;
#endif
        }
    }
}

std::vector<std::valarray<int> > * BoundaryParser::getData(size_t face_index)
{
    if(face_index == 0)
    {
        return tri0 ;
    }
    else if(face_index == 1)
    {
        return tri1 ;
    }
    else if(face_index == 2)
    {
        return tri2 ;
    }
    else if(face_index == 3)
    {
        return tri3 ;
    }
    else if(face_index == 4)
    {
        return tri4 ;
    }
    else if(face_index == 5)
    {
        return tri5 ;
    }
    else
    {
        assert(false) ;
        return nullptr ; //shut up the compiler ;
    }
}

void HexahedronBoundaryParser::readData()
{
    if(!file.fail())
    {
        carre0->clear() ;
        carre1->clear() ;
        carre2->clear() ;
        carre3->clear() ;
        carre4->clear() ;
        carre5->clear() ;
        size_t numberOfCarres ;
        file >> numberOfCarres ;
        for (size_t i = 0 ;  i < numberOfCarres ; i++)
        {
            if( file.eof())
                return ;

            std::valarray<int> vertex_id(4) ;
            size_t face_index ;
            size_t index ;
            size_t to_discard ;
            file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3]>>  face_index >> to_discard ;
            vertex_id[0]-- ;
            vertex_id[1]-- ;
            vertex_id[2]-- ;
            vertex_id[3]-- ;
            if(face_index == 1)
            {
                carre0->push_back(vertex_id) ;
            }
            else if(face_index == 2)
            {
                carre1->push_back(vertex_id) ;
            }
            else if(face_index == 3)
            {
                carre2->push_back(vertex_id) ;
            }
            else if(face_index == 4)
            {
                carre3->push_back(vertex_id) ;
            }
            else if(face_index == 5)
            {
                carre4->push_back(vertex_id) ;
            }
            else if(face_index == 6)
            {
                carre5->push_back(vertex_id) ;
            }
            else
            {
                assert(false) ;
            }
#ifdef DEBUG
            std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2] << ", point 3 : " << vertex_id[3] <<  ", face : "<< face_index<< std::endl ;
#endif
        }
    }
}

std::vector<std::valarray<int> > * HexahedronBoundaryParser::getData(size_t face_index)
{
    if(face_index == 0)
    {
        return carre0 ;
    }
    else if(face_index == 1)
    {
        return carre1 ;
    }
    else if(face_index == 2)
    {
        return carre2 ;
    }
    else if(face_index == 3)
    {
        return carre3 ;
    }
    else if(face_index == 4)
    {
        return carre4 ;
    }
    else if(face_index == 5)
    {
        return carre5 ;
    }
    else
    {
        assert(false) ;
        return nullptr ; //shut up the compiler ;
    }
}


