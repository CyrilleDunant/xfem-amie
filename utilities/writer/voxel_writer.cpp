//
// C++ Interface: voxel writer
//
// Description:
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "voxel_writer.h"
#include "../../physics/stiffness.h"
#include "../../geometry/geometry_3D.h"
#include <iostream>
#include <fstream>
#include <sstream>


namespace Amie
{

VoxelWriter::VoxelWriter(std::string f, int n)  : filename(f), nVoxelX(n), nVoxelY(n), nVoxelZ(n), fullSample(true) { }

VoxelWriter::VoxelWriter(std::string f, int nx, int ny, int nz): filename(f), nVoxelX(nx), nVoxelY(ny), nVoxelZ(nz), fullSample(true) { }


VoxelWriter::VoxelWriter(std::string f, Point bl, Point tr, int n)
{
    filename = f ;
    nVoxelX = n ;
    nVoxelY = n ;
    nVoxelZ = n ;
    fullSample = false ;
    bottom_left = bl ;
    top_right = tr ;
}

VoxelWriter::VoxelWriter(std::string f, Point bl, Point tr, int nx, int ny, int nz)
{
    filename = f ;
    nVoxelX = nx ;
    nVoxelY = ny ;
    nVoxelZ = nz ;
    fullSample = false ;
    bottom_left = bl ;
    top_right = tr ;
}

void VoxelWriter::write()
{
    std::string originalName = filename ;
    for(size_t j = 0 ; j < values.size() ; j++)
    {
        std::ostringstream newname ;
        newname << originalName << "_"<< nVoxelX << "_" << j << ".vox";
        filename = newname.str() ;
        writeHeader() ;
        std::fstream outbin ;
        outbin.open(filename.c_str(), std::ios::out|std::ios::binary|std::ios::app) ;
        // 	for(size_t j = 0 ; j < values.size() ; j++)
        // 	{
        for(int i = 0 ; i < nPoints() ; i++)
        {
            unsigned char val = values[j][i] ;
            outbin.put(val) ;
        }
        // 	}
        outbin.close() ;
    }
    filename = originalName;
}

void VoxelWriter::getField(FeatureTree * F, VWFieldType field)
{
    values.clear() ;
    std::vector<std::valarray<double> > val = getDoubleValues(F, field) ;
    for(int i = 0 ; i < numberOfFields(field) ; i++)
    {
        values.push_back(normalizeArray(val.back(), voids)) ;
        val.pop_back() ;
    }
}


std::vector<std::valarray<double> > VoxelWriter::getDoubleValues(FeatureTree * F, VWFieldType field)
{

    int max = nVoxelX*nVoxelY*nVoxelZ ;
    int count = 0 ;
    voids.resize(max, false);
    std::vector<std::valarray<double> > ret ;
    for(int i = 0 ; i < numberOfFields(field) ; i++)
    {
        Vector reti(-1e9, max) ;
        ret.push_back(reti) ;
    }

    if(fullSample)
    {
        if(F->getFeature(0))
        {
            Hexahedron * box = dynamic_cast<Hexahedron *>(F->getFeature(0)) ;

            Point c = box->getCenter() ;
            double sx = box->getXSize() ;
            double sy = box->getYSize() ;
            double sz = box->getZSize() ;
            Point vec(sx*(nVoxelX-1)/(nVoxelX),sy*(nVoxelY-1)/(nVoxelY),sz*(nVoxelZ-1)/(nVoxelZ)) ;
            vec *= 0.5 ;
            bottom_left = c-vec ;
            top_right = c+vec ;
        }
        else
        {
            if ( F->get3DMesh()->begin().size() == 0)
            {
                std::cout << "no elements in assembly" << std::endl ;
                return ret;
            }

            double minx = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getX() ;
            double maxx = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getX() ;

            double miny = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getY() ;
            double maxy = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getY() ;

            double minz = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getZ() ;
            double maxz = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getZ() ;

            double mint = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getT() ;
            double maxt = F->get3DMesh()->begin()->getBoundingPoint ( 0 ).getT() ;

            for ( auto i = F->get3DMesh()->begin() ; i != F->get3DMesh()->end()  ; i++ )
            {
                for ( size_t j = 0 ;  j < i->getBoundingPoints().size() ; ++j )
                {
                    if ( i->getBoundingPoint ( j ).getX() < minx )
                    {
                        minx = i->getBoundingPoint ( j ).getX() ;
                    }

                    if ( i->getBoundingPoint ( j ).getX() > maxx )
                    {
                        maxx = i->getBoundingPoint ( j ).getX() ;
                    }

                    if ( i->getBoundingPoint ( j ).getY() < miny )
                    {
                        miny = i->getBoundingPoint ( j ).getY() ;
                    }

                    if ( i->getBoundingPoint ( j ).getY() > maxy )
                    {
                        maxy = i->getBoundingPoint ( j ).getY() ;
                    }

                    if ( i->getBoundingPoint ( j ).getZ() < minz )
                    {
                        minz = i->getBoundingPoint ( j ).getZ() ;
                    }

                    if ( i->getBoundingPoint ( j ).getZ() > maxz )
                    {
                        maxz = i->getBoundingPoint ( j ).getZ() ;
                    }

                    if ( i->getBoundingPoint ( j ).getT() < mint )
                    {
                        mint = i->getBoundingPoint ( j ).getT() ;
                    }

                    if ( i->getBoundingPoint ( j ).getT() > maxt )
                    {
                        maxt = i->getBoundingPoint ( j ).getT() ;
                    }

                }


                Point c((maxx+minx)*.5, (maxy+miny)*.5, (maxz+minz)*.5) ;
                double sx = maxx-minx ;
                double sy = maxy-miny ;
                double sz = maxz-minz ;
                Point vec(sx*(nVoxelX-1)/(nVoxelX),sy*(nVoxelY-1)/(nVoxelY),sz*(nVoxelZ-1)/(nVoxelZ)) ;
                vec *= 0.5 ;
                bottom_left = c-vec ;
                top_right = c+vec ;
            }
        }

    }

    std::cerr << "generating values ( "<<filename << " )... " << count << "/" << F->get3DMesh()->begin().size() << std::flush ;
    int nfields = numberOfFields(field) ;
    for(auto t = F->get3DMesh()->begin() ; t != F->get3DMesh()->end() ; t++)
    {
        std::vector<Point> bb = t->getBoundingBox() ;
  
        double minx = bb[7].getX()  ;
        double miny = bb[7].getY()  ;
        double minz = bb[7].getZ()  ;

        double maxx = bb[0].getX() ; ;
        double maxy = bb[0].getY() ;
        double maxz = bb[0].getZ() ;

        for(int i = nVoxelX*(minx-bottom_left.getX())/((top_right.getX())-(bottom_left.getX())) ; i < nVoxelX*(maxx-bottom_left.getX())/((top_right.getX())-(bottom_left.getX())) ; i++)
        {
            if(i >= 0 && i < nVoxelX)
            {
                for(int j = nVoxelY*(miny-bottom_left.getY())/((top_right.getY())-(bottom_left.getY())) ; j < nVoxelY*(maxy-bottom_left.getY())/((top_right.getY())-(bottom_left.getY())) ; j++)
                {
                    if(j >= 0 && j < nVoxelY)
                    {
                        for(int k = nVoxelZ*(minz-bottom_left.getZ())/((top_right.getZ())-(bottom_left.getZ())) ; k < nVoxelZ*(maxz-bottom_left.getZ())/((top_right.getZ())-(bottom_left.getZ())) ; k++)
                        {
                            if(k >= 0 && k < nVoxelZ)
                            {
                                Point p(bottom_left) ;
                                p.getX() += ((top_right.getX())-(bottom_left.getX()))*((double)(i))/(double(nVoxelX-1)) ;
                                p.getY() += ((top_right.getY())-(bottom_left.getY()))*((double)(j))/(double(nVoxelY-1)) ;
                                p.getZ() += ((top_right.getZ())-(bottom_left.getZ()))*((double)(k))/(double(nVoxelZ-1)) ;
                                if(!t->in(p) || !t->getBehaviour())
                                    continue ;
                                if( t->getBehaviour()->type != VOID_BEHAVIOUR)
                                {
                                    std::pair<bool, std::vector<double> > val = getDoubleValue(t,p,field) ;
                                    std::pair<bool, std::vector<double> > valAlternate ;

                                    if(val.first)
                                    {
                                        for(int m = 0 ; m < nfields ; m++)
                                        {
                                            ret[m][k+nVoxelY*j+nVoxelY*nVoxelX*i] = val.second[m] ;
                                        }
                                    }
                                }
                                else if(t->getBehaviour()->type == VOID_BEHAVIOUR)
                                {
                                    for(int m = 0 ; m < nfields ; m++)
                                    {
                                        ret[m][k+nVoxelY*j+nVoxelY*nVoxelX*i] = 0 ;
                                        voids[k+nVoxelY*j+nVoxelY*nVoxelX*i] = true ;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if(t.getPosition() %100 == 0)
            std::cerr << "\rgenerating values ( "<<filename << " )... " << t.getPosition() <<"/" << t.size() << std::flush ;
    }

    std::cerr << "\rgenerating values ( "<<filename << " )... " << F->get3DMesh()->begin().size()<< "/" << F->get3DMesh()->begin().size() << " ...done." << std::endl ;
    return ret ;
}

std::pair<bool,std::vector<double> > VoxelWriter::getDoubleValue(DelaunayTetrahedron * tet, const Point & p, VWFieldType field)
{
    std::vector<double> ret(numberOfFields(field)) ;

    bool found = false ;
    switch(field)
    {
    case VWFT_PRINCIPAL_ANGLE:
    {
        Vector v(0.,1) ;
        tet->getState().getField( PRINCIPAL_STRESS_ANGLE_FIELD, tet->inLocalCoordinates(p), v, true) ;
        ret[0]=v[0] ;
        found = true ;
        break ;
    }

    case VWFT_STIFFNESS:
    {
        Matrix m = tet->getBehaviour()->getTensor(tet->inLocalCoordinates(p)) ;
        if(!m.isNull())
            ret[0] = m[0][0] ;
        else
            ret[0] = 0 ;
        found = true ;
        break ;
    }

    case VWFT_STRAIN:
    {
        Vector tmp(0.,6) ;
        tet->getState().getField( STRAIN_FIELD, tet->inLocalCoordinates(p), tmp, true) ;
        for(int i = 0 ; i < 6 ; i++)
            ret[i] = tmp[5-i] ;
        found = true ;
        break ;
    }
    

    case VWFT_STRAIN_AND_STRESS:
    {
        Vector tmp1(0.,6) ;
        Vector tmp2(0.,6) ;
        tet->getState().getField( STRAIN_FIELD, REAL_STRESS_FIELD, tet->inLocalCoordinates(p), tmp1, tmp2, true) ;
        for(int i = 0 ; i < 6 ; i++)
        {
            ret[i] = tmp1[5-i] ;
            ret[6+i] = tmp2[5-i] ;
        }
        found = true ;
        break ;
    }

    case VWFT_CONCENTRATION:
    {
        Vector v(0.,1) ;
        tet->getState().getField( DISPLACEMENT_FIELD, tet->inLocalCoordinates(p), v, true) ;
        ret[0]=v[0] ;
        found = true ;
        break ;
    }

    case VWFT_STRESS:
    {
        Vector tmp(0.,6) ;
        tet->getState().getField( REAL_STRESS_FIELD, tet->inLocalCoordinates(p), tmp, true) ;
        for(int i = 0 ; i < 6 ; i++)
            ret[i] = tmp[5-i] ;
        found = true ;
        break ;
    }
        case VWFT_PRINCIPAL_STRESS:
    {
        Vector tmp(0.,3) ;
        tet->getState().getField( PRINCIPAL_REAL_STRESS_FIELD, tet->inLocalCoordinates(p), tmp, true) ;
        for(int i = 0 ; i < 3 ; i++)
            ret[i] = tmp[2-i] ;
        found = true ;
        break ;
    }
        case VWFT_PRINCIPAL_STRAIN:
    {
        Vector tmp(0.,3) ;
        tet->getState().getField( PRINCIPAL_STRAIN_FIELD, tet->inLocalCoordinates(p), tmp, true) ;
        for(int i = 0 ; i < 3 ; i++)
            ret[i] = tmp[2-i] ;
        found = true ;
        break ;
    }

    case VWFT_GRADIENT:
    {
        Vector tmp(0.,3) ;
        tet->getState().getField( GRADIENT_FIELD, tet->inLocalCoordinates(p), tmp, true) ;
        for(int i = 0 ; i < 3 ; i++)
            ret[i] = tmp[2-i] ;
        found = true ;
        break ;
    }

    case VWFT_GRADIENT_AND_FLUX:
    {
        Vector tmp1(0.,3) ;
        Vector tmp2(0.,3) ;
        tet->getState().getField( GRADIENT_FIELD, FLUX_FIELD, tet->inLocalCoordinates(p), tmp1, tmp2, true) ;
        for(int i = 0 ; i < 3 ; i++)
        {
            ret[i] = tmp1[2-i] ;
            ret[3+i] = tmp2[2-i] ;
        }
        found = true ;
        break ;
    }

    case VWFT_FLUX:
    {
        Vector tmp(0.,3) ;
        tet->getState().getField( FLUX_FIELD, tet->inLocalCoordinates(p), tmp, true) ;
        for(int i = 0 ; i < 3 ; i++)
            ret[i] = tmp[2-i] ;
        found = true ;
        break ;
    }

    case VWFT_VON_MISES:
    {
        Vector v(0.,1) ;
        tet->getState().getField( VON_MISES_REAL_STRESS_FIELD, tet->inLocalCoordinates(p), v, true) ;
        ret[0]=v[0] ;
        found = true ;
        break ;
    }
    case VWFT_ENRICHEMENT:
    {
        ret[0]=tet->getEnrichmentFunctions().size() ;
        found = true ;
        break ;
    }
    case VWFT_DAMAGE:
    {
        if(tet->getBehaviour()->getDamageModel())
        {
            Vector s = tet->getBehaviour()->getDamageModel()->getState() ;
            double v = std::inner_product(&s[0], &s[s.size()], &s[0], double(0)) ;
            ret[0]= v ;
            found = true ;

        }
        else
        {
            ret[0]= 0 ;
            found = true ;
        }
        break ;
    }
    }
    return std::make_pair(found, ret) ;
}

void VoxelWriter::writeHeader()
{
    std::fstream outstream ;
    outstream.open(filename.c_str(), std::ios::out) ;
    outstream << "VOXELS" << std::endl ;
    outstream << 1 << std::endl ;
    outstream << nVoxelX << std::endl ;
    outstream << nVoxelY << std::endl ;
    outstream << nVoxelZ << std::endl ;
    outstream.close() ;
}

void VoxelWriter::writeMap(std::string filename, FeatureTree * F, Variable axis, double pos, int n, VWFieldType field, int k, int min, int max)
{
    std::valarray<double> vals((n+1)*(n+1)) ;
    voids.resize(false, (n+1)*(n+1));

    Hexahedron * box = dynamic_cast<Hexahedron *>(F->getFeature(0)) ;
    Point c = box->getCenter() ;
    double sx = box->getXSize() ;
    double sy = box->getYSize() ;
    double sz = box->getZSize() ;

    Point origin = c - Point(sx*0.5,sy*0.5,sz*0.5) ;
    Point xlocal ;
    Point ylocal ;

    switch(axis)
    {
    case XI:
        origin.getX() = pos ;
        xlocal = Point(0.,sy/n,0.) ;
        ylocal = Point(0.,0.,sz/n) ;
        break ;

    case ETA:
        origin.getY() = pos ;
        xlocal = Point(0.,0.,sz/n) ;
        ylocal = Point(sx/n,0.,0.) ;
        break ;

    case ZETA:
        origin.getZ() = pos ;
        xlocal = Point(sx/n,0.,0.) ;
        ylocal = Point(0.,sy/n,0.) ;
        break ;
    }

    int count = 0 ;
    VoxelWriter * dummy = new VoxelWriter("", 0) ;
    for(int i = 0 ; i < n+1 ; i++)
    {
        for(int j = 0 ; j < n+1 ; j++)
        {
            Point p(origin) ;
            p += (xlocal*(double) i) ;
            p += (ylocal*(double) j) ;

            std::vector<DelaunayTetrahedron *> tris = F->get3DMesh()->getConflictingElements(&p) ;
            bool done = false ;
            if(!tris.empty())
            {
                for(size_t l = 0 ; l < tris.size() ; l++)
                {
                    if(tris[l]->in(p) && tris[l]->getBehaviour()->type == VOID_BEHAVIOUR)
                    {
                        std::pair<bool, std::vector<double> > val = dummy->getDoubleValue(tris[l],p,field) ;
                        if(val.first)
                        {
                            vals[count] = val.second[k] ;
                            count++ ;
                            done = true ;
                            break ;
                        }
                    }
                    else if(tris[l]->in(p) && tris[l]->getBehaviour()->type == VOID_BEHAVIOUR)
                    {
                        vals[count] = 0 ;
                        voids[count] = 0 ;
                        count++ ;

                        done = true ;
                        break ;
                    }
                }
            }
            if(!done)
            {
                vals[count] = 0 ;
                count++ ;
            }
        }
    }
    delete dummy ;

    std::valarray<unsigned char> val_int = normalizeArray(vals, voids, min,max) ;
    vals.resize(0) ;

    std::fstream outfile ;
    outfile.open(filename.c_str(), std::ios::out) ;
    count = 0 ;
    for(int i = 0 ; i < n+1 ; i++)
    {
        for(int j = 0 ; j < n+1 ; j++)
        {
            outfile << val_int[count] << " " ;
            count++ ;
        }
        outfile << std::endl ;
    }
}

std::valarray<unsigned char> normalizeArray(const std::valarray<double> & val, const std::valarray<bool> & voids, unsigned short int min, unsigned short int max)
{
    Vector sortedArray = val ;
    std::sort(&sortedArray[0], &sortedArray[sortedArray.size()]) ;
    double vmax = sortedArray[std::min(sortedArray.size()*1-1, (size_t)(sortedArray.size()*.99))] ;
    double vmin = sortedArray[sortedArray.size()*.01] ;
    std::valarray<unsigned char> norm(val.size()) ;
    for(size_t i = 0 ; i < val.size() ; i++)
    {
        if(!voids[i])
            norm[i] = (unsigned char) std::min(std::max(round((double) min + (double)(max-min)*((val[i]-vmin)/(vmax-vmin))), (double)min), (double)max) ;
        else
            norm[i] = 0 ;
    }
    return norm ;
}

int numberOfFields(VWFieldType field)
{
    switch(field)
    {
    case VWFT_PRINCIPAL_ANGLE:
        return 1 ;
    case VWFT_STIFFNESS:
        return 1 ;
    case VWFT_STRAIN:
        return 6 ;
    case VWFT_STRAIN_AND_STRESS:
        return 12 ;
    case VWFT_STRESS:
        return 6 ;
    case VWFT_PRINCIPAL_STRESS:
        return 3 ;
    case VWFT_PRINCIPAL_STRAIN:
        return 3 ;
    case VWFT_CONCENTRATION:
        return 1 ;
    case VWFT_GRADIENT:
        return 3 ;
    case VWFT_GRADIENT_AND_FLUX:
        return 6 ;
    case VWFT_FLUX:
        return 3 ;
    case VWFT_VON_MISES:
        return 1 ;
    default:
        return 1 ;
    }
    return 1 ;
}


} ;




























