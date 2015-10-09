// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2013
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_3D.h"
#include "geometry_2D.h"
#include "../utilities/itoa.h"
#include <fstream>
#include <iomanip>

namespace Amie
{



RegularOctahedron::RegularOctahedron(double s, double x, double y, double z) : ConvexGeometry(6), length(s)
{
    gType = REGULAR_OCTAHEDRON ;
    getCenter().set(x, y, z) ;
    double h = s/std::sqrt(2.) ;
    //the four points of the basis are
    boundingPoints[0] = new Point(x-h, y/*-s*.5*/, z) ;
    boundingPoints[1] = new Point(x/*-s*.5*/, y+h, z) ;
    boundingPoints[2] = new Point(x+h, y/*+s*.5*/, z) ;
    boundingPoints[3] = new Point(x/*+s*.5*/, y-h, z) ;
    //the two points are
    boundingPoints[4] = new Point(x, y, z+h) ;
    boundingPoints[5] = new Point(x, y, z-h) ;

}

RegularOctahedron::RegularOctahedron(double s, Point c) : ConvexGeometry(6), length(s)
{
    gType = REGULAR_OCTAHEDRON ;
    getCenter().set(c) ;
    double h = s/std::sqrt(2.) ;
    //the four points of the basis are
    boundingPoints[0] = new Point(c.getX()-h, c.getY()/*-s*.5*/, c.getZ()) ;
    boundingPoints[1] = new Point(c.getX()/*-s*.5*/, c.getY()+h, c.getZ()) ;
    boundingPoints[2] = new Point(c.getX()+h, c.getY()/*+s*.5*/, c.getZ()) ;
    boundingPoints[3] = new Point(c.getX()/*+s*.5*/, c.getY()-h, c.getZ()) ;
    //the two points are
    boundingPoints[4] = new Point(c.getX(), c.getY(), c.getZ()+h) ;
    boundingPoints[5] = new Point(c.getX(), c.getY(), c.getZ()-h) ;
}


void RegularOctahedron::sampleBoundingSurface(size_t num_points)
{
    std::vector<Point> samplingPoints = getSamplingBoundingPoints( num_points) ;

    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
    {
        delete boundingPoints[i] ;
    }

    random_shuffle(samplingPoints.begin(), samplingPoints.end()) ;

    boundingPoints.resize(samplingPoints.size()) ;
    for(size_t i = 0 ; i < samplingPoints.size() ; i++)
        boundingPoints[i] = new Point(samplingPoints[i]) ;
}

std::vector<Point> RegularOctahedron::getSamplingBoundingPoints(size_t num_points) const
{
    size_t pointsOnEquator = (size_t)round(pow(num_points, 0.33333333)*6.) ;
    Rectangle(length, length, getCenter().getX(), getCenter().getY()) ;
    std::vector<Point> samplingPoints ;
    std::vector<Point> newPoints = Rectangle(length, length, getCenter().getX(), getCenter().getY()).getSamplingBoundingPoints(pointsOnEquator) ;
    size_t realPointsOnEquator = newPoints.size() ;

    Matrix rot(3,3) ;
    rot[0][0] = cos(M_PI*.25) ;
    rot[0][1] = sin(M_PI*.25) ;
    rot[1][0] = -sin(M_PI*.25) ;
    rot[1][1] = cos(M_PI*.25) ;
    rot[2][2] = 1 ;
    for(size_t i = 0 ; i < newPoints.size() ; i++)
    {
        newPoints[i] -= Point(center.getX(), center.getY()) ;
        newPoints[i] *= rot ;
        newPoints[i] += Point(center.getX(), center.getY()) ;
        newPoints[i].getZ() = center.getZ() ;
    }

    samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;

    int iterations = 1 ;
    while(true)
    {
        double factor = (double)(realPointsOnEquator-2.*iterations)/(double)realPointsOnEquator ;
        if(factor < 0)
            break ;
        iterations++ ;
    }
    double totalIterations = iterations ;

    iterations = 1 ;

    double sq2 = std::sqrt(2.) ;
    while(newPoints.size() >= 4)
    {
        double factor = (totalIterations-iterations)/totalIterations ;
        if(factor <= 0)
            break ;
        newPoints = Rectangle(length*factor, length*factor, getCenter().getX(), getCenter().getY()).getSamplingBoundingPoints(realPointsOnEquator-2.*iterations) ;
        if(newPoints.size() >= 4)
        {
            for(size_t i = 0 ; i < newPoints.size() ; i++)
            {
                newPoints[i] -= Point(center.getX(), center.getY()) ;
                newPoints[i] *= rot ;
                newPoints[i] += Point(center.getX(), center.getY()) ;
                newPoints[i].getZ() = center.getZ() + (1.-factor)*length/sq2 ;
            }
            samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;

            for(size_t i = 0 ; i < newPoints.size() ; i++)
            {
                newPoints[i].getZ() = center.getZ() - (1.-factor)*length/sq2;
            }
            samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;

            iterations++ ;
        }
        else
            break ;
    }
    samplingPoints.push_back(*boundingPoints[4]) ;
    samplingPoints.push_back(*boundingPoints[5]) ;

    return samplingPoints ;
}

void RegularOctahedron::sampleSurface(size_t num_points)
{

    sampleBoundingSurface(num_points) ;

    int realPointsOnEquator = Rectangle(length, length, getCenter().getX(), getCenter().getY()).getSamplingBoundingPoints(round(sqrt(num_points))).size() ;
    int iterations = 1 ;
    while(true)
    {
        double factor = (double)(realPointsOnEquator-1.*iterations)/(double)realPointsOnEquator ;
        if(factor < 0)
            break ;
        iterations++ ;
    }

    size_t numberOfShells = round(iterations) ;
    iterations = 1 ;
    std::vector<Point> samplingPoints ;

    for(size_t i = 0 ; i < numberOfShells ; i++)
    {
        double factor = (double)(numberOfShells-1.*iterations)/(double)numberOfShells ;
        if(factor <= 0)
            break ;
        std::vector<Point> newPoints = RegularOctahedron(length*factor, getCenter()).getSamplingBoundingPoints(num_points*factor) ;


        if(newPoints.size() <= 6)
            break ;
        samplingPoints.insert(samplingPoints.end(), newPoints.begin(), newPoints.end()) ;

        iterations++ ;
    }
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    random_shuffle(samplingPoints.begin(), samplingPoints.end()) ;

    inPoints.resize(samplingPoints.size()) ;
    for(size_t i = 0 ; i < samplingPoints.size() ; i++)
        inPoints[i] = new Point(samplingPoints[i]) ;

}

bool RegularOctahedron::in(const Point &p) const
{

    Point local = p-center ;
    double test = std::abs(local.getX()) + std::abs(local.getY()) + std::abs(local.getZ()) ;
    return (test-length/std::sqrt(2.)) < POINT_TOLERANCE ;

// 	return true ;
    Matrix rot(3,3) ;
    rot[0][0] = cos(-M_PI*.25) ;
    rot[0][1] = sin(-M_PI*.25) ;
    rot[1][0] = -sin(-M_PI*.25) ;
    rot[1][1] = cos(-M_PI*.25) ;
    rot[2][2] = 1 ;

    Point v = p ;//*rot ;
    v -=center ;
    v*=rot ;
    v+=center ;
    if(std::abs(v.getX()-center.getX()) > length)
        return false ;
    if(std::abs(v.getY()-center.getY()) > length)
        return false ;
    Point v_(v.getX(), v.getZ()) ;
    if(!Triangle(Point(getCenter().getX()-.5*length, center.getZ()), Point(getCenter().getX()+.5*length, getCenter().getZ()), Point(getCenter().getX(), getCenter().getZ()+0.70711*length)).in(v_) && !Triangle(Point(getCenter().getX()-.5*length, getCenter().getZ()), Point(getCenter().getX()+.5*length, getCenter().getZ()), Point(getCenter().getX(), getCenter().getZ()-0.70711*length)).in(v_))
        return false ;

    v_.set(v.getY(), v.getZ()) ;
    if(!Triangle(Point(getCenter().getY()-.5*length, getCenter().getZ()), Point(getCenter().getY()+.5*length, getCenter().getZ()), Point(getCenter().getY(), getCenter().getZ()+0.70711*length)).in(v_) && !Triangle(Point(getCenter().getY()-.5*length, getCenter().getZ()), Point(getCenter().getY()+.5*length, getCenter().getZ()), Point(getCenter().getY(), getCenter().getZ()-0.70711*length)).in(v_))
        return false ;

    return true;

}

double RegularOctahedron::area() const
{
    return 2.*sqrt(3.)*length*length ;
}

double RegularOctahedron::volume() const
{
    return sqrt(2.)/3.*length*length*length ;
}

void RegularOctahedron::project(Point * p) const
{
}

void RegularOctahedron::computeCenter()
{
}

double RegularOctahedron::getRadius() const
{
    return .5*sqrt(2)*length ;
}


std::vector<Point> RegularOctahedron::getBoundingBox() const
{
    std::vector<Point> ret ;
    ret.push_back(Point(getCenter().getX()+0.5*length, getCenter().getY()+0.5*length, getCenter().getZ()+0.5*length)) ;
    ret.push_back(Point(getCenter().getX()+0.5*length, getCenter().getY()+0.5*length, getCenter().getZ()-0.5*length)) ;
    ret.push_back(Point(getCenter().getX()+0.5*length, getCenter().getY()-0.5*length, getCenter().getZ()+0.5*length)) ;
    ret.push_back(Point(getCenter().getX()+0.5*length, getCenter().getY()-0.5*length, getCenter().getZ()-0.5*length)) ;
    ret.push_back(Point(getCenter().getX()-0.5*length, getCenter().getY()+0.5*length, getCenter().getZ()+0.5*length)) ;
    ret.push_back(Point(getCenter().getX()-0.5*length, getCenter().getY()+0.5*length, getCenter().getZ()-0.5*length)) ;
    ret.push_back(Point(getCenter().getX()-0.5*length, getCenter().getY()-0.5*length, getCenter().getZ()+0.5*length)) ;
    ret.push_back(Point(getCenter().getX()-0.5*length, getCenter().getY()-0.5*length, getCenter().getZ()-0.5*length)) ;
    return ret ;
}



Tetrahedron::Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3): ConvexGeometry(4)
{
    gType = TETRAHEDRON ;

    boundingPoints.resize(4) ;
    boundingPoints[0] = p0 ;
    boundingPoints[1] = p1 ;
    boundingPoints[2] = p2 ;
    boundingPoints[3] = p3 ;

    if(computeVolume() < 0 )
    {
        boundingPoints[0] = p1;
        boundingPoints[1] = p0;
    }

    computeCircumCenter() ;
    Vector r(4) ;
    r[0] = squareDist3D(*p0, circumCenter) ;
    r[1] = squareDist3D(*p1, circumCenter) ;
    r[2] = squareDist3D(*p2, circumCenter) ;
    r[3] = squareDist3D(*p3, circumCenter) ;
    std::sort(&r[0], &r[4]) ;
    sqradius = r[0] ;
    radius = sqrt(sqradius);
    cachedvolume= computeVolume() ;
    computeCenter() ;
    cachedarea = computeArea() ; 
    
}

Tetrahedron::Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7): ConvexGeometry(8)
{
    gType = TETRAHEDRON ;

    boundingPoints.resize(8) ;
    boundingPoints[0] = p0 ;
    boundingPoints[1] = p1 ;
    boundingPoints[2] = p2 ;
    boundingPoints[3] = p3 ;
    boundingPoints[4] = p4 ;
    boundingPoints[5] = p5 ;
    boundingPoints[6] = p6 ;
    boundingPoints[7] = p7 ;

    if(computeVolume() < 0 )
    {
        boundingPoints[0] = p1;
        boundingPoints[1] = p0;
    }
    computeCircumCenter() ;
    Vector r(4) ;
    r[0] = squareDist3D(*p0, circumCenter) ;
    r[1] = squareDist3D(*p1, circumCenter) ;
    r[2] = squareDist3D(*p2, circumCenter) ;
    r[3] = squareDist3D(*p3, circumCenter) ;
    std::sort(&r[0], &r[4]) ;
    sqradius = r[0] ;
    radius = sqrt(sqradius);
    cachedvolume= computeVolume() ;
    computeCenter() ;
    cachedarea = computeArea() ; 
    
}

Tetrahedron::Tetrahedron(const Point &p0, const Point &p1, const Point &p2, const Point &p3): ConvexGeometry(4)
{
    gType = TETRAHEDRON ;

    boundingPoints.resize(4) ;
    boundingPoints[0] = new Point(p0) ;
    boundingPoints[1] = new Point(p1) ;
    boundingPoints[2] = new Point(p2) ;
    boundingPoints[3] = new Point(p3) ;
    if(computeVolume() < 0 )
    {
        std::swap(boundingPoints[0], boundingPoints[1]) ;
    }
    computeCircumCenter() ;
    Vector r(4) ;
    r[0] = squareDist3D(p0, circumCenter) ;
    r[1] = squareDist3D(p1, circumCenter) ;
    r[2] = squareDist3D(p2, circumCenter) ;
    r[3] = squareDist3D(p3, circumCenter) ;
    std::sort(&r[0], &r[4]) ;
    sqradius = r[0] ;
    radius = sqrt(sqradius);
    computeCenter() ;
    cachedarea = computeArea() ; 
    cachedvolume= computeVolume() ;
}

Tetrahedron::Tetrahedron(): ConvexGeometry(4)
{
    gType = TETRAHEDRON ;

    boundingPoints.resize(4) ;
    boundingPoints[2] = new Point(1,0,0) ;
    boundingPoints[3] = new Point(0,1,0) ;
    boundingPoints[0] = new Point(0,0,1) ;
    boundingPoints[1] = new Point(0,0,0) ;
    if(computeVolume() < 0 )
    {
        std::swap(boundingPoints[0], boundingPoints[1]) ;
    }
    computeCircumCenter() ;
    Vector r(4) ;
    r[0] = squareDist3D(*boundingPoints[0], circumCenter) ;
    r[1] = squareDist3D(*boundingPoints[1], circumCenter) ;
    r[2] = squareDist3D(*boundingPoints[2], circumCenter) ;
    r[3] = squareDist3D(*boundingPoints[3], circumCenter) ;
    std::sort(&r[0], &r[4]) ;
    sqradius = r[0] ;
    radius = sqrt(sqradius);
    computeCenter() ;
    cachedarea = computeArea() ; 
    cachedvolume= computeVolume() ;
}

std::vector<Point> Tetrahedron::getSamplingBoundingPoints(size_t num_points) const
{
    return std::vector<Point>() ;
}

void Tetrahedron::sampleBoundingSurface(size_t num_points)
{
    //! \todo make it do something
    return ;
}

void Tetrahedron::computeCenter()
{
    center = Point(0,0,0) ;
    for(size_t i = 0  ; i < boundingPoints.size() ; i++)
    {
        center += (*boundingPoints[i])/boundingPoints.size() ;
    }
}

double Tetrahedron::getRadius() const
{
    return radius ;
}

std::vector<Point> Tetrahedron::getBoundingBox() const
{

    double min_x_0 = 0, min_y_0 = 0, max_x_0 = 0, max_y_0 = 0, max_z_0 = 0, min_z_0 = 0;
    min_y_0 = getBoundingPoint(0).getY() ;
    max_y_0 = getBoundingPoint(0).getY() ;
    min_x_0 = getBoundingPoint(0).getX() ;
    max_x_0 = getBoundingPoint(0).getX() ;
    min_z_0 = getBoundingPoint(0).getZ() ;
    max_z_0 = getBoundingPoint(0).getZ() ;

    for(size_t k  =  1 ; k <  getBoundingPoints().size() ; k++)
    {
        if(getBoundingPoint(k).getY() < min_y_0)
            min_y_0 = getBoundingPoint(k).getY() ;
        if(getBoundingPoint(k).getY() > max_y_0)
            max_y_0 = getBoundingPoint(k).getY() ;

        if(getBoundingPoint(k).getX() < min_x_0)
            min_x_0 = getBoundingPoint(k).getX() ;
        if(getBoundingPoint(k).getX() > max_x_0)
            max_x_0 = getBoundingPoint(k).getX() ;

        if(getBoundingPoint(k).getZ() < min_z_0)
            min_z_0 = getBoundingPoint(k).getZ() ;
        if(getBoundingPoint(k).getZ() > max_z_0)
            max_z_0 = getBoundingPoint(k).getZ() ;
    }


    std::vector<Point> ret ;
    ret.push_back(Point(max_x_0, max_y_0, max_z_0)) ;
    ret.push_back(Point(max_x_0, max_y_0, min_z_0)) ;
    ret.push_back(Point(max_x_0, min_y_0, max_z_0)) ;
    ret.push_back(Point(max_x_0, min_y_0, min_z_0)) ;
    ret.push_back(Point(min_x_0, max_y_0, max_z_0)) ;
    ret.push_back(Point(min_x_0, max_y_0, min_z_0)) ;
    ret.push_back(Point(min_x_0, min_y_0, max_z_0)) ;
    ret.push_back(Point(min_x_0, min_y_0, min_z_0)) ;
    return ret ;
}

const Amie::Point& Tetrahedron::getCircumCenter() const
{
    return circumCenter ;
}

void Tetrahedron::sampleSurface(size_t num_points)
{
    //! \todo make it do something
    return ;
}

bool Tetrahedron::in(const Point & p) const
{
    if(!inCircumSphere(p))
        return false ;

    Point v_ ;

    if(boundingPoints.size() == 4)
    {
        Matrix S(4,4) ;
        S[0][0] = getBoundingPoint(2).getX();
        S[0][1] = getBoundingPoint(3).getX();
        S[0][2] = getBoundingPoint(0).getX() ;
        S[0][3] = getBoundingPoint(1).getX() ;

        S[1][0] = getBoundingPoint(2).getY() ;
        S[1][1] = getBoundingPoint(3).getY();
        S[1][2] = getBoundingPoint(0).getY() ;
        S[1][3] = getBoundingPoint(1).getY() ;

        S[2][0] = getBoundingPoint(2).getZ() ;
        S[2][1] = getBoundingPoint(3).getZ();
        S[2][2] = getBoundingPoint(0).getZ() ;
        S[2][3] = getBoundingPoint(1).getZ() ;

        S[3][0] = 1 ;
        S[3][1] = 1 ;
        S[3][2] = 1 ;
        S[3][3]= 1;

        Vector v(4) ;
        v[0] = p.getX() ;
        v[1] = p.getY() ;
        v[2] = p.getZ() ;
        v[3] = 1 ;

        Vector coeff = inverse4x4Matrix(S) * v ;

        v_ = Point(coeff[0],coeff[1],coeff[2], p.getT());
    }
    else
    {
        Matrix S(4,4) ;
        S[0][0] = getBoundingPoint(4).getX();
        S[0][1] = getBoundingPoint(6).getX();
        S[0][2] = getBoundingPoint(0).getX() ;
        S[0][3] = getBoundingPoint(2).getX() ;

        S[1][0] = getBoundingPoint(4).getY() ;
        S[1][1] = getBoundingPoint(6).getY();
        S[1][2] = getBoundingPoint(0).getY() ;
        S[1][3] = getBoundingPoint(2).getY() ;

        S[2][0] = getBoundingPoint(4).getZ() ;
        S[2][1] = getBoundingPoint(6).getZ();
        S[2][2] = getBoundingPoint(0).getZ() ;
        S[2][3] = getBoundingPoint(2).getZ() ;

        S[3][0] = 1 ;
        S[3][1] = 1 ;
        S[3][2] = 1 ;
        S[3][3]= 1;

        Vector v(4) ;
        v[0] = p.getX() ;
        v[1] = p.getY() ;
        v[2] = p.getZ() ;
        v[3] = 1 ;

        Vector coeff = inverse4x4Matrix(S) * v ;

        v_ = Point(coeff[0],coeff[1],coeff[2],p.getT() );
    }
    if(v_.getX() < -POINT_TOLERANCE)
        return false ;
    if(v_.getY() < -POINT_TOLERANCE)
        return false ;
    if(v_.getZ() < -POINT_TOLERANCE)
        return false ;
    if(v_.getX()+v_.getY()+v_.getZ() > 1+3.*POINT_TOLERANCE)
        return false ;
    return true ;

}

double Tetrahedron::area() const
{
    return cachedarea ;
}

double Tetrahedron::computeArea() 
{
    if(getBoundingPoints().size() == 4)
    {
        Segment s0(getBoundingPoint(1), getBoundingPoint(0)) ;
        Segment s1(getBoundingPoint(1), getBoundingPoint(2)) ;
        Segment s2(getBoundingPoint(1), getBoundingPoint(3)) ;

        return 0.5*((s0.vector()^s1.vector()).norm()+
                    (s0.vector()^s2.vector()).norm()+
                    (s1.vector()^s2.vector()).norm()+
                    (((s1.vector()-s0.vector())^s2.vector())-s0.vector()).norm());
    }
    else
    {
        Segment s0(getBoundingPoint(2), getBoundingPoint(0)) ;
        Segment s1(getBoundingPoint(2), getBoundingPoint(4)) ;
        Segment s2(getBoundingPoint(2), getBoundingPoint(6)) ;
        return 0.5*((s0.vector()^s1.vector()).norm()+
                    (s0.vector()^s2.vector()).norm()+
                    (s1.vector()^s2.vector()).norm()+
                    (((s1.vector()-s0.vector())^s2.vector())-s0.vector()).norm());
    }
}

double Tetrahedron::volume() const
{
    return cachedvolume ;
}

double Tetrahedron::computeVolume() 
{
    if(getBoundingPoints().size() == 4 || timePlanes() > 1)
    {
        return ((getBoundingPoint(1)- getBoundingPoint(2))^(getBoundingPoint(1)- getBoundingPoint(0)))*(getBoundingPoint(1)- getBoundingPoint(3))/6. ;
    }
    else
    {
        return ((getBoundingPoint(2)- getBoundingPoint(4))^(getBoundingPoint(2)- getBoundingPoint(0)))*(getBoundingPoint(2)- getBoundingPoint(6))/6. ;
    }
}

void Tetrahedron::project(Point * p) const
{
//! \todo make it do something

    return ;
}

void Tetrahedron::computeCircumCenter()
{
    //Pick three points from a face
    Point * A = boundingPoints[0] ;
    Point * B = boundingPoints[1] ;
    Point * C = boundingPoints[2] ;
    //Pick last point from a face
    Point * A_ = boundingPoints[3] ;
    Point BC = (*C-*B) ;

    Point BC_mid = (*B+*C)*.5 ;

    Point BA = (*A-*B) ;
    Point BA_mid = (*B+*A)*.5 ;

    Point BA_ = (*A_-*B) ;
    Point BA__mid = (*B+*A_)*.5 ;

    Matrix planes(3,3) ;
    planes[0][0] = BC.getX() ;
    planes[0][1] = BC.getY() ;
    planes[0][2] = BC.getZ() ;
    planes[1][0] = BA.getX() ;
    planes[1][1] = BA.getY() ;
    planes[1][2] = BA.getZ() ;
    planes[2][0] = BA_.getX() ;
    planes[2][1] = BA_.getY() ;
    planes[2][2] = BA_.getZ() ;

    Vector coord(3) ;
    coord[0] = - BC.getX()*BC_mid.getX()  - BC.getY()*BC_mid.getY()  - BC.getZ()*BC_mid.getZ() ;
    coord[1] = - BA.getX()*BA_mid.getX()  - BA.getY()*BA_mid.getY()  - BA.getZ()*BA_mid.getZ() ;
    coord[2] = - BA_.getX()*BA__mid.getX()  - BA_.getY()*BA__mid.getY()  - BA_.getZ()*BA__mid.getZ() ;
    invert3x3Matrix( planes) ;
    coord= planes*coord ;

    circumCenter = Point(-coord[0], -coord[1], -coord[2]) ;


}

bool Tetrahedron::inCircumSphere(const Point & p) const
{
    return inCircumSphere(&p) ;
}

bool Tetrahedron::inCircumSphere(const Point *p) const
{
    if(p->getX() > circumCenter.getX()+1.01*radius)
        return false ;
    if(p->getX() < circumCenter.getX()-1.01*radius)
        return false ;
    if(p->getY() > circumCenter.getY()+1.01*radius)
        return false ;
    if(p->getY() < circumCenter.getY()-1.01*radius)
        return false ;
    if(p->getZ() > circumCenter.getZ()+1.01*radius)
        return false ;
    if(p->getZ() < circumCenter.getZ()-1.01*radius)
        return false ;

    double d = dist(&circumCenter, p) ;

    return  d-radius < POINT_TOLERANCE ;
}

Hexahedron::Hexahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7)
{
    gType = HEXAHEDRON;

    std::vector<Point *> pts ;
    pts.push_back(p0) ;
    pts.push_back(p1) ;
    pts.push_back(p2) ;
    pts.push_back(p3) ;
    pts.push_back(p4) ;
    pts.push_back(p5) ;
    pts.push_back(p6) ;
    pts.push_back(p7) ;

// 	std::stable_sort(pts.begin(), pts.end(), PointLessThan()) ;

    boundingPoints.resize(8) ;
    for(size_t i = 0 ;  i < 8 ; i++)
    {
        boundingPoints[i] = pts[i] ;
    }

    center = Point(*p0 + *p1 + *p2 + *p3 + *p4+ *p5+ *p6+ *p7)*0.125;
    size_x = std::abs(boundingPoints[7]->getX() - boundingPoints[0]->getX());
    size_y = std::abs(boundingPoints[7]->getY() - boundingPoints[0]->getY());
    size_z = std::abs(boundingPoints[7]->getZ() - boundingPoints[0]->getZ());

}

double Hexahedron::getXSize() const
{
    return size_x ;
}

double Hexahedron::getYSize() const
{
    return size_y ;
}

double Hexahedron::getZSize() const
{
    return size_z ;
}

Hexahedron::Hexahedron(Point p0, Point p1, Point p2, Point p3, Point p4, Point p5, Point p6, Point p7)
{
    gType =HEXAHEDRON  ;
    std::vector<Point> pts ;
    pts.push_back(p0) ;
    pts.push_back(p1) ;
    pts.push_back(p2) ;
    pts.push_back(p3) ;
    pts.push_back(p4) ;
    pts.push_back(p5) ;
    pts.push_back(p6) ;
    pts.push_back(p7) ;

    std::stable_sort(pts.begin(), pts.end()) ;

    boundingPoints.resize(8) ;

    for(size_t i = 0 ;  i < 8 ; i++)
    {
        boundingPoints[i] = new Point(pts[i]) ;
    }


    center =  Point(p0 + p1 + p2 + p3 + p4+ p5+ p6+ p7)*0.125;

    size_x = std::abs(boundingPoints[6]->getX() - boundingPoints[0]->getX());
    size_y = std::abs(boundingPoints[6]->getY() - boundingPoints[0]->getY());
    size_z = std::abs(boundingPoints[6]->getZ() - boundingPoints[0]->getZ());

}

Hexahedron::Hexahedron(double x, double y, double z, double originX, double originY, double originZ) {
    gType =HEXAHEDRON  ;
    boundingPoints.resize(8) ;
    boundingPoints[0] = new Point(originX -0.5*x, originY-0.5*y, originZ-0.5*z) ;
    boundingPoints[1] = new Point(originX -0.5*x, originY-0.5*y, originZ+0.5*z) ;
    boundingPoints[2] = new Point(originX -0.5*x, originY+0.5*y, originZ-0.5*z) ;
    boundingPoints[3] = new Point(originX -0.5*x, originY+0.5*y, originZ+0.5*z) ;
    boundingPoints[4] = new Point(originX +0.5*x, originY-0.5*y, originZ-0.5*z) ;
    boundingPoints[5] = new Point(originX +0.5*x, originY-0.5*y, originZ+0.5*z) ;
    boundingPoints[6] = new Point(originX +0.5*x, originY+0.5*y, originZ-0.5*z) ;
    boundingPoints[7] = new Point(originX +0.5*x, originY+0.5*y, originZ+0.5*z) ;

    center = Point(originX, originY, originZ ) ;

    size_x = x ;
    size_y = y ;
    size_z = z ;
}

Hexahedron::Hexahedron(double x, double y, double z, const Point & c)
{
    gType =HEXAHEDRON  ;
    boundingPoints.resize(8) ;
    boundingPoints[0] = new Point(c.getX() -0.5*x, c.getY()-0.5*y, c.getZ()-0.5*z) ;
    boundingPoints[1] = new Point(c.getX() -0.5*x, c.getY()-0.5*y, c.getZ()+0.5*z) ;
    boundingPoints[2] = new Point(c.getX() -0.5*x, c.getY()+0.5*y, c.getZ()-0.5*z) ;
    boundingPoints[3] = new Point(c.getX() -0.5*x, c.getY()+0.5*y, c.getZ()+0.5*z) ;
    boundingPoints[4] = new Point(c.getX() +0.5*x, c.getY()-0.5*y, c.getZ()-0.5*z) ;
    boundingPoints[5] = new Point(c.getX() +0.5*x, c.getY()-0.5*y, c.getZ()+0.5*z) ;
    boundingPoints[6] = new Point(c.getX() +0.5*x, c.getY()+0.5*y, c.getZ()-0.5*z) ;
    boundingPoints[7] = new Point(c.getX() +0.5*x, c.getY()+0.5*y, c.getZ()+0.5*z) ;

    center = c ;

    size_x = x ;
    size_y = y ;
    size_z = z ;
}

Hexahedron::Hexahedron()
{
    gType = HEXAHEDRON ;

    boundingPoints.resize(8) ;
    boundingPoints[0] = new Point(-1,-1,-1) ;
    boundingPoints[1] = new Point(-1,-1,1) ;
    boundingPoints[2] = new Point(-1,1,-1) ;
    boundingPoints[3] = new Point(-1,1,1) ;
    boundingPoints[4] = new Point(1,-1,-1) ;
    boundingPoints[5] = new Point(1,-1,1) ;
    boundingPoints[6] = new Point(1,1,-1) ;
    boundingPoints[7] = new Point(1,1,1) ;

    center =  Point(*boundingPoints[0] +
                          *boundingPoints[1] +
                          *boundingPoints[2] +
                          *boundingPoints[3] +
                          *boundingPoints[4] +
                          *boundingPoints[5] +
                          *boundingPoints[6] +
                          *boundingPoints[7] )*0.125;

    size_x = 2 ;
    size_y = 2 ;
    size_z = 2 ;
}

std::vector<Point> Hexahedron::getSamplingBoundingPoints(size_t num_points) const
{
    std::vector<Point> points ;

    Point point000(*boundingPoints[0]) ;

    Point point111(*boundingPoints[7]) ;

    size_t numPointsPerDirection = std::max(static_cast<size_t>(ceil(pow((double)num_points, 1./3.))),(size_t)2) ;

    std::vector<double> ds ;

    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        double v = (static_cast<double>(i)/static_cast<double>(numPointsPerDirection-1)) ;

        ds.push_back(v) ;
    }

    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            points.push_back(Point(point000.getX(), point111.getY()*ds[i]+point000.getY()*(1.-ds[i]), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
        }
    }

    //face [000] [010] [100]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            Point candidate(point111.getX()*ds[i]+point000.getX()*(1.-ds[i]), point111.getY()*ds[j]+point000.getY()*(1.-ds[j]),point000.getZ()) ;
            bool in = false ;
            for(size_t k = 0 ; k < points.size() ; k++)
            {
                if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                {
                    in = true ;
                    break ;
                }
            }
            if(!in)
                points.push_back(candidate) ;
        }
    }

    //face [000] [001] [100]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            Point candidate(point111.getX()*ds[i]+point000.getX()*(1.-ds[i]),point000.getY(), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j])) ;
            bool in = false ;
            for(size_t k = 0 ; k < points.size() ; k++)
            {
                if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                {
                    in = true ;
                    break ;
                }
            }
            if(!in)
                points.push_back(candidate) ;
        }
    }

    //face [001] [011] [101]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            Point candidate(point000.getX()*ds[i]+point111.getX()*(1.-ds[i]), point111.getY()*ds[j]+point000.getY()*(1.-ds[j]), point111.getZ()) ;
            bool in = false ;
            for(size_t k = 0 ; k < points.size() ; k++)
            {
                if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                {
                    in = true ;
                    break ;
                }
            }
            if(!in)
                points.push_back(candidate) ;
        }
    }

    //face [100] [110] [101]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            Point candidate(point111.getX(), point111.getY()*ds[i]+point000.getY()*(1.-ds[i]), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j])) ;
            bool in = false ;
            for(size_t k = 0 ; k < points.size() ; k++)
            {
                if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                {
                    in = true ;
                    break ;
                }
            }
            if(!in)
                points.push_back(candidate) ;
        }
    }

    //face [010] [011] [110]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            Point candidate(point111.getX()*(1.-ds[i])+point000.getX()*ds[i], point111.getY(), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j])) ;
            bool in = false ;
            for(size_t k = 0 ; k < points.size() ; k++)
            {
                if(squareDist3D(points[k], candidate)< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                {
                    in = true ;
                    break ;
                }
            }
            if(!in)
                points.push_back(candidate) ;
        }
    }

    return points ;
}

void Hexahedron::sampleBoundingSurface(size_t num_points)
{
    std::vector<Point> points = getSamplingBoundingPoints(4*num_points);
    for (size_t i = 0 ; i < boundingPoints.size() ; i++)
        delete boundingPoints[i] ;

    boundingPoints.resize(points.size()) ;
    for(size_t i = 0 ; i < points.size() ; i++)
    {
        boundingPoints[i] = new Point(points[i]) ;
    }

}

void Hexahedron::computeCenter()
{
    for(size_t i = 0  ; i < boundingPoints.size() ; i++)
    {
        center += (*boundingPoints[i])/boundingPoints.size() ;
    }
}

double Hexahedron::getRadius() const
{
    return dist(getCenter(), getCenter()+Point(size_x, size_y, size_z)*.5) ;
}

void Hexahedron::sampleSurface(size_t num_points)
{
    Point point000_(*boundingPoints[0]) ;

    Point point111_(*boundingPoints[7]) ;

    size_t factor = 4 ;

    sampleBoundingSurface(num_points*factor) ;

    int numPointsPerDirection = std::max(static_cast<size_t>(ceil(pow((double)num_points*factor, 1./3.))),(size_t)2) -2;
    int numPointsPerDirectionOnBoundary = numPointsPerDirection+2 ;
    std::vector<Point> points ;

    while(numPointsPerDirection > 0)
    {
        Point point000(point000_ - center) ;
        point000 *= (double)(numPointsPerDirection-1)/(double)(numPointsPerDirectionOnBoundary-1) ;
        point000 += center ;

        Point point111(point111_ - center) ;
        point111 *= (double)(numPointsPerDirection-1)/(double)(numPointsPerDirectionOnBoundary-1) ;
        point111 += center ;

        std::vector<double> ds ;
        ds.push_back(0) ;
        for(int i = 1 ; i < numPointsPerDirection-1 ; i++)
        {
            double v = (static_cast<double>(i)/static_cast<double>(numPointsPerDirection-1)) ;

            ds.push_back(v) ;
        }
        ds.push_back(1) ;

        for(int i = 0 ; i < numPointsPerDirection ; i++)
        {
            for(int j = 0 ; j < numPointsPerDirection ; j++)
            {
                points.push_back(Point(point000.getX(),
                                       point111.getY()*ds[i]+point000.getY()*(1.-ds[i]),
                                       point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
            }
        }

        //face [000] [010] [100]
        for(int i = 0 ; i < numPointsPerDirection ; i++)
        {
            for(int j = 0 ; j < numPointsPerDirection ; j++)
            {
                points.push_back(Point(point111.getX()*ds[i]+point000.getX()*(1.-ds[i]), point111.getY()*ds[j]+point000.getY()*(1.-ds[j]),point000.getZ())) ;
            }
        }

        //face [000] [001] [100]
        for(int i = 0 ; i < numPointsPerDirection ; i++)
        {
            for(int j = 0 ; j < numPointsPerDirection ; j++)
            {
                points.push_back(Point(point111.getX()*ds[i]+point000.getX()*(1.-ds[i]),point000.getY(), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
            }
        }

        //face [001] [011] [101]
        for(int i = 0 ; i < numPointsPerDirection ; i++)
        {
            for(int j = 0 ; j < numPointsPerDirection ; j++)
            {
                points.push_back(Point(point000.getX()*ds[i]+point111.getX()*(1.-ds[i]), point111.getY()*ds[j]+point000.getY()*(1.-ds[j]), point111.getZ())) ;
            }
        }

        //face [100] [110] [101]
        for(int i = 0 ; i < numPointsPerDirection ; i++)
        {
            for(int j = 0 ; j < numPointsPerDirection ; j++)
            {
                points.push_back(Point(point111.getX(), point111.getY()*ds[i]+point000.getY()*(1.-ds[i]), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
            }
        }

        //face [010] [011] [110]
        for(int i = 0 ; i < numPointsPerDirection ; i++)
        {
            for(int j = 0 ; j < numPointsPerDirection ; j++)
            {
                points.push_back(Point(point111.getX()*(1.-ds[i])+point000.getX()*ds[i], point111.getY(), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
            }
        }

        numPointsPerDirection -= 2 ;

    }
    std::sort(points.begin(), points.end()) ;
    auto e = std::unique(points.begin(), points.end()) ;
    points.erase(e, points.end()) ;

    for (size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    inPoints.resize(points.size()) ;
    for(size_t i = 0 ; i < points.size() ; i++)
    {
        Point randomize(.2*std::min(size_z, std::min(size_x, size_y))/numPointsPerDirectionOnBoundary*(double)rand()/RAND_MAX,
                        .2*std::min(size_z, std::min(size_x, size_y))/numPointsPerDirectionOnBoundary*(double)rand()/RAND_MAX,
                        .2*std::min(size_z, std::min(size_x, size_y))/numPointsPerDirectionOnBoundary*(double)rand()/RAND_MAX
                       ) ;
        inPoints[i] = new Point(points[i]+randomize) ;
    }

}

bool Hexahedron::in(const Point & v) const
{
    return v.getX() >= (center.getX() - size_x*.5 - POINT_TOLERANCE) &&
           v.getX() <= (center.getX() + size_x*.5 + POINT_TOLERANCE) &&
           v.getY() >= (center.getY() - size_y*.5 - POINT_TOLERANCE) &&
           v.getY() <= (center.getY() + size_y*.5 + POINT_TOLERANCE) &&
           v.getZ() >= (center.getZ() - size_z*.5 - POINT_TOLERANCE) &&
           v.getZ() <= (center.getZ() + size_z*.5 + POINT_TOLERANCE) ;
}

double Hexahedron::area() const
{
    return 2.*size_x*size_y + 2.*size_x*size_z + 2.*size_y*size_z ;
}

double Hexahedron::volume() const
{
    return size_x*size_y*size_z ;
}

std::vector<Point> Hexahedron::getBoundingBox() const
{
    std::vector<Point> ret ;
    ret.push_back(Point(center.getX()+0.5*size_x, center.getY()+0.5*size_y, center.getZ()+0.5*size_z)) ;
    ret.push_back(Point(center.getX()+0.5*size_x, center.getY()+0.5*size_y, center.getZ()-0.5*size_z)) ;
    ret.push_back(Point(center.getX()+0.5*size_x, center.getY()-0.5*size_y, center.getZ()+0.5*size_z)) ;
    ret.push_back(Point(center.getX()+0.5*size_x, center.getY()-0.5*size_y, center.getZ()-0.5*size_z)) ;
    ret.push_back(Point(center.getX()-0.5*size_x, center.getY()+0.5*size_y, center.getZ()+0.5*size_z)) ;
    ret.push_back(Point(center.getX()-0.5*size_x, center.getY()+0.5*size_y, center.getZ()-0.5*size_z)) ;
    ret.push_back(Point(center.getX()-0.5*size_x, center.getY()-0.5*size_y, center.getZ()+0.5*size_z)) ;
    ret.push_back(Point(center.getX()-0.5*size_x, center.getY()-0.5*size_y, center.getZ()-0.5*size_z)) ;
    return ret ;
}

void Hexahedron::project(Point * p) const
{
    std::vector<Point> bbox = getBoundingBox() ;

    Plane p0(bbox[0], bbox[1], bbox[2]) ;
    Plane p1(bbox[4], bbox[5], bbox[6]) ;
    Plane p2(bbox[0], bbox[1], bbox[4]) ;
    Plane p3(bbox[2], bbox[3], bbox[6]) ;
    Plane p4(bbox[0], bbox[2], bbox[4]) ;
    Plane p5(bbox[1], bbox[3], bbox[5]) ;
    Line l0(bbox[0], bbox[1]-bbox[0]) ;
    Line l1(bbox[2], bbox[2]-bbox[3]) ;
    Line l2(bbox[4], bbox[4]-bbox[5]) ;
    Line l3(bbox[6], bbox[6]-bbox[7]) ;
    Line l4(bbox[0], bbox[2]-bbox[0]) ;
    Line l5(bbox[1], bbox[1]-bbox[3]) ;
    Line l6(bbox[4], bbox[4]-bbox[6]) ;
    Line l7(bbox[5], bbox[5]-bbox[7]) ;
    Line l8(bbox[0], bbox[0]-bbox[4]) ;
    Line l9(bbox[1], bbox[1]-bbox[5]) ;
    Line l10(bbox[2], bbox[2]-bbox[6]) ;
    Line l11(bbox[3], bbox[3]-bbox[7]) ;

    std::map<double, Point> targets ;
    targets[squareDist3D(p0.projection(*p), *p)] = p0.projection(*p) ;
    targets[squareDist3D(p1.projection(*p), *p)] = p1.projection(*p) ;
    targets[squareDist3D(p2.projection(*p), *p)] = p2.projection(*p) ;
    targets[squareDist3D(p3.projection(*p), *p)] = p3.projection(*p) ;
    targets[squareDist3D(p4.projection(*p), *p)] = p4.projection(*p) ;
    targets[squareDist3D(p5.projection(*p), *p)] = p5.projection(*p) ;
    targets[squareDist3D(l0.projection(*p), *p)] = l0.projection(*p) ;
    targets[squareDist3D(l1.projection(*p), *p)] = l1.projection(*p) ;
    targets[squareDist3D(l2.projection(*p), *p)] = l2.projection(*p) ;
    targets[squareDist3D(l3.projection(*p), *p)] = l3.projection(*p) ;
    targets[squareDist3D(l4.projection(*p), *p)] = l4.projection(*p) ;
    targets[squareDist3D(l5.projection(*p), *p)] = l5.projection(*p) ;
    targets[squareDist3D(l6.projection(*p), *p)] = l6.projection(*p) ;
    targets[squareDist3D(l7.projection(*p), *p)] = l7.projection(*p) ;
    targets[squareDist3D(l8.projection(*p), *p)] = l8.projection(*p) ;
    targets[squareDist3D(l9.projection(*p), *p)] = l9.projection(*p) ;
    targets[squareDist3D(l10.projection(*p), *p)] = l10.projection(*p) ;
    targets[squareDist3D(l11.projection(*p), *p)] = l11.projection(*p) ;
    for(size_t i = 0 ; i < 8 ; i++)
        targets[squareDist3D(*p, bbox[i])] = bbox[i] ;

    for(auto i = targets.begin() ; i!=targets.end() ; ++i)
    {
        if(in(i->second))
        {
            *p = i->second ;
            return ;
        }
    }
    *p = targets.begin()->second ;
}

Sphere::Sphere(double r,double x, double y, double z ) : radius(r), sqradius(r*r)
{
    center = Point(x,y,z) ;
    gType = SPHERE ;
}

Sphere::Sphere(double r,const Point * p0 ) : radius(r) , sqradius(r*r)
{
    center = Point(*p0) ;
    gType = SPHERE ;
}

Sphere::Sphere(double r,const Point p0 ) :radius(r) , sqradius(r*r)
{
    center = Point(p0) ;
    gType = SPHERE ;
}

void Sphere::dumpSampleBoundingPoints(size_t n, size_t iter)
{
    Sphere sph(1.,0.,0.,0.) ;
    std::vector<Point> p = sph.getSamplingPointsOnSphere(n,1., iter) ;
    std::valarray<double> d(p.size()*3) ;
    for(size_t i = 0 ; i < p.size() ; i++)
    {
        d[3*i+0] = p[i].getX() ;
        d[3*i+1] = p[i].getY() ;
        d[3*i+2] = p[i].getZ() ;
    }

    std::string file = "../geometry/sphere_mesh_" ;
    file += itoa(n,10) ;
    std::fstream out ;
    out.open(file.c_str(), std::ios::out|std::ios::binary) ;
//	out << p.size() << std::endl ;
    for(size_t i = 0 ; i < d.size() ; i++)
    {
//		unsigned char * buf = (unsigned char *) & d[i] ;
//		for(int j = 0 ; j < sizeof(double) ; j++)
//			out << buf[j] ;
        out << std::setprecision(16) << d[i] ;
    }
    p[0].print() ;
    out.close() ;
}

std::vector<Point> Sphere::importStandardBoundingPoints(size_t n)
{
    std::vector<Point> ret ;
    std::string file = "../geometry/sphere_mesh_" ;
    file += itoa(n,10) ;
    std::fstream in ;
    in.open(file.c_str(), std::ios::in|std::ios::binary) ;
    double dx, dy, dz ;
    for(size_t i = 0 ; i < n ; i++)
    {
        in >> dx ;
        in >> dy ;
        in >> dz ;
        ret.push_back(Point(dx,dy,dz)) ;
    }
    return ret ;
}


std::vector<Point> Sphere::getBoundingBox() const
{
    double r = getRadius() ;
    std::vector<Point> ret ;
    ret.push_back(Point(center.getX()+0.5*r, center.getY()+0.5*r, center.getZ()+0.5*r)) ;
    ret.push_back(Point(center.getX()+0.5*r, center.getY()+0.5*r, center.getZ()-0.5*r)) ;
    ret.push_back(Point(center.getX()+0.5*r, center.getY()-0.5*r, center.getZ()+0.5*r)) ;
    ret.push_back(Point(center.getX()+0.5*r, center.getY()-0.5*r, center.getZ()-0.5*r)) ;
    ret.push_back(Point(center.getX()-0.5*r, center.getY()+0.5*r, center.getZ()+0.5*r)) ;
    ret.push_back(Point(center.getX()-0.5*r, center.getY()+0.5*r, center.getZ()-0.5*r)) ;
    ret.push_back(Point(center.getX()-0.5*r, center.getY()-0.5*r, center.getZ()+0.5*r)) ;
    ret.push_back(Point(center.getX()-0.5*r, center.getY()-0.5*r, center.getZ()-0.5*r)) ;
    return ret ;
}

std::vector<Point> Sphere::getSamplingPointsOnSphere(size_t num_points, double r, size_t iter, size_t threshold) const
{

    size_t import = 128 ;
    while(import < num_points)
        import *= 2 ;


    double ns = ((double) import-(double)num_points)/(double) (import/2) ;
    ns *= import ;
//	ns = std::max(ns,(double)import/4) ;

//	std::cerr << iter << "-" << ns << "-" << import << "-" << num_points << std::endl ;

// 	if(iter > ns)
// 	{
// 		std::vector<Point> standard = getStandardSamplingBoundingPointsOnSphere(num_points) ;
// 		for(size_t i = 0 ; i < standard.size() ; i++)
// 		{
// 			standard[i] -= getCenter() ;
// 			standard[i] *= r ;
// 			standard[i] += getCenter() ;
// 		}
//
// 		return standard ;
// 	}

    //equation of a sphere =
    //x = r sin(theta) cos(phi)
    //y = r sin(theta) sin(phi)
    //z = r cos(theta)

    std::vector<Point> points ;
    if(r  < POINT_TOLERANCE)
        return points ;

// 		first we sample a cube.

    Point point000(-r + center.getX(),
                   -r + center.getY(),
                   -r + center.getZ()) ;
    Point point111(r + center.getX(),
                   r + center.getY(),
                   r + center.getZ()) ;


    size_t numPointsPerDirection = static_cast<size_t>(pow((double)0.25*num_points, 1./3.)) ;
    //face [000] [010] [001]

    std::vector<double> ds ;

    if(numPointsPerDirection > 1)
    {
        for(size_t i = 0 ; i < numPointsPerDirection ; i++)
        {
            double v = (static_cast<double>(i)/static_cast<double>(numPointsPerDirection-1)) ;
            ds.push_back(v) ;
        }
    }
    else
    {
        ds.push_back(1) ;
    }

    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            points.push_back(Point(point000.getX(), point111.getY()*ds[i]+point000.getY()*(1.-ds[i]), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
        }
    }

    //face [000] [010] [100]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            points.push_back(Point(point111.getX()*ds[i]+point000.getX()*(1.-ds[i]), point111.getY()*ds[j]+point000.getY()*(1.-ds[j]),point000.getZ())) ;
        }
    }

    //face [000] [001] [100]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            points.push_back(Point(point111.getX()*ds[i]+point000.getX()*(1.-ds[i]),point000.getY(), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
        }
    }

    //face [001] [011] [101]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            points.push_back(Point(point000.getX()*ds[i]+point111.getX()*(1.-ds[i]), point111.getY()*ds[j]+point000.getY()*(1.-ds[j]), point111.getZ())) ;
        }
    }

    //face [100] [110] [101]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            points.push_back(Point(point111.getX(), point111.getY()*ds[i]+point000.getY()*(1.-ds[i]), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
        }
    }

    //face [010] [011] [110]
    for(size_t i = 0 ; i < numPointsPerDirection ; i++)
    {
        for(size_t j = 0 ; j < numPointsPerDirection ; j++)
        {
            points.push_back(Point(point111.getX()*(1.-ds[i])+point000.getX()*ds[i], point111.getY(), point111.getZ()*ds[j]+point000.getZ()*(1.-ds[j]))) ;
        }
    }



    for(size_t i = points.size() ; i < num_points ; i++)
    {
        double rx = ((double)rand()/(RAND_MAX))*2.-1. ;
        double ry = ((double)rand()/(RAND_MAX))*2.-1. ;
        double rz = ((double)rand()/(RAND_MAX))*2.-1. ;
        points.push_back(Point(rx, ry, rz)+center) ;
        project(&points.back(),r) ;
    }

    bool haveDuplicates = true ;
    while(haveDuplicates)
    {
        haveDuplicates = false ;
        for(size_t i  = 0 ; i < points.size() ; i++)
        {
            for(size_t j  = i+1 ; j < points.size() ; j++)
            {
                if(squareDist3D(points[i], points[j])< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                {
                    haveDuplicates = true ;
                    points.erase(points.begin()+j) ;
                    break ;
                }
            }

            if(haveDuplicates)
                break ;
        }
    }

    smooth(points, r, iter*10) ;

    return points ;
}

void Sphere::smooth(std::vector<Point> & points,double r, size_t iter) const
{
    if(points.empty())
        return ;
    std::valarray<Point> speeds(points.size()) ;
    Point vec ;
    double error = 2. ;
    double last_error = 1. ;
    int count = 0 ;
    for(size_t i = 0 ; /*(i < iter) &&*/ std::abs(error-last_error)/last_error > 0.01*0.01*points.size()*points.size() && (count == 0); i++)
    {

        last_error = error ;
        error = 0. ;
        for(size_t j = 0 ; j < points.size() ; j++)
        {
            for(size_t k = j+1 ; k < points.size() ; k++)
            {
                if(squareDist3D( points[j], points[k]) > 512.*POINT_TOLERANCE*POINT_TOLERANCE)
                {
                    vec.set(points[j].getX()-points[k].getX(),points[j].getY()-points[k].getY(),points[j].getZ()-points[k].getZ()) ;
                    double n = vec.sqNorm() ;
                    error += n ;
                    vec *= (sqradius/n) ;

                    speeds[j] += vec ;
                    speeds[k] -= vec ;
                }
                else
                {
                    count++ ;
                    Point vec(r,0,0) ;

                    speeds[j] += vec ;
                    speeds[k] -= vec ;
                }
            }
        }


        for(size_t j = 0 ; j < points.size() ; j++)
        {
            points[j] += speeds[j];
            project(&points[j],r) ;
        }
        if(last_error < 0.01*getRadius())
            break ;

        speeds = Point() ;
    }

}

std::vector<Point> Sphere::getStandardSamplingBoundingPointsOnSphere(size_t n) const
{
    size_t import = 128 ;
    while(import < n)
        import *= 2 ;

    if(import > 16384)
        import = 16384 ;


    std::vector<Point> p = Sphere::importStandardBoundingPoints(import) ;
    std::random_shuffle(p.begin(), p.end());
    if(n <= p.size())
        p.erase(p.begin()+n, p.end()) ;


    Point c = getCenter() ;
    for(size_t i = 0 ; i < p.size() ; i++)
        p[i] += c ;

    double ns = ((double) import-n)/(double) (import-import/2) ;
    ns *= import ;

    smooth(p,1., int(ns*2)) ;

    return p ;
}

std::vector<Point> Sphere::getSamplingBoundingPoints(size_t num_points) const
{
    return getSamplingPointsOnSphere(num_points, radius) ;
}

void Sphere::sampleBoundingSurface(size_t num_points)
{
    std::vector<Point> points = getSamplingPointsOnSphere(num_points, radius) ;
    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        delete boundingPoints[i] ;

    boundingPoints.resize(points.size()) ;
    for(size_t i = 0 ; i < points.size() ; i++)
    {
        boundingPoints[i] = new Point(points[i]) ;
    }


}

void Sphere::sampleSurface(size_t num_points)
{

    if(num_points < 4)
        num_points = 4 ;

    sampleBoundingSurface(num_points) ;

    std::vector<Point> points ;

    size_t numPointsOnSurface = boundingPoints.size() ;
    size_t numberOfRadii = static_cast<size_t>(round(sqrt(numPointsOnSurface/6))) ;

    for(size_t i = 0 ; i < numberOfRadii ; i++)
    {

        if(i == 0)
            numPointsOnSurface = 6*numPointsOnSurface/7 ;
        if(i == 1)
            numPointsOnSurface = 5*numPointsOnSurface/6 ;
        if(i == 2)
            numPointsOnSurface = 4*numPointsOnSurface/5 ;
        if(i == 3)
            numPointsOnSurface = 3*numPointsOnSurface/4 ;
        if(i == 4)
            numPointsOnSurface = 2*numPointsOnSurface/3 ;
        if(i == 5)
            numPointsOnSurface = 1*numPointsOnSurface/2 ;

        double r = radius*((double)(numberOfRadii-i-1)/(double)numberOfRadii) ;
        num_points = (size_t)((double)numPointsOnSurface*((r*r)/(radius*radius)));
        if(num_points < 8 )
            break ;

        std::vector<Point> newPoints = getSamplingPointsOnSphere(num_points, r) ;

        for(size_t j = 0 ; j < newPoints.size() ; j++)
        {
// 			double r = ((double)rand()/(double)(RAND_MAX+1)*2. -1.)*.0005*(sqrt(radius)/(double)(numberOfRadii+1)) ;
// 			Point dr = (newPoints[j] - center)*r ;
            points.push_back(newPoints[j]/*+ dr*/) ;
        }
    }
    points.push_back(center) ;
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    inPoints.resize(points.size()) ;

    for(size_t i = 0 ; i < points.size() ; i++)
    {
        inPoints[i] = new Point(points[i]) ;
    }
}

bool Sphere::in(const Point & v) const
{
    if(v.getX() < center.getX()-1.0001*getRadius())
        return false ;
    if(v.getX() > center.getX()+1.0001*getRadius())
        return false ;
    if(v.getY() < center.getY()-1.0001*getRadius())
        return false ;
    if(v.getY() > center.getY()+1.0001*getRadius())
        return false ;
    if(v.getZ() < center.getZ()-1.0001*getRadius())
        return false ;
    if(v.getZ() > center.getZ()+1.0001*getRadius())
        return false ;
    return squareDist3D(v, getCenter()) < getRadius()*getRadius() + 2.*getRadius()*POINT_TOLERANCE + POINT_TOLERANCE*POINT_TOLERANCE;
}

double Sphere::area() const
{
    return 4.*M_PI*getRadius()*getRadius() ;
}

double Sphere::volume() const
{
    return 1.333333333333333333*M_PI*sqradius*radius ;
}

void Sphere::project(Point * p) const
{
    project(p, getRadius()) ;
}

void Sphere::project(Point * p, double r) const
{
    //x = r sin(theta) cos(phi)
    //y = r sin(theta) sin(phi)
    //z = r cos(theta)
    int id = p->getId() ;
    if(squareDist3D(*p, getCenter() ) < POINT_TOLERANCE*POINT_TOLERANCE)
    {
        p->getX() += r ;
        return ;
    }

    Point p_prime = *p-getCenter() ;
    double n = p_prime.norm() ;
    if(std::abs(n - r) < POINT_TOLERANCE)
        return ;
    p_prime *= r/n ;
    *p = getCenter()+p_prime ;
    p->getId() = id ;
    return ;

}

void Sphere::computeCenter() { } 

double Sphere::getRadius() const
{
    return radius ;
}

void Sphere::setRadius(double newr)
{

    double ratio = newr/radius ;

    for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
    {
        getBoundingPoint(i).getX() = (getBoundingPoint(i).getX() - center.getX())*ratio + center.getX() ;
        getBoundingPoint(i).getY() = (getBoundingPoint(i).getY() - center.getY())*ratio + center.getY() ;
    }

    for(size_t i = 0 ; i < getInPoints().size() ; i++)
    {
        getInPoint(i).getX() = (getInPoint(i).getX() - center.getX())*ratio + center.getX() ;
        getInPoint(i).getY() = (getInPoint(i).getY() - center.getY())*ratio + center.getY() ;
    }

    sqradius = newr*newr ;
    radius = newr ;
}

void TriangulatedSurface::addPoint(Point * p)
{
//     const Point * nearest = boundary[0] ;
    for(size_t i = 0 ; i < boundary.size() ; i++)
    {

    }
}

TriangulatedSurface::TriangulatedSurface(const Point * p0, const Point * p1, const Point * p2)
{
    mesh.push_back(TriPoint(p0, p1, p2)) ;
    boundary.push_back(p0) ;
    boundary.push_back(p1) ;
    boundary.push_back(p2) ;

    center = (*p0 + *p1 + *p2)/3 ;
}

TriangulatedSurface::TriangulatedSurface()
{
    mesh.push_back(TriPoint(new Point(0,0,1), new Point(1,0,0), new Point(0,0,0))) ;
    boundary.push_back(&mesh[0].point[0]) ;
    boundary.push_back(&mesh[0].point[1]) ;
    boundary.push_back(&mesh[0].point[2]) ;
    center = Point(1./3., 0, 1./3.) ;
}

double TriangulatedSurface::area() const
{
    double ret = 0 ;
    for(size_t i = 0 ; i < mesh.size() ; i++)
    {
        ret += mesh[i].area() ;
    }

    return ret ;
}

void TriangulatedSurface::project(Point * p) const
{

}

void TriangulatedSurface::computeCenter()
{

}

double TriangulatedSurface::getRadius() const {
    return 0 ;
}

std::vector<Point> TriangulatedSurface::getBoundingBox() const {
    return std::vector<Point>(0) ;
}


void PolygonPrism::computeCenter()
{
    center = axis*.5 + origin ;
}
void LoftedPolygonPrism::computeCenter()
{
    center = interpolatingPointAndTangent(0.5).first ;
}


PolygonPrism::PolygonPrism(const std::valarray<Point *> & points, const Point & vector, const Point & origin) : NonConvexGeometry(0), base(points), axis(vector), origin(origin), rotationMatrix(3,3)
{
    transform(&base,TRANSLATE, -base.getCenter());
    gType = POLYGON_PRISM ;
    computeCenter();
    
    Point renormAxis = vector/vector.norm() ;
    Matrix rotx(3,3) ;
    Matrix rotx_(3,3) ;
    Matrix roty(3,3) ;

    if(std::abs(renormAxis.getY()) > POINT_TOLERANCE)
    {
        double angle = atan2( renormAxis.getY(),renormAxis.getZ()) ;
        double cost = cos(angle) ;
        double sint = sin(angle) ;

        rotx[0][0] = 1 ;
        rotx[1][1] = cost ;
        rotx[1][2] = -sint ;
        rotx[2][1] = sint ;
        rotx[2][2] = cost ;

        Point inter = rotx*renormAxis ;
        angle = atan2( inter.getX(),inter.getZ()) ;

        cost = cos(angle) ;
        sint = sin(angle) ;
        roty[0][0] = cost ;
        roty[0][2] = -sint ;
        roty[1][1] = 1 ;
        roty[2][0] = sint ;
        roty[2][2] = cost ;
        inter = roty*inter ;
        
        angle = atan2( inter.getY(),inter.getZ()) ;
        cost = cos(angle) ;
        sint = sin(angle) ;

        rotx_[0][0] = 1 ;
        rotx_[1][1] = cost ;
        rotx_[1][2] = -sint ;
        rotx_[2][1] = sint ;
        rotx_[2][2] = cost ;

        rotationMatrix = inverse3x3Matrix(rotx*roty*rotx_)  ;
//         std::swap(rotationMatrix[2][0],rotationMatrix[0][0]) ;
//         std::swap(rotationMatrix[2][1],rotationMatrix[0][1]) ;
//         std::swap(rotationMatrix[2][2],rotationMatrix[0][2]) ;
    }
    else
    {
        double angle = atan2( renormAxis.getX(),renormAxis.getZ()) ;
        double cost = cos(angle) ;
        double sint = sin(angle) ;

        rotx[0][0] = cost ;
        rotx[0][2] = -sint ;
        rotx[1][1] = 1 ;
        rotx[2][0] = sint ;
        rotx[2][2] = cost ;

        Point inter = rotx*renormAxis ;
        angle = atan2( inter.getY(),inter.getZ()) ;

        cost = cos(angle) ;
        sint = sin(angle) ;
        roty[0][0] = 1 ;    
        roty[1][1] = cost ; ;
        roty[1][2] = -sint ;
        roty[2][1] = sint ; 
        roty[2][2] = cost ; 
        inter = roty*inter ;
        
        angle = atan2( inter.getX(),inter.getZ()) ;
        cost = cos(angle) ;
        sint = sin(angle) ;

        rotx_[0][0] = cost ;
        rotx_[0][2] = -sint ;
        rotx_[1][1] = 1 ;
        rotx_[2][0] = sint ;
        rotx_[2][2] = cost ;
        
        rotationMatrix = inverse3x3Matrix(rotx*roty*rotx_)  ;
//         std::swap(rotationMatrix[2][0],rotationMatrix[0][0]) ;
//         std::swap(rotationMatrix[2][1],rotationMatrix[0][1]) ;
//         std::swap(rotationMatrix[2][2],rotationMatrix[0][2]) ;
    }
}

LoftedPolygonPrism::LoftedPolygonPrism(const std::valarray< Point* >& points, const std::vector< Point >& ip) :  NonConvexGeometry(0),base(points),interpolationPoints(ip)
{
    vstart = 2.*interpolationPoints[0] - interpolationPoints[1] ;
    vend = 2.*interpolationPoints[interpolationPoints.size()-1]-interpolationPoints[interpolationPoints.size()-2] ;

    transform(&base,TRANSLATE, -base.getCenter());
    gType = LOFTED_POLYGON ;
    computeCenter();    
    centerInXYPlane = true ;
    centerInXZPlane = false ;
    centerInYZPlane = false ; 
    tanAtCentre = interpolatingPointAndTangent(0.5).second ;
    if(std::abs(tanAtCentre*Point(0,1,0)) > std::abs(tanAtCentre*Point(0,0,1)))
    {
        centerInXYPlane = false ;
        centerInXZPlane = true ;
        centerInYZPlane = false ;
        
        if(std::abs(tanAtCentre*Point(1,0,0)) > std::abs(tanAtCentre*Point(0,1,0)))
        {
            centerInXYPlane = false ;
            centerInXZPlane = false ;
            centerInYZPlane = true ;
        }
    }

   
//     for(double t = 0. ; t < 1 ; t += 0.01)
//     {
//         interpolatingPointAndTangent(t, Point(-10, -10)) ;
//         interpolatingPointAndTangent(t, Point(-10,  10)) ;
//         interpolatingPointAndTangent(t, Point(-0 ,  20)) ;
//         interpolatingPointAndTangent(t, Point( 10,  10)) ;
//         interpolatingPointAndTangent(t, Point( 10, -10)) ;
//     }
//     exit(0) ;
}

Point catmullRom(double t_, const std::vector<Point> & controlPoints)
{
        double t0 = 0 ; 
        double t1 = dist(controlPoints[0], controlPoints[1]) ;
        double t2 = t1+dist(controlPoints[1], controlPoints[2]) ;
        double t3 = t2+dist(controlPoints[2], controlPoints[3]) ;
        t_ *= t2-t1 ;
        t_ += t1 ;
        
        
        Point A1 = controlPoints[0]*(t1 - t_)/(t1-t0) + controlPoints[1]*(t_-t0)/(t1-t0);
        Point A2 = controlPoints[1]*(t2 - t_)/(t2-t1) + controlPoints[2]*(t_-t1)/(t2-t1);
        Point A3 = controlPoints[2]*(t3 - t_)/(t3-t2) + controlPoints[3]*(t_-t2)/(t3-t2);
        Point B1 = A1*(t2 - t_)/(t2-t0) +  A2*(t_-t0)/(t2-t0) ;
        Point B2 = A2*(t3 - t_)/(t3-t1) +  A3*(t_-t1)/(t3-t1) ;
        
        return  B1*(t2-t_)/(t2-t1) + B2*(t_-t1)/(t2-t1)  ;
}

Point catmullRomTangent(double t_, const std::vector<Point> & controlPoints)
{
        double t0 = 0 ; 
        double t1 = dist(controlPoints[0], controlPoints[1]) ;
        double t2 = t1+dist(controlPoints[1], controlPoints[2]) ;
        double t3 = t2+dist(controlPoints[2], controlPoints[3]) ;
        t_ *= t2-t1 ;
        t_ += t1 ;
        
        
        Point A1 = controlPoints[0]*(t1 - t_)/(t1-t0) + controlPoints[1]*(t_-t0)/(t1-t0);
        Point A1_ = -controlPoints[0]/(t1-t0) + controlPoints[1]/(t1-t0);
        Point A2 = controlPoints[1]*(t2 - t_)/(t2-t1) + controlPoints[2]*(t_-t1)/(t2-t1);
        Point A2_ = -controlPoints[1]/(t2-t1) + controlPoints[2]/(t2-t1);
        Point A3 = controlPoints[2]*(t3 - t_)/(t3-t2) + controlPoints[3]*(t_-t2)/(t3-t2);
        Point A3_ = -controlPoints[2]/(t3-t2) + controlPoints[3]/(t3-t2);
        Point B1 = A1*(t2 - t_)/(t2-t0) +  A2*(t_-t0)/(t2-t0) ;
        Point B1_ = A1_*(t2 - t_)/(t2-t0) - A1/(t2-t0)+  A2_*(t_-t0)/(t2-t0) + A2/(t2-t0) ;
        Point B2 = A2*(t3 - t_)/(t3-t1) +  A3*(t_-t1)/(t3-t1) ;
        Point B2_ = A2_*(t3 - t_)/(t3-t1) - A2/(t3-t1) +  A3_*(t_-t1)/(t3-t1) +  A3/(t3-t1);
       
        Point tan = B1_*(t2-t_)/(t2-t1) - B1/(t2-t1) + B2_*(t_-t1)/(t2-t1) + B2/(t2-t1)  ; 
        return tan / tan.norm() ;                                   
}

std::pair<Point,Point> LoftedPolygonPrism::interpolatingPointAndTangent(double t) const
{
    
    size_t endIndex = ceil(t*(interpolationPoints.size()-1)) ;
    if(endIndex > interpolationPoints.size()-1)
        endIndex = interpolationPoints.size()-1 ;
    if(endIndex < 1 )
        endIndex = 1 ;
    size_t startindex = endIndex-1 ;
    
    std::vector<Point> controlPoints ;
    if(startindex == 0)
    {
        controlPoints.push_back(vstart);
        
        if(controlPoints.size() == 2)
        {
            controlPoints.push_back(interpolationPoints[0]);
            controlPoints.push_back(interpolationPoints[1]);
            controlPoints.push_back(vend);
        }
        else
        {
            controlPoints.push_back(interpolationPoints[0]);
            controlPoints.push_back(interpolationPoints[1]);
            controlPoints.push_back(interpolationPoints[2]);
        }
        
    }
    else if(startindex == interpolationPoints.size()-2)
    {
        if(controlPoints.size() == 2)
        {
            controlPoints.push_back(vstart);
            controlPoints.push_back(interpolationPoints[0]);
            controlPoints.push_back(interpolationPoints[1]);
            controlPoints.push_back(vend);
        }
        else
        {
            controlPoints.push_back(interpolationPoints[startindex-1]);
            controlPoints.push_back(interpolationPoints[startindex]);
            controlPoints.push_back(interpolationPoints[startindex+1]);
            controlPoints.push_back(vend);
        }
    }
    else
    {
        controlPoints.push_back(interpolationPoints[startindex-1]);
        controlPoints.push_back(interpolationPoints[startindex]);
        controlPoints.push_back(interpolationPoints[startindex+1]);
        controlPoints.push_back(interpolationPoints[startindex+2]);
    }

    double t_ =  t*(interpolationPoints.size()-1.) - startindex ;
    t_ = std::min(std::max(t_, 0.), 1.) ;

    Point ipoint = catmullRom(t_, controlPoints) ;
    Point itan = catmullRomTangent(t_, controlPoints); 

    return std::make_pair(ipoint,itan) ;
}

std::pair<Point,Point> LoftedPolygonPrism::interpolatingPointAndTangent(double t, const Point & offset) const
{
    size_t endIndex = ceil(t*(interpolationPoints.size()-1)) ;
    if(endIndex > interpolationPoints.size()-1)
        endIndex = interpolationPoints.size()-1 ;
    if(endIndex < 1 )
        endIndex = 1 ;
    size_t startindex = endIndex-1 ;
    
    std::vector<Point> controlPoints ;
    if(startindex == 0)
    {
        controlPoints.push_back(vstart);
        
        if(controlPoints.size() == 2)
        {
            controlPoints.push_back(interpolationPoints[0]);
            controlPoints.push_back(interpolationPoints[1]);
            controlPoints.push_back(vend);
        }
        else
        {
            controlPoints.push_back(interpolationPoints[0]);
            controlPoints.push_back(interpolationPoints[1]);
            controlPoints.push_back(interpolationPoints[2]);
        }
        
    }
    else if(startindex == interpolationPoints.size()-2)
    {
        if(controlPoints.size() == 2)
        {
            controlPoints.push_back(vstart);
            controlPoints.push_back(interpolationPoints[0]);
            controlPoints.push_back(interpolationPoints[1]);
            controlPoints.push_back(vend);
        }
        else
        {
            controlPoints.push_back(interpolationPoints[startindex-1]);
            controlPoints.push_back(interpolationPoints[startindex]);
            controlPoints.push_back(interpolationPoints[startindex+1]);
            controlPoints.push_back(vend);
        }
    }
    else
    {
        controlPoints.push_back(interpolationPoints[startindex-1]);
        controlPoints.push_back(interpolationPoints[startindex]);
        controlPoints.push_back(interpolationPoints[startindex+1]);
        controlPoints.push_back(interpolationPoints[startindex+2]);
    }

    double t_ =  t*(interpolationPoints.size()-1.) - startindex ;
    t_ = std::min(std::max(t_, 0.), 1.) ;
     
    Point ipoint = catmullRom(t_, controlPoints) ;
    Point itan = catmullRomTangent(t_, controlPoints);
    Matrix transform = rotateFromVector(Point(0,0,1),itan) ;
    Point rof = offset*transform ;
    ipoint  = ipoint + rof ;
    
    return std::make_pair(ipoint,itan) ;
}

Matrix LoftedPolygonPrism::rotateToVector(const Point & tovector,const Point & fromvector) const
{   
//     inverse3x3Matrix(rotateFromVector(tovector,fromvector)).print() ;
    return inverse3x3Matrix(rotateFromVector(tovector,fromvector)) ;
}

Matrix LoftedPolygonPrism::rotateFromVector(const Point & fromvector, const Point & tovector) const
{
    double nprod = (fromvector.norm()*tovector.norm()) ;

    Point renormAxis = fromvector^tovector ;
    double norm = renormAxis.norm() ;
    if(nprod < POINT_TOLERANCE*POINT_TOLERANCE || norm < POINT_TOLERANCE)
    {
        double angle = acos(tanAtCentre*tovector) + M_PI *.5;
        if(tovector.getY()*tovector.getX() > 0)
            angle = -angle + M_PI;
        Matrix rot(3,3) ;
        rot[0][0] = cos(-angle) ; rot[0][1] = sin(-angle) ; 
        rot[1][0] = -sin(-angle) ; rot[1][1] = cos(-angle) ; rot[2][2] = 1 ;
        return rot ;
    }
    double angle = asin(norm/nprod) ;
    renormAxis /= norm ;

    double q0, q1, q2, q3 ;
    q1 = renormAxis.getX()*sin(angle*.5) ;
    q2 = renormAxis.getY()*sin(angle*.5) ;
    q3 = renormAxis.getZ()*sin(angle*.5) ;
    q0 = cos(angle*.5) ;
    
    double n = sqrt(q1*q1+q0*q0+q2*q2+q3*q3) ;
    q0 /= n ;
    q1 /= n ;
    q2 /= n ;
    q3 /= n ;
    Matrix transform(3,3) ;
    transform[0][0] = q0*q0+q1*q1-q2*q2-q3*q3 ; transform[0][1] = 2.*(q1*q2-q0*q3)        ; transform[0][2] = 2.*(q0*q2+q1*q3) ;
    transform[1][0] = 2.*(q1*q2+q0*q3)        ; transform[1][1] = q0*q0-q1*q1+q2*q2-q3*q3 ; transform[1][2] = 2.*(q3*q2-q1*q0) ;
    transform[2][0] = 2.*(q1*q3-q0*q2)        ; transform[2][1] = 2.*(q0*q1+q2*q3)        ; transform[2][2] = q0*q0-q1*q1-q2*q2+q3*q3 ;

    double nn = tovector.norm()*tanAtCentre.norm() ;
    if(nn < POINT_TOLERANCE*POINT_TOLERANCE || std::abs((tanAtCentre*tovector)/nn - 1.) < POINT_TOLERANCE)
        angle = .5*M_PI ;
    else
    {
        angle = acos((tanAtCentre*tovector)/nn) + .5*M_PI;
        if(tovector.getY()*tovector.getX() > 0)
            angle = -angle + M_PI ;
    }
    Matrix rot(3,3) ;
    rot[0][0] = cos(-angle) ; rot[0][1] = sin(-angle) ; 
    rot[1][0] = -sin(-angle) ; rot[1][1] = cos(-angle) ; rot[2][2] = 1 ;

    return rot*transform ;
 
}

PolygonPrism::~PolygonPrism() { }

LoftedPolygonPrism::~LoftedPolygonPrism() { }

void PolygonPrism::sampleBoundingSurface(size_t num_points)
{
    std::vector<Point> newPoints = getSamplingBoundingPoints(num_points) ;
    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        delete boundingPoints[i] ;

    boundingPoints.resize(newPoints.size());

    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        boundingPoints[i] = new Point(newPoints[i]) ;

}

void LoftedPolygonPrism::sampleBoundingSurface(size_t num_points)
{
    std::vector<Point> newPoints = getSamplingBoundingPoints(num_points) ;
    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        delete boundingPoints[i] ;

    boundingPoints.resize(newPoints.size());

    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        boundingPoints[i] = new Point(newPoints[i]) ;

}

std::vector<Point> PolygonPrism::getSamplingBoundingPoints(size_t num_points) const
{
    std::vector<Point> newPoints ;

    std::valarray<Point *> pts(base.originalPoints.size()) ;
    for(size_t i = 0 ; i < pts.size() ; i++)
        pts[i] = new Point(base.originalPoints[i]) ;

    Polygon basel(pts) ;
    basel.sampleSurface(num_points) ;
    for(size_t i = 0 ; i < basel.getInPoints().size() ; i++)
        newPoints.push_back(basel.getInPoint(i));
    for(size_t i = 0 ; i < basel.getBoundingPoints().size() ; i++)
        newPoints.push_back(basel.getBoundingPoint(i));

    for(size_t j = 0 ; j< newPoints.size() ; j++)
    {
        newPoints[j] = rotationMatrix * newPoints[j];
        newPoints[j] += origin ;
    }

    double numSlices = std::max(round(axis.norm()/basel.getPerimeter()*num_points), 2.) ;
    std::vector<Point> ret ;
    for(size_t j = 0 ; j< newPoints.size() ; j++)
    {
        ret.push_back(newPoints[j]);
    }
    for(size_t j = 0 ; j< newPoints.size() ; j++)
    {
        ret.push_back( newPoints[j] + axis);
    }

    for(double i = 1 ; i< numSlices-1.1 ; i++)
    {
        for(size_t j = 0 ; j< basel.getBoundingPoints().size() ; j++)
        {
            Point p = basel.getBoundingPoint(j) ;
            p = rotationMatrix*p ;
            p += origin ;
            ret.push_back( axis * i/(numSlices-1.) + p);
        }
    }

    return ret ;
}

std::vector<Point> LoftedPolygonPrism::getSamplingBoundingPoints(size_t num_points) const
{
    std::vector<Point> newPoints ;

    std::valarray<Point *> pts(base.originalPoints.size()) ;
    for(size_t i = 0 ; i < pts.size() ; i++)
        pts[i] = new Point(base.originalPoints[i]) ;

    Polygon basel(pts) ;
    basel.sampleSurface(num_points) ;

    std::vector<Point> ret ;
    for(size_t j = 0 ; j< basel.getInPoints().size() ; j++)
    {
        Point p = basel.getInPoint(j) ;
        p = p*rotateFromVector(Point(0,0,1), interpolatingPointAndTangent(0.).second);
        p += interpolationPoints.front() ;
        ret.push_back(p);
    }

    for(size_t j = 0 ; j< basel.getBoundingPoints().size() ; j++)
    {
        double numSlices = std::max(ceil(sweepNorm(basel.getBoundingPoint(j))/basel.getPerimeter()*basel.getBoundingPoints().size()), 2.) ;
        std::vector<double> dt = isoDistribute(basel.getBoundingPoint(j),numSlices) ;
        for(double i = 0 ; i< dt.size() ; i++)
        {
            std::pair<Point, Point> pt = interpolatingPointAndTangent(dt[i], basel.getBoundingPoint(j)) ;
            ret.push_back( pt.first );
        }
    }    
    
    for(size_t j = 0 ; j< basel.getInPoints().size() ; j++)
    {
        Point p = basel.getInPoint(j) ;
        p = p*rotateFromVector(Point(0,0,1), interpolatingPointAndTangent(1.).second) ;
        p += interpolationPoints.back() ;
        ret.push_back(p);
    }

    return ret ;
}

void PolygonPrism::sampleSurface(size_t num_points)
{
    sampleBoundingSurface(num_points);
    base.sampleSurface(num_points) ;

    double numSlices = std::max(round(axis.norm()/base.getPerimeter()*num_points), 2.) ;

    inPoints.resize((numSlices-2)*base.getInPoints().size());
    int piterator = 0  ;
    for(double i = 1 ; i< numSlices-1.1 ; i++)
    {
        for(size_t j = 0 ; j< base.getInPoints().size() ; j++)
        {
            Point p = base.getInPoint(j) ;
            p = rotationMatrix*p ;
            p += origin ;
            inPoints[piterator++] = new Point( axis * i/(numSlices-1.) + p );
        }
    }
}

void LoftedPolygonPrism::sampleSurface(size_t num_points)
{
//     num_points *=1.5 ;
    sampleBoundingSurface(num_points*1.5);
    base.sampleSurface(num_points) ;

    std::vector<Point> newPoints ;

    for(size_t j = 0 ; j< base.getInPoints().size() ; j++)
    {
        double numSlices = std::max(ceil(sweepNorm(base.getInPoint(j))/base.getPerimeter()*base.getBoundingPoints().size()), 2.);
        std::vector<double> dt = isoDistribute(base.getInPoint(j),numSlices) ;
        for(size_t i = 1 ; i< dt.size()-1 ; i++)
        {
            std::pair<Point, Point> pt = interpolatingPointAndTangent(dt[i], base.getInPoint(j)) ;
            newPoints.push_back( pt.first );
        }
    }

    inPoints.resize(newPoints.size());
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        inPoints[i] = new Point(newPoints[i]) ;
}

bool PolygonPrism::in(const Point & v) const
{
    Point tpoint(v) ;
    tpoint -= origin ;
    tpoint = inverse3x3Matrix(rotationMatrix)*tpoint;

    if(tpoint.getZ() < 0 || tpoint.getZ() > axis.norm())
        return false ;
    tpoint.getZ() = 0 ;

    return base.in(tpoint) ;
}

bool LoftedPolygonPrism::in(const Point & v) const
{
    double coordinate ; 
    std::pair<Point,Point> toCenter = projectToCenterLine(v,&coordinate) ;
    
    if(squareDist3D(v, toCenter.first) > base.getRadius()*base.getRadius())
        return false ;
    
    if(std::abs((v-toCenter.first)*toCenter.second) > 1e-6)
        return false ;
    
    Point tpoint(v) ;
    tpoint -= toCenter.first ;
    tpoint = tpoint*rotateToVector(Point(0,0,1), toCenter.second);
    tpoint.getZ() = 0 ;

    return base.in(tpoint) ;
    
}

double PolygonPrism::area() const
{
    return base.getPerimeter()*axis.norm() +2.*base.area();
}

double LoftedPolygonPrism::area() const
{
    return base.getPerimeter()*sweepNorm() +2.*base.area();
}

std::vector<double> LoftedPolygonPrism::isoDistribute(int npoints) const
{
    double delta = sweepNorm()/(npoints-1) ;
    std::vector<double> ret ;
    ret.push_back(0.);
    double dt = 0.005 ;
    double d = 0 ;
    double i = 0 ;
    for( ; i < 1.-dt ; i += dt)
    {
        d += dist(interpolatingPointAndTangent(i).first, interpolatingPointAndTangent(i+dt).first) ;
        if(d >= delta)
        {
            ret.push_back(i);
            d = 0 ;
        }
    }
    ret.push_back(1.);
    return ret ;
}

std::vector<double> LoftedPolygonPrism::isoDistribute(const Point & offset, int npoints) const
{
    double delta = sweepNorm(offset)/(npoints-1) ;
    std::vector<double> ret ;
    ret.push_back(0.);
    double dt = 0.005 ;
    double d = 0 ;
    double i = 0 ;
    for( ; i < 1.-dt ; i += dt)
    {
        d += dist(interpolatingPointAndTangent(i,offset).first, interpolatingPointAndTangent(i+dt,offset).first) ;
        if(d >= delta)
        {
            ret.push_back(i);
            d = 0 ;
        }
    }
    ret.push_back(1.);
    return ret ;
}

double LoftedPolygonPrism::sweepNorm(double end) const 
{
    double d = 0 ;
    double delta = 0.005 ;
    for(double i = 0 ; i < end ; i += delta)
    {
        d += dist(interpolatingPointAndTangent(i).first, interpolatingPointAndTangent(i+delta).first) ;
    }
//     exit(0) ;
    return d ;
}

double LoftedPolygonPrism::sweepNorm( const Point & offset, double end) const 
{
    double d = 0 ;
    double delta = 0.005 ;
    double lastd = 0 ;
    for(double i = 0 ; i < end ; i += delta)
    {
        d += dist(interpolatingPointAndTangent(i, offset).first, interpolatingPointAndTangent(i+delta, offset).first) ;
        lastd = i+delta ;
    }
    d += dist(interpolatingPointAndTangent(lastd, offset).first, interpolatingPointAndTangent(end, offset).first) ;
    return d ;
}
          
double PolygonPrism::volume() const {
    return base.area()*axis.norm() ;
}

double LoftedPolygonPrism::volume() const {
    return base.area()*sweepNorm() ;
}

std::pair<Point, Point> LoftedPolygonPrism::getEndNormals() const
{
    std::pair<Point, Point> n0 = interpolatingPointAndTangent(0.) ;
    std::pair<Point, Point> n1 = interpolatingPointAndTangent(1.) ;
    if(n0.second*(vstart-n0.first) < 0)
        n0.second *= -1 ;
    if(n1.second*(vend-n1.first) < 0)
        n1.second *= -1 ;
    return std::make_pair( n0.second, n1.second) ;
}

Point LoftedPolygonPrism::projectTest(const Amie::Point& init, const std::pair< Amie::Point, Amie::Point >& toCenter) const
{
    if(squareDist3D(toCenter.first, interpolationPoints.front()) < POINT_TOLERANCE*POINT_TOLERANCE )
    {
        Point tpoint(init) ;
        tpoint -= interpolationPoints.front() ;
        tpoint = tpoint*rotateToVector(Point(0,0,1), toCenter.second);
        tpoint.getZ() = 0 ;
        if(base.in(tpoint))
        {
            tpoint = tpoint*rotateFromVector(Point(0,0,1), toCenter.second);
            tpoint += interpolationPoints.front() ; 
            return tpoint;
        }
        else
        {
            base.project(&tpoint);
            tpoint = tpoint*rotateFromVector(Point(0,0,1), toCenter.second);
            tpoint += interpolationPoints.front() ; 
            return tpoint;
        }
    } 
    if(squareDist3D(toCenter.first, interpolationPoints.back()) < POINT_TOLERANCE*POINT_TOLERANCE)
    {
        Point tpoint(init) ;
        tpoint -= interpolationPoints.back() ;
        
        tpoint = tpoint*rotateToVector(Point(0,0,1), toCenter.second);
        tpoint.getZ() = 0 ;
        if(base.in(tpoint))
        {
            tpoint = tpoint*rotateFromVector(Point(0,0,1), toCenter.second);
            tpoint += interpolationPoints.back() ; 
            return tpoint;
        }
        else
        {
            base.project(&tpoint);
            tpoint = tpoint*rotateFromVector(Point(0,0,1), toCenter.second);
            tpoint += interpolationPoints.back() ; 
            return tpoint;
        }
    }
    
    Point tpoint(init) ;
    tpoint -= toCenter.first ;
    tpoint = tpoint*rotateToVector(Point(0,0,1), toCenter.second);
    tpoint.getZ() = 0 ;
    base.project(&tpoint);
    tpoint = tpoint*rotateFromVector(Point(0,0,1), toCenter.second);
    tpoint += toCenter.first ; 
    return tpoint ;
}

void LoftedPolygonPrism::project(Point * init) const
{
     init->set(projectTest(*init, projectToCenterLine(*init))) ;
}

std::pair<Point,Point> LoftedPolygonPrism::projectToCenterLine(const Point& p, double * coordinate) const
{
    std::map<double, double> distances ;
    double nsegs = (interpolationPoints.size()-1)*2 ;
    for(size_t i = 1 ; i < nsegs ; i++)
    {
        std::pair<Point, Point> pp = interpolatingPointAndTangent((double)i/nsegs) ;
        distances[squareDist3D(pp.first,p)] = (double)i/nsegs ;
    }
    double closestCoordinateDistance = distances.begin()->first ;
    double closestCoordinate = distances.begin()->second ;
    double minDistance = closestCoordinateDistance ;
    double bestCoordinate = closestCoordinate ;

    // step 2, Newton descent to get accurate result
    double previousClosestCoordinate = 0.5 ;
    int numit = 0 ;
    while(std::abs(closestCoordinate-previousClosestCoordinate) > 1e-4 || !numit)
    {
        double dh = 1e-6 ;
        double distDerivative = (squareDist3D(interpolatingPointAndTangent(closestCoordinate+dh).first,p)- squareDist3D(interpolatingPointAndTangent(closestCoordinate-dh).first,p))/(2.*dh) ;
        if(std::abs(distDerivative) < POINT_TOLERANCE*POINT_TOLERANCE)
            break ;
        double distSecondDerivative = (squareDist3D(interpolatingPointAndTangent(closestCoordinate+dh).first,p) -2.*closestCoordinateDistance + squareDist3D(interpolatingPointAndTangent(closestCoordinate-dh).first,p))/(dh*dh) ;
        if(std::abs(distSecondDerivative) < POINT_TOLERANCE*POINT_TOLERANCE)
            break ;
        closestCoordinate = closestCoordinate - distDerivative/distSecondDerivative ;
        if(closestCoordinate > 1.)
            closestCoordinate = 1 ;
        if(closestCoordinate < 0)
            closestCoordinate = 0 ;
        closestCoordinateDistance = squareDist3D(interpolatingPointAndTangent(closestCoordinate).first,p) ;
        if(closestCoordinateDistance < minDistance)
        {
            minDistance = closestCoordinateDistance ;
            bestCoordinate = closestCoordinate ;
        }
        if(numit++ > 8)
        { 
            closestCoordinate = bestCoordinate ;
            break ;
        }
        
    }
    if(closestCoordinate > 1.)
        closestCoordinate = 1 ;
    if(closestCoordinate < 0)
        closestCoordinate = 0 ;
    if(coordinate)
        *coordinate = closestCoordinate ;
    return interpolatingPointAndTangent(closestCoordinate) ;
}

void PolygonPrism::project(Point * init) const
{
    Point tpoint(*init) ;
    tpoint -= origin ;
    tpoint = inverse3x3Matrix(rotationMatrix)*tpoint;

    double originalz = tpoint.getZ() ;
    double anorm = axis.norm() ;
    tpoint.getZ() = 0 ;
    std::map<double, Point> findIt ;
    if(originalz > anorm || originalz < 0)
    {
        if(base.in(tpoint) )
        {
            if(originalz > anorm*.99)
            {
                tpoint = rotationMatrix*tpoint ;
                tpoint += origin + axis;
                findIt[squareDist3D(*init, tpoint)] = tpoint ;
            }
            else
            {
                tpoint = rotationMatrix*tpoint ;
                tpoint += origin ;
                findIt[squareDist3D(*init, tpoint)] = tpoint ;
            }
        }
        else
        {
            base.project(&tpoint);
            if(originalz > anorm*.99)
            {
                tpoint = rotationMatrix*tpoint ;
                tpoint += origin + axis;
                findIt[squareDist3D(*init, tpoint)] = tpoint ;
            }
            else
            {
                
                tpoint = rotationMatrix*tpoint ;
                tpoint += origin ;
                findIt[squareDist3D(*init, tpoint)] = tpoint ;
            }
        }
    }
    else
    {    
        base.project(&tpoint);
        tpoint = rotationMatrix*tpoint ;
        tpoint += origin + axis*originalz/anorm;
        findIt[squareDist3D(*init, tpoint)] = tpoint ;
    }

    init->set(findIt.begin()->second) ;

}

double PolygonPrism::getRadius() const
{
    double rb = base.getRadius() ;
    double ra = axis.norm() ;
    return sqrt(ra*ra+rb*rb) ;
}

double LoftedPolygonPrism::getRadius() const
{
    std::vector<Point> bb = getBoundingBox() ;
    return 0.5*dist(bb[0], bb[7]) ;
}

SpaceDimensionality PolygonPrism::spaceDimensions() const
{
    return SPACE_THREE_DIMENSIONAL ;
}

SpaceDimensionality LoftedPolygonPrism::spaceDimensions() const
{
    return SPACE_THREE_DIMENSIONAL ;
}

std::vector<Point> LoftedPolygonPrism::getBoundingBox() const
{
    std::valarray<Point> baseCopy = base.originalPoints;
    for(size_t i = 0 ; i < baseCopy.size() ; i++ )
    {
        baseCopy[i] = baseCopy[i]*rotateToVector(Point(0,0,1), interpolationPoints.front());

        baseCopy[i] += interpolationPoints.front() ;
    }

    double maxx = baseCopy[0].getX() ;
    double minx = baseCopy[0].getX() ;
    double maxy = baseCopy[0].getY() ;
    double miny = baseCopy[0].getY() ;
    double maxz = baseCopy[0].getZ() ;
    double minz = baseCopy[0].getZ() ;    
    
    for(size_t i = 1 ; i < baseCopy.size() ; i++ )
    {
        maxx = std::max(maxx,baseCopy[i].getX()) ;
        minx = std::min(minx,baseCopy[i].getX()) ;
        maxy = std::max(maxy,baseCopy[i].getY()) ;
        miny = std::min(miny,baseCopy[i].getY()) ;
        maxz = std::max(maxz,baseCopy[i].getZ()) ;
        minz = std::min(minz,baseCopy[i].getZ()) ;
    }
    
    for(double i = 0 ; i < interpolationPoints.size()+1 ;  i++)
    {
        std::valarray<Point> baseCopy = base.originalPoints;
        for(size_t j = 0 ; j < baseCopy.size() ; j++ )
        {
            baseCopy[j] = interpolatingPointAndTangent(i/interpolationPoints.size(), baseCopy[j]).first ;
        }
        for(size_t j = 0 ; j < baseCopy.size() ; j++ )
        {
            maxx = std::max(maxx,baseCopy[j].getX()) ;
            minx = std::min(minx,baseCopy[j].getX()) ;
            maxy = std::max(maxy,baseCopy[j].getY()) ;
            miny = std::min(miny,baseCopy[j].getY()) ;
            maxz = std::max(maxz,baseCopy[j].getZ()) ;
            minz = std::min(minz,baseCopy[j].getZ()) ;
        }

    }
    
    std::vector<Point> ret ;
    ret.push_back(Point(maxx, maxy, maxz)) ;
    ret.push_back(Point(maxx, maxy, minz)) ;
    ret.push_back(Point(maxx, miny, maxz)) ;
    ret.push_back(Point(maxx, miny, minz)) ;
    ret.push_back(Point(minx, maxy, maxz)) ;
    ret.push_back(Point(minx, maxy, minz)) ;
    ret.push_back(Point(minx, miny, maxz)) ;
    ret.push_back(Point(minx, miny, minz)) ;
    return ret ;    
}

std::vector<Point> PolygonPrism::getBoundingBox() const
{
    std::valarray<Point> baseCopy = base.originalPoints;
    for(size_t i = 0 ; i < baseCopy.size() ; i++ )
    {
        baseCopy[i] = rotationMatrix*baseCopy[i];
        baseCopy[i] += origin ;
    }

    double maxx = baseCopy[0].getX() ;
    double minx = baseCopy[0].getX() ;
    double maxy = baseCopy[0].getY() ;
    double miny = baseCopy[0].getY() ;
    double maxz = baseCopy[0].getZ() ;
    double minz = baseCopy[0].getZ() ;

    for(size_t i = 1 ; i < baseCopy.size() ; i++ )
    {
        maxx = std::max(maxx,baseCopy[i].getX()) ;
        minx = std::min(minx,baseCopy[i].getX()) ;
        maxy = std::max(maxy,baseCopy[i].getY()) ;
        miny = std::min(miny,baseCopy[i].getY()) ;
        maxz = std::max(maxz,baseCopy[i].getZ()) ;
        minz = std::min(minz,baseCopy[i].getZ()) ;
    }
    for(size_t i = 0 ; i < baseCopy.size() ; i++ )
        baseCopy[i] += axis ;

    for(size_t i = 0 ; i < baseCopy.size() ; i++ )
    {
        maxx = std::max(maxx,baseCopy[i].getX()) ;
        minx = std::min(minx,baseCopy[i].getX()) ;
        maxy = std::max(maxy,baseCopy[i].getY()) ;
        miny = std::min(miny,baseCopy[i].getY()) ;
        maxz = std::max(maxz,baseCopy[i].getZ()) ;
        minz = std::min(minz,baseCopy[i].getZ()) ;
    }

    std::vector<Point> ret ;
    ret.push_back(Point(maxx, maxy, maxz)) ;
    ret.push_back(Point(maxx, maxy, minz)) ;
    ret.push_back(Point(maxx, miny, maxz)) ;
    ret.push_back(Point(maxx, miny, minz)) ;
    ret.push_back(Point(minx, maxy, maxz)) ;
    ret.push_back(Point(minx, maxy, minz)) ;
    ret.push_back(Point(minx, miny, maxz)) ;
    ret.push_back(Point(minx, miny, minz)) ;
    return ret ;

}


}
