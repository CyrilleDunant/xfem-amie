// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_2D.h"
#include "../polynomial/vm_base.h"

namespace Amie
{

Triangle::Triangle() : ConvexGeometry(3)
{
    gType = TRIANGLE ;
    assert(this->size() == 3) ;

    boundingPoints[0] = new Point(0, 1) ;
    boundingPoints[1] = new Point(0, 0) ;
    boundingPoints[2] = new Point(1, 0) ;

    computeCenter() ;
    computeCircumCenter() ;
    radius = 0.471405;
    sqradius = radius*radius ;
    cachedarea = computeArea() ;
}

Triangle::Triangle( const Point & p0,  const Point & p1,  const Point & p2) : ConvexGeometry(3)
{
    gType = TRIANGLE ;

    assert(this->size() == 3) ;

    boundingPoints[0] = new Point(p0) ;
    boundingPoints[1] = new Point(p1) ;
    boundingPoints[2] = new Point(p2) ;


    if((p0.getZ() == p1.getZ())  && (p0.getZ() == p2.getZ()) &&  (p0.getZ() == 0))
    {
        if(!isTrigoOriented())
        {
            std::swap(boundingPoints[1], boundingPoints[2]) ;
            if(!isTrigoOriented())
            {
                std::cout << "arrrgh !" << std::endl ;
                boundingPoints[0]->print() ;
                boundingPoints[1]->print() ;
                boundingPoints[2]->print() ;
            }
        }
    }

    computeCircumCenter() ;
    computeCenter() ;

    radius = sqrt((squareDist2D(p1, circumCenter)+ squareDist2D(p0, circumCenter))*.5);
    sqradius = radius*radius ;
    cachedarea = computeArea() ;

}

void Triangle::setCenter(const Point & newCenter)
{
    Geometry::setCenter(newCenter);
    computeCircumCenter();
}


OrientedRectangle::OrientedRectangle() : ConvexGeometry(4)
{
    gType =PARALLELOGRAMME;
    assert(this->size() == 4);

    boundingPoints[0] = new Point(-1,1);
    boundingPoints[1] = new Point(-1,-1);
    boundingPoints[2] = new Point(1, -1);
    boundingPoints[3] = new Point(1, 1);
    computeCenter() ;
}

OrientedRectangle::OrientedRectangle( const Point & p0,  const Point & p1,  const Point & p2,  const Point & p3): ConvexGeometry(4)
{
    gType = PARALLELOGRAMME ;

    assert(this->size() == 4);

    boundingPoints[0] = new Point(p0) ;
    boundingPoints[1] = new Point(p1) ;
    boundingPoints[2] = new Point(p2) ;
    boundingPoints[3] = new Point(p3) ;

    if((p0.getZ() == p1.getZ())  && (p0.getZ() == p2.getZ()) &&  (p0.getZ() == 0))
    {
        if(!isTrigoOriented())
        {
            std::swap(boundingPoints[1], boundingPoints[2]) ;
        }
    }

    computeCircumCenter() ;
    computeCenter() ;

    if((p0.getZ() == p1.getZ())  && (p0.getZ() == p2.getZ()) &&  (p0.getZ() == 0))
    {
        if(!this->in(this->getCenter()))
        {
            assert(false) ;
        }
    }

    radius = (squareDist2D(p1, circumCenter) + squareDist2D(p0, circumCenter) + squareDist2D(p2, circumCenter)+squareDist2D(p3, circumCenter))/4.;

}

OrientedRectangle::OrientedRectangle( const Point *p0,  const Point *p1,  const Point *p2,  const Point *p3): ConvexGeometry(4)
{
    gType = PARALLELOGRAMME ;

    assert(this->size() == 4);

    boundingPoints[0] = new Point(*p0) ;
    boundingPoints[1] = new Point(*p1) ;
    boundingPoints[2] = new Point(*p2) ;
    boundingPoints[3] = new Point(*p3) ;

    if((p0->getZ() == p1->getZ())  && (p0->getZ() == p2->getZ()) &&  (p0->getZ() == 0))
    {
        if(!isTrigoOriented())
        {
            std::swap(boundingPoints[1], boundingPoints[2]) ;
        }
    }

    computeCircumCenter() ;
    computeCenter() ;

    radius = (squareDist2D(*p1, circumCenter) + squareDist2D(*p0, circumCenter) + squareDist2D(*p2, circumCenter)+squareDist2D(*p3, circumCenter))/4.;

}

void OrientedRectangle::computeCenter()
{
    for(size_t i = 0 ; i < this->size() ; i++)
    {
        this->center += this->getPoint(i) ;
    }

    this->center = this->center/this->size() ;
}

double OrientedRectangle::getRadius() const
{
    return sqrt(radius) ;
}

Point OrientedRectangle::getCircumCenter() const
{
    return this->circumCenter ;
}

std::vector<Point> OrientedRectangle::getBoundingBox() const
{
    double maxx = boundingPoints[0]->getX() ;
    double minx = boundingPoints[0]->getX() ;
    double maxy = boundingPoints[0]->getY() ;
    double miny = boundingPoints[0]->getY() ;

    for(size_t i = boundingPoints.size()/4 ; i < boundingPoints.size() ; i += boundingPoints.size()/4 )
    {
        double x = boundingPoints[i]->getX() ;
        double y = boundingPoints[i]->getY() ;
        if(x > maxx)
            maxx = x ;
        if(x < minx)
            minx = x ;
        if(y > maxy)
            maxy = y ;
        if(y < miny)
            miny = y ;
    }

    std::vector<Point> ret ;
    ret.push_back(Point(minx,miny)) ;
    ret.push_back(Point(maxx,miny)) ;
    ret.push_back(Point(maxx,maxy)) ;
    ret.push_back(Point(minx,maxy)) ;

    return  ret ;
}

void OrientedRectangle::computeCircumCenter()
{
    if (fabs(boundingPoints[1]->getY()-boundingPoints[0]->getY()) < 20*POINT_TOLERANCE)
    {
        double m2 = - (boundingPoints[2]->getX()-boundingPoints[1]->getX()) / (boundingPoints[2]->getY()-boundingPoints[1]->getY());
        double mx2 = (boundingPoints[1]->getX() + boundingPoints[2]->getX()) / 2.0;
        double my2 = (boundingPoints[1]->getY() + boundingPoints[2]->getY()) / 2.0;
        double xc = (boundingPoints[1]->getX() + boundingPoints[0]->getX()) / 2.0;
        double yc = fma(m2, (xc - mx2), my2);

        circumCenter = Point(xc, yc) ;
    }
    else if (fabs(boundingPoints[2]->getY()-boundingPoints[1]->getY()) < 20*POINT_TOLERANCE)
    {
        double m1 = - (boundingPoints[1]->getX()-boundingPoints[0]->getX()) / (boundingPoints[1]->getY()-boundingPoints[0]->getY());
        double mx1 = (boundingPoints[0]->getX() + boundingPoints[1]->getX()) / 2.0;
        double my1 = (boundingPoints[0]->getY() + boundingPoints[1]->getY()) / 2.0;
        double xc = (boundingPoints[2]->getX() + boundingPoints[1]->getX()) / 2.0;
        double yc = fma(m1, (xc - mx1), my1);

        circumCenter = Point(xc, yc) ;
    }
    else
    {
        double m1 = - (boundingPoints[1]->getX()-boundingPoints[0]->getX()) / (boundingPoints[1]->getY()-boundingPoints[0]->getY());
        double m2 = - (boundingPoints[2]->getX()-boundingPoints[1]->getX()) / (boundingPoints[2]->getY()-boundingPoints[1]->getY());
        double mx1 = (boundingPoints[0]->getX() + boundingPoints[1]->getX()) / 2.0;
        double mx2 = (boundingPoints[1]->getX() + boundingPoints[2]->getX()) / 2.0;
        double my1 = (boundingPoints[0]->getY() + boundingPoints[1]->getY()) / 2.0;
        double my2 = (boundingPoints[1]->getY() + boundingPoints[2]->getY()) / 2.0;
        double xc = fma(m1, mx1, fma(- m2 , mx2, my2 - my1)) / (m1 - m2);
        double yc = fma(m1, (xc - mx1), my1);

        circumCenter = Point(xc, yc) ;
    }
}

bool OrientedRectangle::inCircumCircle(const Point p) const
{
    return  squareDist2D(circumCenter, p) < radius  ;
}

bool OrientedRectangle::inCircumCircle(const Point *p) const
{
    return  squareDist2D(circumCenter, (*p)) < radius  ;
}

bool OrientedRectangle::is1D() const
{
    return false ;
}

double OrientedRectangle::area() const
{
    assert(this->boundingPoints.size() == 4) ;

    Segment s0((*boundingPoints[0]), (*boundingPoints[1])) ;
    Segment s1((*boundingPoints[0]), (*boundingPoints[2])) ;

    return ((s0.vector())^(s1.vector())).norm() ;
}


void OrientedRectangle::project(Point * p) const
{
    Segment s(getCenter(), *p) ;

    Point c0 = getBoundingPoint(0) ;
    Point c1 = getBoundingPoint(boundingPoints.size()/4) ;
    Point c2 = getBoundingPoint(boundingPoints.size()*2/4) ;
    Point c3 = getBoundingPoint(boundingPoints.size()*3/4) ;

    Segment s0(c0, c1) ;
    Segment s1(c1, c2) ;
    Segment s2(c2, c3) ;
    Segment s3(c3, c0) ;

    Point p0(*p) ;
    Point p1(*p) ;
    Point p2(*p) ;
    Point p3(*p) ;

    p0 = s0.project(p0) ;
    p1 = s1.project(p1) ;
    p2 = s2.project(p2) ;
    p3 = s3.project(p3) ;

    double d0 = squareDist3D(p0, *p) ;
    double d1 = squareDist3D(p1, *p) ;
    double d2 = squareDist3D(p2, *p) ;
    double d3 = squareDist3D(p3, *p) ;

    int i = 0 ;
    double d = d0 ;
    if(d1 < d)
    {
        d = d1 ;
        i = 1 ;
    }
    if(d2 < d)
    {
        d = d2 ;
        i = 2 ;
    }
    if(d3 < d)
    {
        d = d3 ;
        i = 3 ;
    }

    switch(i)
    {
    case 0:
        p->getX() = p0.getX() ;
        p->getY() = p0.getY() ;
        break ;
    case 1:
        p->getX() = p1.getX() ;
        p->getY() = p1.getY() ;
        break ;
    case 2:
        p->getX() = p2.getX() ;
        p->getY() = p2.getY() ;
        break ;
    case 3:
        p->getX() = p3.getX() ;
        p->getY() = p3.getY() ;
        break ;
    }

    return ;

}


bool OrientedRectangle::in(const Point &p) const
{

    Point c0 = getBoundingPoint(0) ;
    Point c1 = getBoundingPoint(boundingPoints.size()/4) ;
    Point c2 = getBoundingPoint(boundingPoints.size()*2/4) ;
    Point c3 = getBoundingPoint(boundingPoints.size()*3/4) ;

    Segment s0(c0, c1) ;
    Segment s1(c1, c2) ;
    Segment s2(c2, c3) ;
    Segment s3(c3, c0) ;

    Segment t(getCenter(), p) ;

    return !(t.intersects(s0) || t.intersects(s1) || t.intersects(s2) || t.intersects(s3)) ;

    bool in = false ;

    for (size_t i = 0, j  =  boundingPoints.size()-1; i <  boundingPoints.size(); j = i++)
    {
        if (
            (((boundingPoints[i]->getY() <= p.getY() )
              && (p.getY()<boundingPoints[j]->getY()))
             || ((boundingPoints[j]->getY() <= p.getY())
                 && (p.getY()<boundingPoints[i]->getY())))
            && (p.getX() < (boundingPoints[j]->getX() - boundingPoints[i]->getX()) * (p.getY() - boundingPoints[i]->getY()) / (boundingPoints[j]->getY() - boundingPoints[i]->getY()) + boundingPoints[i]->getX()))
            in = !in;
    }

    return in ;

}
std::vector<Point> OrientedRectangle::getSamplingBoundingPoints(double linearDensity) const
{
    std::vector<Point> ret ;

    Point  v0(*boundingPoints[0]) ;
    Point  v1(*boundingPoints[boundingPoints.size()/4]) ;
    Point  v2(*boundingPoints[2*boundingPoints.size()/4]) ;
    Point  v3(*boundingPoints[3*boundingPoints.size()/4]) ;

    ret.push_back(v0) ;
    size_t num_points = round(dist(v0,v1)*linearDensity+.25) ;
    for(size_t i = 1 ; i < num_points; i++)
    {
        ret.push_back(v0*i/num_points + v1*(1.-i/num_points)) ;
    }
    
    ret.push_back(v1) ;
    num_points = round(dist(v2,v1)*linearDensity+.25) ;
    for(size_t i = 1 ; i < num_points; i++)
    {
        ret.push_back(v1*i/num_points + v2*(1.-i/num_points)) ;
    }

    ret.push_back(v2) ;
    num_points = round(dist(v2,v3)*linearDensity+.25) ;
    for(size_t i = 1 ; i < num_points; i++)
    {
        ret.push_back(v2*i/num_points + v3*(1.-i/num_points)) ;
    }

    ret.push_back(v3) ;
    num_points = round(dist(v3,v0)*linearDensity+.25) ;
    for(size_t i = 1 ; i < num_points; i++)
    {
        ret.push_back(v3*i/num_points + v0*(1.-i/num_points)) ;
    }

    return ret ;
}


void OrientedRectangle::sampleBoundingSurface(double linearDensity)
{

    std::vector<Point> newPoints = getSamplingBoundingPoints(linearDensity) ;
    
    Point * v0(boundingPoints[0]) ;
    Point * v1(boundingPoints[boundingPoints.size()/4]) ;
    Point * v2(boundingPoints[2*boundingPoints.size()/4]) ;
    Point * v3(boundingPoints[3*boundingPoints.size()/4]) ;

    for(size_t i = 1 ; i < boundingPoints.size() ; i++)
    {
        if(i != boundingPoints.size()/4 && i != 2*boundingPoints.size()/4 && i != 3*boundingPoints.size()/4)
            delete boundingPoints[i] ;
    }
    boundingPoints.resize(newPoints.size()) ;

    boundingPoints[0] = v0 ;

    for(size_t i = 1 ; i < newPoints.size() ; i++)
    {
        if(newPoints[i] == *v0)
            boundingPoints[i] = v0 ;
        else if(newPoints[i] == *v1)
            boundingPoints[i] = v1 ;
        else if(newPoints[i] == *v2)
            boundingPoints[i] = v2 ;
        else if(newPoints[i] == *v3)
            boundingPoints[i] = v3 ;
        else
            boundingPoints[i] = new Point(newPoints[i]) ;
    }

}

void OrientedRectangle::sampleSurface(double linearDensity)
{
//	size_t n = 2*num_points ;
    sampleBoundingSurface(linearDensity) ;

    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    std::vector<Point> newInPoints ;

    size_t numberOfPointsAlongX = boundingPoints.size()/4 ;
    size_t numberOfPointsAlongY = boundingPoints.size()/4 ;
    size_t numberOfBoundingPoints = boundingPoints.size() ;

    Point c = getBoundingPoint(boundingPoints.size()/4)-getBoundingPoint(0) ;

    std::default_random_engine generator;
    std::uniform_real_distribution< double > distribution(-0.1/numberOfPointsAlongY, 0.1/numberOfPointsAlongY);
    for(size_t i = 1 ; i < numberOfPointsAlongX ; i++)
    {
        for(size_t j = 1 ; j < numberOfPointsAlongY ; j++)
        {
            double dj = (double) j / (double) numberOfPointsAlongY ;
            dj += distribution(generator) ;
            double di = distribution(generator) ;
            Point p1(*boundingPoints[i]) ;
            Point p2(*boundingPoints[numberOfBoundingPoints-numberOfPointsAlongX-i]) ;
            p1 = p1*dj + p2*(1.-dj) + c*di ;
            if(in(p1))
                newInPoints.push_back(p1) ;
        }
    }

    inPoints.resize(newInPoints.size()) ;
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        inPoints[i] = new Point(newInPoints[i]) ;

}

Triangle::Triangle( Point *p0,  Point *p1,  Point *p2): ConvexGeometry(3)
{
    gType = TRIANGLE ;

    assert(this->size() == 3) ;


    boundingPoints[0] = p0 ;
    boundingPoints[1] = p1 ;
    boundingPoints[2] = p2 ;


    if((p0->getZ() == p1->getZ())  && (p0->getZ() == p2->getZ()) &&  (p0->getZ() == 0))
    {
        if(!isTrigoOriented())
        {
            std::swap(boundingPoints[1], boundingPoints[2]) ;
            if(!isTrigoOriented())
                std::cout << "arrrgh !" << std::endl ;
        }
    }

    computeCircumCenter() ;
    computeCenter() ;

    radius = (dist(p1, &circumCenter)+ dist(p0, &circumCenter))/2.;
    sqradius = radius*radius ;
    cachedarea = computeArea() ;
}


void Triangle::computeCenter()
{
    for(size_t i = 0 ; i < this->size() ; i++)
    {
        this->center += this->getPoint(i) ;
    }

    this->center = this->center/this->size() ;
}

double Triangle::getRadius() const
{
    return radius ;
}

const Point & Triangle::getCircumCenter() const
{
    return this->circumCenter ;
}



std::vector<Point> Triangle::getBoundingBox() const
{
    std::vector<Point> box ;
    box.push_back(getCircumCenter()+Point(-getRadius(), getRadius())) ;
    box.push_back(getCircumCenter()+Point(getRadius(), getRadius())) ;
    box.push_back(getCircumCenter()+Point(getRadius(), -getRadius())) ;
    box.push_back(getCircumCenter()+Point(-getRadius(), -getRadius())) ;

    return box ;
}

void Triangle::computeCircumCenter()
{

    if (std::abs(getBoundingPoint(1).getY()-getBoundingPoint(0).getY()) < 100.*POINT_TOLERANCE)
    {
        double m2  =  (getBoundingPoint(1).getX()-getBoundingPoint(2).getX() ) / (getBoundingPoint(2).getY()-getBoundingPoint(1).getY());
        double mx2 = (getBoundingPoint(1).getX() + getBoundingPoint(2).getX()) ;
        double my2 = (getBoundingPoint(1).getY() + getBoundingPoint(2).getY()) ;
        double xc  = (getBoundingPoint(1).getX() + getBoundingPoint(0).getX()) ;
        double yc  = m2 * (xc - mx2) + my2;

        circumCenter.set(xc/2., yc/2.) ;
    }
    else if (std::abs(getBoundingPoint(2).getY()-getBoundingPoint(1).getY()) < 100.*POINT_TOLERANCE)
    {
        double m1  =  (getBoundingPoint(0).getX() - getBoundingPoint(1).getX() ) / (getBoundingPoint(1).getY()-getBoundingPoint(0).getY());
        double mx1 = (getBoundingPoint(0).getX() + getBoundingPoint(1).getX()) ;
        double my1 = (getBoundingPoint(0).getY() + getBoundingPoint(1).getY()) ;
        double xc  = (getBoundingPoint(2).getX() + getBoundingPoint(1).getX()) ;
        double yc  = m1 * (xc - mx1) + my1;

        circumCenter.set(xc/2., yc/2.) ;
    }
    else
    {
        double m1  = (getBoundingPoint(0).getX()-getBoundingPoint(1).getX()) / (getBoundingPoint(1).getY()-getBoundingPoint(0).getY());
        double m2  = (getBoundingPoint(1).getX()-getBoundingPoint(2).getX()) / (getBoundingPoint(2).getY()-getBoundingPoint(1).getY());
        double mx1 = (getBoundingPoint(0).getX() + getBoundingPoint(1).getX()) ;
        double mx2 = (getBoundingPoint(1).getX() + getBoundingPoint(2).getX()) ;
        double my1 = (getBoundingPoint(0).getY() + getBoundingPoint(1).getY()) ;
        double my2 = (getBoundingPoint(1).getY() + getBoundingPoint(2).getY()) ;
        double xc  = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
        double yc  = m1 * (xc - mx1) + my1;

        circumCenter.set(xc/2., yc/2.) ;
    }
}

bool Triangle::inCircumCircle(const Point & p) const
{
    if(p.getX() > circumCenter.getX()+1.01*radius)
        return false ;
    if(p.getX() < circumCenter.getX()-1.01*radius)
        return false ;
    if(p.getY() > circumCenter.getY()+1.01*radius)
        return false ;
    if(p.getY() < circumCenter.getY()-1.01*radius)
        return false ;

    if(squareDist2D(circumCenter, p) < .99*sqradius)
        return true ;

    double delta = POINT_TOLERANCE ;
    Point a(p) ;
    a.getX() += delta ;
    a.getY() += delta ;
    Point c(p) ;
    c.getX() += delta ;
    c.getY() -= delta ;
    Point e(p) ;
    e.getX() -= delta ;
    e.getY() += delta ;
    Point g(p) ;
    g.getX() -= delta ;
    g.getY() -= delta ;

    return  squareDist2D(circumCenter, a) < sqradius
            &&  squareDist2D(circumCenter, c) < sqradius
            &&  squareDist2D(circumCenter, e) < sqradius
            &&  squareDist2D(circumCenter, g) < sqradius;
    double x = circumCenter.getX() -p.getX() ;
    double y = circumCenter.getY() -p.getY() ;
    return  fma(x, x, y*y)< sqradius*(1. - 100.*POINT_TOLERANCE)  ;
}

bool Triangle::inCircumCircle(const Point *p) const
{
    if(p->getX() > circumCenter.getX()+1.01*radius)
        return false ;
    if(p->getX() < circumCenter.getX()-1.01*radius)
        return false ;
    if(p->getY() > circumCenter.getY()+1.01*radius)
        return false ;
    if(p->getY() < circumCenter.getY()-1.01*radius)
        return false ;

    if(squareDist2D(circumCenter, *p) < .99*sqradius)
        return true ;

    double delta = POINT_TOLERANCE ;
    Point a(*p) ;
    a.getX() += delta ;
    a.getY() += delta ;
    Point c(*p) ;
    c.getX() += delta ;
    c.getY() -= delta ;
    Point e(*p) ;
    e.getX() -= delta ;
    e.getY() += delta ;
    Point g(*p) ;
    g.getX() -= delta ;
    g.getY() -= delta ;

    return  squareDist2D(circumCenter, a) < sqradius
            &&  squareDist2D(circumCenter, c) < sqradius
            &&  squareDist2D(circumCenter, e) < sqradius
            &&  squareDist2D(circumCenter, g) < sqradius;
    double x = circumCenter.getX() -p->getX() ;
    double y = circumCenter.getY() -p->getY() ;
    return  fma(x, x, y*y) < sqradius*(1. - 100.*POINT_TOLERANCE)  ;
}


double Triangle::area() const
{
    return cachedarea ;
}


double Triangle::computeArea() 
{

    int pointsInTimePlane = this->boundingPoints.size()/timePlanes() ;

    return 0.5*std::abs((getBoundingPoint(0).getX()- getBoundingPoint(pointsInTimePlane/3).getX())*(getBoundingPoint(0).getY()- getBoundingPoint(2*pointsInTimePlane/3).getY()) - (getBoundingPoint(0)- getBoundingPoint(pointsInTimePlane/3)).getY()*(getBoundingPoint(0).getX()- getBoundingPoint(2*pointsInTimePlane/3).getX())) ;

}

void Triangle::project(Point * p) const
{
    //find closest vertex
    std::map<double, Point> vert ;
    vert[squareDist2D(*p, getBoundingPoint(0))] =  getBoundingPoint(0) ;
    vert[squareDist2D(*p, getBoundingPoint(getBoundingPoints().size()/(3*timePlanes())))] =  getBoundingPoint(getBoundingPoints().size()/(3*timePlanes())) ;
    vert[squareDist2D(*p, getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes())))] =  getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes())) ;
    if(vert.begin()->second == getBoundingPoint(0) )
    {
        Segment s0(getBoundingPoint(0), getBoundingPoint(getBoundingPoints().size()/(3*timePlanes()))) ;
        Segment s2( getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes())), getBoundingPoint(0)) ;
        std::map<double, Point> proj ;
        Point p0 = s0.project(*p) ;
        Point p2 = s2.project(*p) ;
        proj[squareDist2D(p0, *p)] = p0 ;
        proj[squareDist2D(p2, *p)] = p2 ;
        *p = proj.begin()->second ;
        return ;
    }
    else if(vert.begin()->second == getBoundingPoint(getBoundingPoints().size()/(3*timePlanes())))
    {
        Segment s0(getBoundingPoint(0), getBoundingPoint(getBoundingPoints().size()/(3*timePlanes()))) ;
        Segment s1(getBoundingPoint(getBoundingPoints().size()/(3*timePlanes())), getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes()))) ;
        std::map<double, Point> proj ;
        Point p0 = s0.project(*p) ;
        Point p1 = s1.project(*p) ;
        proj[squareDist2D(p0, *p)] = p0 ;
        proj[squareDist2D(p1, *p)] = p1 ;
        *p = proj.begin()->second ;
        return ;
    }
    
    Segment s1(getBoundingPoint(getBoundingPoints().size()/(3*timePlanes())), getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes()))) ;
    Segment s2( getBoundingPoint(2*getBoundingPoints().size()/(3*timePlanes())), getBoundingPoint(0)) ;
    std::map<double, Point> proj ;
    Point p1 = s1.project(*p) ;
    Point p2 = s2.project(*p) ;
    proj[squareDist2D(p1, *p)] = p1 ;
    proj[squareDist2D(p2, *p)] = p2 ;
    *p = proj.begin()->second ;
    return ;

}

bool Triangle::in(const Point &p) const
{

    if(!inCircumCircle(p))
        return false ;

    for (size_t i = 0; i <  getBoundingPoints().size(); i++)
    {
        if(p == getBoundingPoint(i))
        {
            return true ;
        }
    }

    Segment s(p, getCenter()) ;

    size_t npts = getBoundingPoints().size()/timePlanes() ;

    if(s.on(getBoundingPoint(0)) || s.on(getBoundingPoint(npts/3)) || s.on(getBoundingPoint(npts*2/3)))
        return false ;

    Segment sa(getBoundingPoint(0),getBoundingPoint(npts/3)) ;
    if(sa.intersects(s))
        return false ;
    Segment sb(getBoundingPoint(npts/3),getBoundingPoint(npts*2/3)) ;
    if(sb.intersects(s))
        return false ;
    Segment sc(getBoundingPoint(npts*2/3),getBoundingPoint(0)) ;
    if(sc.intersects(s))
        return false ;
    
    return true ;

}


std::vector<Point> Triangle::getSamplingBoundingPoints(double linearDensity) const
{
    std::vector<Point> ret ;

    int n = getBoundingPoints().size() ;

    Point v0 = *boundingPoints[0] ;
    Point v1 = *boundingPoints[n/3] ;
    Point v2 = *boundingPoints[2*n/3] ;

    ret.push_back(v0) ;
    size_t num_points = round(linearDensity*dist(v0,v1)) ;
    for(size_t i = 1 ; i < num_points ; i++)
    {
        ret.push_back(v0*i/num_points + v1*(1.-i/num_points)) ;
    }

    ret.push_back(v1) ;
    num_points = round(linearDensity*dist(v1,v2)) ;
    for(size_t i = 1 ; i < num_points ; i++)
    {
        ret.push_back(v1*i/num_points + v2*(1.-i/num_points)) ;
    }

    ret.push_back(v2) ;
    num_points = round(linearDensity*dist(v2,v0)) ;
    for(size_t i = 1 ; i < num_points ; i++)
    {
        ret.push_back(v2*i/num_points + v0*(1.-i/num_points)) ;
    }

    return ret ;
}

void Triangle::sampleBoundingSurface(double linearDensity) 
{
    std::vector<Point> newPoints = getSamplingBoundingPoints(linearDensity) ;
    
    Point * v0(boundingPoints[0]) ;
    Point * v1(boundingPoints[boundingPoints.size()/3]) ;
    Point * v2(boundingPoints[2*boundingPoints.size()/3]) ;

    for(size_t i = 1 ; i < boundingPoints.size() ; i++)
    {
        if(i != boundingPoints.size()/3 && i != 2*boundingPoints.size()/3 && i != 0)
            delete boundingPoints[i] ;
    }
    boundingPoints.resize(newPoints.size()) ;

    boundingPoints[0] = v0 ;
    for(size_t i = 1 ; i < newPoints.size() ; i++)
    {
        if(newPoints[i] == *v0)
            boundingPoints[i] = v0 ;
        else if(newPoints[i] == *v1)
            boundingPoints[i] = v1 ;
        else if(newPoints[i] == *v2)
            boundingPoints[i] = v2 ;
        else
            boundingPoints[i] = new Point(newPoints[i]) ;
    }

}

void Triangle::sampleSurface(double linearDensity)
{
    size_t num_points = round((dist(boundingPoints[0], boundingPoints[boundingPoints.size()/3]) + 
                         dist(boundingPoints[boundingPoints.size()/3],boundingPoints[2*boundingPoints.size()/3]) +
                         dist(boundingPoints[0], boundingPoints[2*boundingPoints.size()/3]) )*linearDensity) ;
    num_points += 3-num_points%3 ;

    sampleBoundingSurface(linearDensity*3) ;

    std::vector<Point> newPoints ;

    size_t end_i = boundingPoints.size()/3 ;
    Point p1 = getBoundingPoint( 0 ) ;
    Point p2 = getBoundingPoint( end_i )-p1 ;
    Point p3 = getBoundingPoint( end_i*2 )-p1 ;
    double d = (p2.norm()+p3.norm())*0.01 ;


    for(int i = 0 ; i < ((int)num_points-1) ; i+=1)
    {
        for(int j = 0 ; j < i+1 ; j+=1)
        {
            double fact = (double)(j+1)/(double)(i+1) ;
            double xrand = ((double)rand()/(double)(RAND_MAX)*2.-1.)*d ;
            double yrand = ((double)rand()/(double)(RAND_MAX)*2.-1.)*d ;
            if(this->in(getBoundingPoint(i+1)*(1.-fact) + getBoundingPoint(end_i-1-i)*fact+Point(xrand, yrand)))
                newPoints.push_back(getBoundingPoint(i+1)*(1.-fact) + getBoundingPoint(end_i-1-i)*fact+Point(xrand, yrand)) ;
        }
    }
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    this->inPoints.resize(newPoints.size()) ;

    for(size_t i = 0 ; i < inPoints.size() ; i++)
    {
        inPoints[i] = new Point(newPoints[i]) ;
    }

}

void Triangle::print() const
{
    std::cout << "inPoints" << std::endl ;

    for(size_t i = 0 ; i < inPoints.size() ; i++ )
    {
        getInPoint(i).print() ;
    }
    std::cout << std::endl ;

    std::cout << "boundingPoints" << std::endl ;

    for(size_t i = 0 ; i < boundingPoints.size() ; i++ )
    {
        getBoundingPoint(i).print() ;
    }
    std::cout << std::endl ;

}

Rectangle::Rectangle(double x, double y, double originX, double originY) : ConvexGeometry(4), size_y(y), size_x(x)
{
    gType = RECTANGLE ;
    this->center = Point(originX,originY) ;
    topLeft = Point(originX-0.5*x, originY+0.5*y) ;
    boundingPoints[0] = new Point(topLeft) ;
    topRight = Point(originX+0.5*x, originY+0.5*y) ;
    boundingPoints[1] = new Point(topRight) ;
    bottomRight = Point(originX+0.5*x, originY-0.5*y) ;
    boundingPoints[2] = new Point(bottomRight) ;
    bottomLeft =  Point(originX-0.5*x, originY-0.5*y) ;
    boundingPoints[3] = new Point(bottomLeft) ;
}

Rectangle::Rectangle(std::vector<Point> box) : ConvexGeometry(4)
{
    gType = RECTANGLE ;
    topLeft = Point(box[0]) ;
    topRight = Point(box[1]) ;
    bottomRight = Point(box[2]) ;
    bottomLeft = Point(box[3]) ;
    boundingPoints[0] = new Point(topLeft) ;
    boundingPoints[1] = new Point(topRight) ;
    boundingPoints[2] = new Point(bottomRight) ;
    boundingPoints[3] = new Point(bottomLeft) ;
    center = (box[0]+box[1]+box[2]+box[3])*0.25 ;
    size_x = std::abs(box[1].getX()-box[0].getX()) ;
    size_y = std::abs(box[2].getY()-box[1].getY()) ;
}

Rectangle::Rectangle(double x, double y, const Point &center) :  ConvexGeometry(4), size_y(y), size_x(x)
{
    gType = RECTANGLE ;
    this->center = center ;
    topLeft = Point(center.getX()-0.5*x, center.getY()+0.5*y) ;
    boundingPoints[0] = new Point(topLeft) ;
    topRight = Point(center.getX()+0.5*x, center.getY()+0.5*y) ;
    boundingPoints[1] = new Point(topRight) ;
    bottomRight = Point(center.getX()+0.5*x, center.getY()-0.5*y) ;
    boundingPoints[2] = new Point(bottomRight) ;
    bottomLeft =  Point(center.getX()-0.5*x, center.getY()-0.5*y) ;
    boundingPoints[3] = new Point(bottomLeft) ;
}

Rectangle::Rectangle() :  ConvexGeometry(4), size_y(2), size_x(2)
{
    gType = RECTANGLE ;
    this->center = Point(0,0) ;
    topLeft = Point(-1, 1) ;
    boundingPoints[0] = new Point(topLeft) ;
    topRight = Point(1,1) ;
    boundingPoints[1] = new Point(topRight) ;
    bottomRight = Point(1, -1) ;
    boundingPoints[2] = new Point(bottomRight) ;
    bottomLeft = Point(-1, -1) ;
    boundingPoints[3] = new Point(bottomLeft) ;


}

std::vector<Point> Rectangle::getBoundingBox() const
{
    std::vector<Point> box ;
    box.push_back(topLeft) ;
    box.push_back(topRight) ;
    box.push_back(bottomRight) ;
    box.push_back(bottomLeft) ;

    return box ;
}

void Rectangle::computeCenter()
{
    for(size_t i = 0 ; i < this->size() ; i++)
        this->center += this->getPoint(i) ;

    if(this->size() != 0)
        this->center = this->center/this->size() ;
}

double  Rectangle::getRadius() const
{
    return sqrt(width()*width()*0.25 + height()*height()*0.25) ;
}

bool Rectangle::in(const Point & p) const
{
    if(p.getX() < getCenter().getX() - 0.5*width()-POINT_TOLERANCE)
        return false ;
    if(p.getX()  > getCenter().getX() + 0.5*width()+POINT_TOLERANCE)
        return false ;
    if(p.getY() > getCenter().getY() + 0.5*height()+POINT_TOLERANCE)
        return false ;
    if(p.getY()  < getCenter().getY() - 0.5*height()-POINT_TOLERANCE)
        return false ;

    return true ;
}

double Rectangle::area() const
{
    return this->size_x*this->size_y ;
}

double Rectangle::width() const
{
    return size_x ;
}

double Rectangle::height() const
{
    return size_y ;
}


void Rectangle::project(Point * p) const
{
    Segment A(topLeft, bottomLeft) ;

    Segment B(topLeft, topRight) ;

    Segment C(topRight, bottomRight) ;

    Segment D(bottomRight, bottomLeft) ;

    std::map<double, Point> tries ;
    Point p0 = A.project(*p) ;
    Point p1 = B.project(*p) ;
    Point p2 = C.project(*p) ;
    Point p3 = D.project(*p) ;
    tries[dist(p, &p0)] = p0 ;
    tries[dist(p, &p1)] = p1 ;
    tries[dist(p, &p2)] = p2 ;
    tries[dist(p, &p3)] = p3 ;
    *p = tries.begin()->second ;
}

std::vector<Point> Rectangle::getSamplingBoundingPoints(double linearDensity) const
{
    double perimeter = 2*(size_x+size_y) ;

    size_t num_points = round(linearDensity*perimeter) ;
    double distanceBetweenPoints = perimeter/num_points ;
    std::vector<Point> ret ;

    double numberOfPointsAlongX = static_cast<size_t>(std::ceil(size_x/distanceBetweenPoints) + 1);
    double distanceBetweenPointsAlongX = size_x/(numberOfPointsAlongX-1) ;

    double numberOfPointsAlongY = static_cast<size_t>(std::ceil(size_y/distanceBetweenPoints) + 1);
    double distanceBetweenPointsAlongY = size_y/(numberOfPointsAlongY-1) ;

    num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;

    for (size_t i = 0 ; i < numberOfPointsAlongY; i++)
    {
// 		double randx=((2.*rand()/(RAND_MAX+1.0))-1.)*0.15*(size_x/numberOfPointsAlongX) ;
        double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        if(i == 0 || i == numberOfPointsAlongY-1)
            randy = 0 ;
        ret.push_back(Point(center.getX()-0.5*size_x, center.getY() + 0.5*size_y - i*distanceBetweenPointsAlongY+ randy)) ;
    }
    for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
    {
        double randx= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        if(i == numberOfPointsAlongX-1)
            randx = 0 ;
        ret.push_back(Point( center.getX()-0.5*size_x+i*distanceBetweenPointsAlongX+ randx,
                             getCenter().getY()-0.5*size_y));
    }
    for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
    {
        double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        if(i == numberOfPointsAlongY-1)
            randy = 0 ;
        ret.push_back(Point(center.getX()+0.5*size_x,
                            center.getY()-0.5*size_y+i*distanceBetweenPointsAlongY+ randy));
    }
    for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
    {
        double randx=  0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        assert(2*numberOfPointsAlongY+numberOfPointsAlongX+i-3< num_points) ;
        ret.push_back(Point(
                          center.getX() + 0.5*size_x - i*distanceBetweenPointsAlongX +randx ,
                          center.getY() + 0.5*size_y)) ;
    }

    return ret ;
}

void Rectangle::sampleBoundingSurface(double linearDensity)
{
    //	assert(num_points%4 == 0) ;
    double perimeter = 2.*(size_x+size_y) ;
    linearDensity = std::max(linearDensity, 2./std::min(size_x,size_y)) ;
    size_t num_points = round(perimeter*linearDensity) ;

    double distanceBetweenPointsx = std::min(perimeter/num_points, size_x) ;
    double distanceBetweenPointsy = std::min(perimeter/num_points, size_y) ;
    double dy = distanceBetweenPointsy ;
    double dx = distanceBetweenPointsx ;
    distanceBetweenPointsy = std::min(dx*1.5, dy) ;
    distanceBetweenPointsx = std::min(dy*1.5, dx) ;

    numberOfPointsAlongX = round(1.5*static_cast<size_t>(std::ceil(size_x/distanceBetweenPointsx) + 1));
    double distanceBetweenPointsAlongX = size_x/(this->numberOfPointsAlongX-1) ;

    numberOfPointsAlongY = round(1.5*static_cast<size_t>(std::ceil(size_y/distanceBetweenPointsy) + 1));
    double distanceBetweenPointsAlongY = size_y/(this->numberOfPointsAlongY-1) ;

    num_points = ((numberOfPointsAlongX)*2 + (numberOfPointsAlongY)*2 - 4) ;

    for(size_t i = 0 ; i < boundingPoints.size(); i++)
        delete boundingPoints[i] ;

    boundingPoints.resize(num_points) ;

    for (size_t i = 0 ; i < numberOfPointsAlongY; i++)
    {
// 		double randx=((2.*rand()/(RAND_MAX+1.0))-1.)*0.15*(size_x/numberOfPointsAlongX) ;
        double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        if(i == 0 || i == numberOfPointsAlongY-1)
            randy = 0 ;
        boundingPoints[i] = new Point(center.getX()-0.5*size_x, center.getY() + 0.5*size_y - i*distanceBetweenPointsAlongY+ randy) ;
    }
    for (size_t i = 1 ; i < numberOfPointsAlongX ; i++)
    {
        double randx= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        if(i == numberOfPointsAlongX-1)
            randx = 0 ;
        boundingPoints[numberOfPointsAlongY+i-1] = new Point( center.getX()-0.5*size_x + i*distanceBetweenPointsAlongX+ randx,
                getCenter().getY()-0.5*size_y);
    }
    for (size_t i = 1 ; i < numberOfPointsAlongY ; i++)
    {
        double randy= 0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        if(i == numberOfPointsAlongY-1)
            randy = 0 ;
        boundingPoints[numberOfPointsAlongX+numberOfPointsAlongY+i-2] = new Point(center.getX()+0.5*size_x,
                center.getY()-0.5*size_y+i*distanceBetweenPointsAlongY+ randy);
    }
    for (size_t i = 1 ; i < numberOfPointsAlongX-1 ; i++)
    {
        double randx=  0 ;//((2.*rand()/(RAND_MAX+1.0))-1.)*0.22*(size_y/numberOfPointsAlongY) ;
        assert(2*numberOfPointsAlongY+numberOfPointsAlongX+i-3< num_points) ;
        boundingPoints[2*numberOfPointsAlongY+numberOfPointsAlongX+i-3] = new Point(
            center.getX() + 0.5*size_x - i*distanceBetweenPointsAlongX +randx ,
            center.getY() + 0.5*size_y) ;
    }
    numberOfPointsAlongY *= .666 ;
    numberOfPointsAlongX *= .666 ;
}

void Rectangle::sampleSurface(double linearDensity)
{

//     if(std::max(size_x/size_y, size_y/size_x) < 10)
//     {
        sampleBoundingSurface(linearDensity*2.) ;
//     }
//     else if(std::max(size_x/size_y, size_y/size_x) < 60)
//     {
//         sampleBoundingSurface((size_t)round((double)num_points*0.5*std::max(size_x/size_y, size_y/size_x)/(M_PI))) ;
//     }
//     else
//     {
//         sampleBoundingSurface((size_t)round((double)num_points*0.2*std::max(size_x/size_y, size_y/size_x)/(M_PI))) ;
//     }

    size_t nip = static_cast<size_t>((numberOfPointsAlongX-2)*(numberOfPointsAlongY-2)) ;

    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    std::vector<Point *> newInPoints ;

    double distanceBetweenPointsAlongX = size_x/(numberOfPointsAlongX-1) ;
    double distanceBetweenPointsAlongY = size_y/(numberOfPointsAlongY-1) ;

    if(nip > 0)
    {
        for(size_t j = 0 ; j < numberOfPointsAlongY-2 ; j++)
        {
            for(size_t i = 0 ; i < numberOfPointsAlongX-2+(j+1)%2 ; i++)
            {

                double randx= ((2.*rand()/(RAND_MAX+1.0))-1.)*0.25*distanceBetweenPointsAlongX ;
                double randy= ((2.*rand()/(RAND_MAX+1.0))-1.)*0.25*distanceBetweenPointsAlongY ;

                newInPoints.push_back( new Point(center.getX() - 0.5*size_x + (double)(i+0.66)*distanceBetweenPointsAlongX+(double)((j)%2)*distanceBetweenPointsAlongX*.5-(double)((j+1)%2)*distanceBetweenPointsAlongX*.15+randx,
                                                 center.getY() - 0.5*size_y + (double)(j+1)*distanceBetweenPointsAlongY+ randy)) ;
            }
        }
    }

    inPoints.resize(newInPoints.size()) ;
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        inPoints[i] = newInPoints[i] ;
}

Circle::Circle(double r, double originX, double originY)
{
    gType = CIRCLE ;
    this->center = Point(originX, originY) ;
    this->radius = r ;
    this->sqradius = r*r ;
}

Circle::Circle(double r, const Point *center)
{
    gType = CIRCLE ;
    this->center = Point(*center) ;
    this->radius = r ;
    this->sqradius = r*r ;
}

Circle::Circle(double r, const Point & center)
{
    gType = CIRCLE ;
    this->center = center ;
    this->radius = r ;
    this->sqradius = r*r ;
}


void Circle::setRadius(double newr)
{
    double ratio = newr/(radius) ;


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

    this->radius = newr ;
    this->sqradius = newr*newr ;
}

void Circle::computeCenter()
{
}

std::vector<Point> Circle::getBoundingBox() const
{
    std::vector<Point> box ;
    box.push_back(getCenter()+Point(-getRadius(), getRadius())) ;
    box.push_back(getCenter()+Point(getRadius(), getRadius())) ;
    box.push_back(getCenter()+Point(getRadius(), -getRadius())) ;
    box.push_back(getCenter()+Point(-getRadius(), -getRadius())) ;

    return box ;
}

void Circle::project(Point * p) const
{
    if(squareDist2D(p, &getCenter() ) < POINT_TOLERANCE*POINT_TOLERANCE)
    {
        p->getX() +=getRadius() ;
        return ;
    }


    Line l(*p, *p-getCenter()) ;

    std::vector<Point> inter = l.intersection(this) ;
    *p = inter[0] ;
// 	if(inter.empty() || inter.size() == 1)
// 	{
// 		p->print() ;
// 		getCenter().print() ;
// 	}
    if((inter.size() == 2) && squareDist3D(inter[0], *p) > squareDist3D(inter[1], *p))
    {
       *p = inter[1] ;
    }
    
}

std::vector<Point> Circle::getSamplingBoundingPoints(double linearDensity) const
{
    std::vector<Point> ret ;

    size_t num_points = round(linearDensity*2.*M_PI*getRadius()) ;
    double angle = 2.*M_PI/ (num_points) ;
    double start = angle * (0.5-(double) std::rand()/(double) RAND_MAX) ;

    for (size_t i = 0 ; i < num_points ; i ++)
    {
        ret.push_back(Point(getRadius()*cos(start+angle*(double) i) + getCenter().getX(), getRadius()*sin(start+angle*(double) i) + getCenter().getY()));
    }

    return ret ;
}

std::vector<Point> Circle::getSamplingBoundingPointsOnArc(size_t num_points, const Point & start, const Point & finish, bool reverse) const
{
    std::vector<Point> ret ;
    if(num_points == 0)
        return ret ;
    Point init(start) ;
    project(&init) ;
    init -= getCenter() ;
    Point fin(finish) ;
    project(&fin) ;
    fin -= getCenter() ;
    double a = init.angle() ;
    double b = fin.angle() ;
    if(b > a)
    {
        a += 2.*M_PI ;
    }
    double angle = (a-b)/num_points ;
    if(reverse)
        angle = (2.*M_PI-(a-b))/num_points ;

    for (double i = 0 ; i< num_points ; i++)
    {
        Point newPoint(init.getX()*cos(i*angle)+init.getY()*sin(i*angle), -init.getX()*sin(i*angle)+init.getY()*cos(i*angle)) ;
        newPoint+= getCenter() ;
        ret.push_back(newPoint);
    }
    return ret ;
}

void Circle::sampleBoundingSurface(double linearDensity)
{
    size_t num_points = round(linearDensity*2.*M_PI*getRadius()) ;
    num_points = std::max(num_points, (size_t)8) ;
    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        delete boundingPoints[i] ;

    boundingPoints.resize(num_points) ;
    double angle = 2.*M_PI/ (num_points) ;

    for (size_t i = 0 ; i < num_points ; i++)
    {
        boundingPoints[i] = new Point(getRadius()*cos((double)i*(angle)) + getCenter().getX(), getRadius()*sin((double)i*(angle)) + getCenter().getY());
    }
}

void Circle::sampleSurface(double linearDensity)
{
    
    if(!sampled)
    {
        sampleBoundingSurface(linearDensity*4.5) ;
        sampled = true ;
        size_t numberOfRings = static_cast<size_t>((double)getBoundingPoints().size()/(3. * M_PI )) ;
/*        size_t num_points = getBoundingPoints().size()*0.666666 ;

        double angle = 2.*M_PI/ (num_points) ;
        double offset = 0 ;

        size_t num_points_start = num_points ;*/

        std::vector<Point*> temp ;

        for (size_t i = 0 ; i< numberOfRings ; ++i)
        {
            double r = getRadius()*(1. - (double)(i + 1)/(numberOfRings+1)) ;

            Circle c(r, getCenter().getX(), getCenter().getY()) ;
            std::vector<Point> pts = c.getSamplingBoundingPoints(linearDensity*2.) ;

            if(pts.size() > 4)
            {
                for(size_t j = 0 ; j < pts.size() ; j++)
                    temp.push_back(new Point(pts[j].getX(), pts[j].getY())) ;
            }
            else
                break ;

//            for (size_t j = 0 ; j< num_points ; ++j)
//            {
//                temp.push_back(new Point(r*cos((double)(j+0.5*(i))*angle+offset) + getCenter().getX(), r*sin((double)(j+0.5*(i))*angle+offset) + getCenter().getY()));
//            }

//            num_points = (size_t)(/*std::max(*/(double)num_points_start*(r/getRadius())/*, (double)8)*/) ;

//            angle = 2.*M_PI/ (num_points) ;

//            offset = 0.5*(2.*M_PI/ (num_points) -  2.*M_PI/ (num_points*1.1)) ;

//            if(num_points < 5)
//                break ;
        }
        for(size_t i = 0 ; i < inPoints.size() ; i++)
            delete inPoints[i] ;

        inPoints.resize(temp.size() + 1) ;
        inPoints[0] = new Point(center) ;
        std::copy(temp.begin(), temp.end(),&inPoints[1]) ;
    }
}

bool Circle::in(const Point & v) const
{
    if(v.getX() < getCenter().getX()-getRadius()*1.01)
        return false ;
    if(v.getX() > getCenter().getX()+getRadius()*1.01)
        return false ;
    if(v.getY() < getCenter().getY()-getRadius()*1.01)
        return false ;
    if(v.getY() > getCenter().getY()+getRadius()*1.01)
        return false ;

    return squareDist2D(v, getCenter()) < sqradius ;
}

double Circle::getRadius() const
{
    return radius ;
}

double Circle::getSquareRadius() const
{
    return sqradius ;
}

double Circle::area() const
{
    return M_PI*sqradius ;
}


const std::vector<double> & LayeredCircle::getRadii() const
{
    return radiuses ;
}

LayeredCircle::LayeredCircle(std::vector<double> radii,double originX,double originY): Circle(radii[0], originX, originY)
{
    this->gType = LAYERED_CIRCLE ;
    std::sort(radii.begin(), radii.end()) ;
    Circle::setRadius(*radii.rbegin()) ;
    radiuses = radii ;
    center = Point(originX, originY) ;
}


LayeredCircle::LayeredCircle(std::vector<double> radii, const Point c) : Circle(radii[0], c)
{
    this->gType = LAYERED_CIRCLE ;
    std::sort(radii.begin(), radii.end()) ;
    Circle::setRadius(*radii.rbegin()) ;
    radiuses = radii ;
    center = c ;
}

LayeredCircle::LayeredCircle(double r, const Point center) : Circle(r, center)
{
    radiuses.push_back(r) ;
}

void LayeredCircle::sampleSurface(size_t num_points)
{

    std::vector<double> samplingRadiuses = radiuses;
    bool to_add = true ;
    std::vector<Point*> temp ;
    size_t numberOfRings = static_cast<size_t>((double)num_points/(2. * M_PI )) ;
    double meanDelta = getRadius()/(numberOfRings+1);

    while(to_add)
    {
        to_add = false ;
        std::vector<double> newRadii ;
        if(samplingRadiuses[0] > meanDelta*.25)
        {
            newRadii.push_back(samplingRadiuses[0]*.5) ;
            newRadii.push_back(samplingRadiuses[0]) ;
            to_add = true ;
        }
        for(size_t i = 1 ; i < samplingRadiuses.size() ; i++)
        {
            newRadii.push_back(samplingRadiuses[i]) ;
            if(samplingRadiuses[i]-samplingRadiuses[i-1] > meanDelta*.5)
            {
                newRadii.push_back(samplingRadiuses[i-1]*.5+samplingRadiuses[i]*.5) ;
                to_add = true ;
            }
        }
        std::sort(newRadii.begin() , newRadii.end()) ;

        samplingRadiuses = newRadii ;
    }


    num_points = std::max(num_points * 2.* getRadius()/meanDelta, 12.);
    size_t num_points_start = num_points ;
    sampleBoundingSurface(num_points) ;
    sampled = true ;

    double angle = 2.*M_PI/ (num_points) ;
    double offset = 0 ;

    for (int i = samplingRadiuses.size()-2 ; i >= 0 ; --i)
    {
        double r = samplingRadiuses[i] ;
        for (size_t j = 0 ; j< num_points ; ++j)
        {
            double randa= 0 ; //((2.*(double)rand()/(RAND_MAX+1.0))-1.)*0.2*(M_PI/num_points) ;
            double randr= 0 ; //(.2*r/(numberOfRings+1))*((double)rand()/RAND_MAX*2.-1.0) ;
            temp.push_back(new Point((r+randr)*cos((double)(j+0.5*(i))*angle+randa+offset) + getCenter().getX(), (r+randr)*sin((double)(j+0.5*(i))*angle+randa) + getCenter().getY()));
        }

        num_points = (size_t)((double)num_points_start*(r/getRadius())) ;

        angle = 2.*M_PI/ (num_points) ;

        offset = 0.5*(2.*M_PI/ (num_points) -  2.*M_PI/ (num_points*1.1)) ;
    }

    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    inPoints.resize(temp.size() + 1) ;
    inPoints[0] = new Point(center) ;
    std::copy(temp.begin(), temp.end(),&inPoints[1]) ;
}

void LayeredCircle::setRadius(double newr)
{
    double ratio = newr/(*radiuses.rbegin()) ;
    radiuses.pop_back() ;
    radiuses.push_back(newr) ;
    std::sort(radiuses.begin(), radiuses.end()) ;


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
    this->radius = newr ;
    this->sqradius = newr*newr ;
}

void LayeredCircle::addRadius(double newr)
{
    radiuses.push_back(newr) ;
    std::sort(radiuses.begin(), radiuses.end()) ;
    this->sqradius = (*radiuses.rbegin())*(*radiuses.rbegin()) ;
}

SegmentedLine::SegmentedLine(const std::valarray<Point *> & points) : NonConvexGeometry(0)
{
    gType = SEGMENTED_LINE ;

    if(points[points.size()/2])
        this->center = *points[points.size()/2] ;
    else
        this->center = Point() ;
    this->boundingPoints.resize(points.size()) ;
    for(size_t i = 0 ; i < points.size() ; i++)
        boundingPoints[i] = points[i] ;
}

SegmentedLine::SegmentedLine(const Point & p0, const Point & p1) : NonConvexGeometry(0)
{
    gType = SEGMENTED_LINE ;

    this->center = (p0+p1)*.5 ;

    this->boundingPoints.resize(2) ;
    boundingPoints[0] = new Point(p0) ;
    boundingPoints[1] = new Point(p1) ;
}


void SegmentedLine::computeCenter()
{
}

double SegmentedLine::getRadius() const
{
    //! \todo make it do something
    return 0 ;
}

Point * SegmentedLine::getHead() const
{
    return  boundingPoints[0];
}

Point * SegmentedLine::getTail() const
{
    return boundingPoints[this->boundingPoints.size()-1] ;
}

std::vector<Point> SegmentedLine::getSamplingBoundingPoints(double linearDensity) const
{
    std::vector<Point> ret ;

    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        ret.push_back(*boundingPoints[i]) ;

    return ret ;
}

void SegmentedLine::sampleBoundingSurface(double linearDensity)
{
}

void SegmentedLine::project(Point *p) const
{

    std::map<double, Point> projections ;

    for(size_t i = 0 ; i < boundingPoints.size()-1 ; i++)
    {
        Point proj = Segment(getBoundingPoint(i), getBoundingPoint(i+1)).project(*p) ;
        projections[squareDist2D(*p, proj)] = proj ;
    }

    p->getX() = projections.begin()->second.getX() ;
    p->getY() = projections.begin()->second.getY() ;
}

void SegmentedLine::sampleSurface(double linearDensity)
{

}

bool SegmentedLine::in(const Point & v) const
{
    return false ;
    for (size_t i = 0 ; i < boundingPoints.size() ; i++)
        if(*boundingPoints[i] == v)
            return true ;

    return false ;
}


Ellipse::Ellipse(Point center, Point a, Point b) : majorAxis(a), minorAxis(b)
{
    gType = ELLIPSE ;
    this->center = center ;
    this->computeCenter() ;
}

Ellipse::Ellipse(Circle c)
{
    gType = ELLIPSE ;
    this->center = c.getCenter() ;
    majorAxis = Point(c.getRadius(),0.) ;
    minorAxis = Point(0.,c.getRadius()) ;
    this->computeCenter() ;
}

Ellipse::Ellipse(Point center, Point a, double b) : majorAxis(a)
{
    gType = ELLIPSE ;
    double b_ = std::min(std::abs(b), 1.) ;
    this->minorAxis = Point(-a.getY()*b_, a.getX()*b_) ;
    this->center = center ;
    this->computeCenter() ;
}


Ellipse::Ellipse(Point center, double a, double b)
{
    gType = ELLIPSE ;
    this->center = center ;
    double dir_x = (double) rand() / (double) RAND_MAX ;
    double dir_y = (double) rand() / (double) RAND_MAX ;
//     double n = sqrt(dir_x*dir_x+dir_y*dir_y) ;
//     dir_x /=n ;
//     dir_y /=n ;
    double a_ = std::max(std::abs(a),std::abs(b)) ;
    double b_ = std::min(std::abs(a),std::abs(b)) ;
    this->majorAxis = Point(a_*dir_x,a_*dir_y) ;
    this->minorAxis = Point(-b_*dir_y,b_*dir_x) ;
    this->computeCenter() ;
}


const Point Ellipse::getFocus(bool dir) const
{
    int dirsign = 1 ;
    if(dir)
    {
        dirsign = - 1 ;
    }
    Point focus(center + majorAxis * (getExcentricity() * dirsign)) ;
    return focus ;
}

const std::pair<Point, Point> Ellipse::getBothFocus() const
{
    std::pair<Point, Point> focus(getFocus(true),getFocus(false)) ;
    return focus ;
}

void Ellipse::computeCenter()
{
       double h = (getMajorRadius()-getMinorRadius())*(getMajorRadius()-getMinorRadius())/((getMajorRadius()+getMinorRadius())*(getMajorRadius()+getMinorRadius())) ;
       perimeter = M_PI*(getMajorRadius()+getMinorRadius())*(1.+h*h/4.+h*h*h*h/64.+h*h*h*h*h*h/256.+25.*h*h*h*h*h*h*h*h/16384.)  ;
}


const Point Ellipse::getPointOnEllipse(double theta) const
{
    Point p(center) ;
    p += getMajorAxis() * std::cos(theta) ;
    p += getMinorAxis() * std::sin(theta) ;
    return p ;
}


Point Ellipse::project(const Point &p) const
{
    if((p-center).norm() < POINT_TOLERANCE*100.)
    {
        return p+getMajorAxis() ;
    } else {
        Point proj(p.getX(), p.getY()) ;
        this->project(&proj) ;
        return proj ;
    }
}

void Ellipse::project(Point * p) const
{
    Point test(p->getX(),p->getY()) ;
    if(test == center)
    {
//		p->getX() += getMinorAxis().getX() ;
//		p->getY() += getMinorAxis().getY() ;
        p->print() ;
        return ;
    }
    else
    {
        double alpha = majorAxis.angle() ;
        Point prot((*p-center).getX()*cos(-alpha)-(*p-center).getY()*sin(-alpha),
                   +(*p-center).getX()*sin(-alpha)+(*p-center).getY()*cos(-alpha)) ;

        Line majL(this->getCenter(), this->getMajorAxis()) ;
        Line minL(this->getCenter(), this->getMinorAxis()) ;

        if(minL.on(test))
        {
            // case point is colinear to minor axis
            Point A = center + getMinorAxis() ;
            Point B = center - getMinorAxis() ;
            Point pa(p->getX()-A.getX(),p->getY()-A.getY())  ;
            Point pb(p->getX()-B.getX(),p->getY()-B.getY())  ;
            if(pa.norm() < pb.norm())
            {
                p->getX() = A.getX() ;
                p->getY() = A.getY() ;
                return ;
            }
            else
            {
                p->getX() = B.getX() ;
                p->getY() = B.getY() ;
                return ;
            }
        }

        if(majL.on(test))
        {
            // case point is colinear to major axis
            Point C = center + getMajorAxis() ;
            Point D = center - getMajorAxis() ;
            Point pc(p->getX()-C.getX(),p->getY()-C.getY())  ;
            Point pd(p->getX()-D.getX(),p->getY()-D.getY())  ;
            if(pc.norm() < pd.norm())
            {
                p->getX() = C.getX() ;
                p->getY() = C.getY() ;
                return ;
            }
            else
            {
                p->getX() = D.getX() ;
                p->getY() = D.getY() ;
                return ;
            }
        }

        // else, we need brute force

        double tmin = 0. ;
        double tmax = M_PI * 2. ;
        int iter = 0 ;
        while(std::abs(tmin-tmax) > POINT_TOLERANCE && iter < 16)
        {
            iter++ ;
            double dtheta = (tmax - tmin) / 32. ;
            Vector allDistances(33) ;
            for(size_t i = 0 ; i < allDistances.size() ; i++)
                allDistances[i] = dist( test, getPointOnEllipse( tmin+dtheta*i )) ;

            double j = 0 ;
            for(size_t i = 0 ; i < allDistances.size() ; i++)
            {
                if(allDistances[i] == allDistances.min())
                {
                    j = (double) i ;
                    break ;
                }
            }
            tmin = tmin+dtheta*(j-1) ;
            tmax = tmin+dtheta*(j+1) ;
        }

        double found = (tmin+tmax)/2. ;

        p->getX() = getPointOnEllipse(found).getX() ;
        p->getY() = getPointOnEllipse(found).getY() ;

    }

}

std::vector<Point> Ellipse::getSamplingBoundingPoints(double linearDensity) const
{
    std::vector<Point> ret ;

    size_t num_points = 2*round((double) perimeter*(double) linearDensity) ;
    if(num_points < 4)
       return ret ;
    double dist = perimeter/num_points ;
    double theta = (double) std::rand()/RAND_MAX*2.*M_PI/num_points ;
    Vector allTheta(num_points*8) ; allTheta = 0. ;

    std::vector<Point> all ;
    for(size_t i = 0 ; i < num_points*8 ; i++)
    {
        all.push_back( getPointOnEllipse( theta + 2.*M_PI*(double) i/((double) num_points*8) ) ) ;
        allTheta[i] = theta+2.*M_PI*(double) i/((double) num_points*8) ;
        if(allTheta[i] > M_PI) { allTheta[i] -= 2.*M_PI ; }
    }

    ret.push_back( all[0] ) ;
    double finish = 1e9 ;
    size_t i = 0 ;
//    std::cout << all.size() << std::endl ;
    while(finish > dist*dist && i < all.size()-2)
    {
        size_t j = i+1 ;
        double coeff = 1. ;
        if(std::abs(std::cos(allTheta[j])) > 0.9)
            coeff = 0.8 ;
        while(squareDist2D(ret[ret.size()-1], all[j]) < dist*dist*coeff && j < all.size()-1) 
            j++ ;
        if(j != all.size()-1)
            ret.push_back( all[j] ) ;
        else
            break ;
        finish = squareDist2D(ret[0], ret[ret.size()-1])/(coeff*coeff) ;
        i = j ;
    }
    if(finish < dist*dist*0.5)
        ret.pop_back() ;
    
/*    while(ret.size() < num_points || finish > dist*dist)
    {
        Vector next(8) ;
        next = 0 ;
        for(size_t j = 0 ; j < next.size() ; j++)
            next[j] = squareDist2D( ret[ret.size()-1], getPointOnEllipse(theta+angle*(next.size()/2+j)/(next.size())) ) ;

         next -= dist*dist ;
         next = abs(next) ;
         for(size_t j = 0 ; j < next.size() ; j++)
         {
             if(next[j] == next.min())
             {
                 ret.push_back(getPointOnEllipse(theta+angle*(next.size()/2+j)/(next.size()))) ;
                 theta += angle*(next.size()/2+j)/(next.size()) ;
                 break ;
             }
         }

         finish = squareDist2D(ret[0], ret[ret.size()-1]) ;

    }*/


    return ret ;

/*    double multiplier = 2.- getMinorRadius() / getMajorRadius();
    double num_points = std::max(round(multiplier*getPerimeter()*linearDensity),8.) ;
    
    double angle = 2. * M_PI / num_points ;

    double thisangle = angle ;
    double lastangle = 0. ;
    double redfactor = 0.8 ; // factor for angle decrease < 1

    ret.push_back(getPointOnEllipse(0.)) ;


    Point thispoint = getPointOnEllipse(angle) ;
    Point lastpoint = getPointOnEllipse(lastangle) ;
    Point lastlastpoint = getPointOnEllipse(-angle) ;
    double criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;

    int n_iter = 0 ;



    bool accepted = true ;
    while (accepted)
    {
        n_iter = 0 ;
        while(((criteria > angle) || (criteria < angle * redfactor)) && n_iter < 20)
        {
            if(criteria > angle)
            {
                thisangle = lastangle + (thisangle - lastangle) * redfactor ;
            }
            else
            {
                n_iter = 20 ;
            } //thisangle = lastangle + (thisangle - lastangle) / redfactor ; }6290 6298 3.4343e-09 3.10041e-09

            thispoint = getPointOnEllipse(thisangle) ;
            criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;
            n_iter++ ;
        }

        if(thisangle < 2.*M_PI)
            ret.push_back(thispoint) ;
        else
            accepted = false ;

        lastlastpoint = lastpoint ;
        lastpoint = thispoint ;

        lastangle = thisangle ;
        thisangle += angle ;

        thispoint = getPointOnEllipse(thisangle) ;
    }

    if(ret.size() > 1)
    {

        double div = (ret[0]-ret[1]).norm() ;

        while(ret.size() > 1 && (ret[ret.size()-1]-ret[0]).norm() < 0.15*div)
        {
            ret.pop_back() ;
        }
    }

    return ret ;*/

}

Function Ellipse::getEllipseFormFunction() const
{
    double alpha = majorAxis.angle() ;
    Function x("x") ;
    Function y("y") ;
    Function x_((x-center.getX())*cos(alpha)-(y-center.getY())*sin(-alpha)) ;
    Function y_((x-center.getX())*sin(-alpha)+(y-center.getY())*cos(alpha)) ;
    double mm0 = (getMajorRadius()*getMajorRadius()) ;
    double mm1 = (getMinorRadius()*getMinorRadius()) ;
    Function f0 = x_*x_/mm0 ;
    Function f1 = y_*y_/mm1 ;
    Function f2 = f0+f1-1. ;
    return f2 ;
}

bool Ellipse::in(const Point &p) const
{
    Point q = toLocalCoordinates(p) ;
    double a = getMajorRadius() ;
    double b = getMinorRadius() ;
    return q.getX()*q.getX()/(a*a) + q.getY()*q.getY()/(b*b) - 1. < POINT_TOLERANCE ;
}

std::vector<Point> Ellipse::getSampleBoundingPointsOnArc(size_t num_points, double alpha, double beta) const
{
    std::vector<Point> ret ;

    double angle = (beta-alpha) * M_PI / num_points ;

    double thisangle = alpha + angle ;
    double lastangle = alpha ;
    double redfactor = 0.8 ; // factor for angle decrease < 1


    Point thispoint = getPointOnEllipse(angle) ;
    Point lastpoint = getPointOnEllipse(lastangle) ;
    Point lastlastpoint = getPointOnEllipse(alpha-angle) ;
    double criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;

    int n_iter = 0 ;
    std::default_random_engine generator;
    std::uniform_real_distribution< double > distribution(-angle*0.5, angle*0.5);

    for (size_t i = 0 ; i< num_points + 1 ; i++)
    {
        n_iter = 0 ;
        while(((criteria > angle) || (criteria < angle * redfactor)) && n_iter < 20)
        {
            if(criteria > angle)
            {
                thisangle = lastangle + (thisangle - lastangle) * redfactor ;
            }
            else
            {
                thisangle = lastangle + (thisangle - lastangle) / redfactor ;
            }
            thispoint = getPointOnEllipse(thisangle) ;
            criteria = acos (((lastpoint - lastlastpoint) * (thispoint - lastlastpoint) / ((lastpoint - lastlastpoint).norm() * (thispoint - lastlastpoint).norm()))) ;
            n_iter++ ;
        }
        thispoint = getPointOnEllipse(thisangle + distribution(generator)) ;

        ret.push_back(thispoint) ;

        lastlastpoint = lastpoint ;
        lastpoint = thispoint ;

        lastangle = thisangle ;
        thisangle += angle ;

        thispoint = getPointOnEllipse(thisangle) ;
    }

    double div = (ret[0]-ret[1]).norm() ;

    while((ret[ret.size()-1]-ret[0]).norm() < 0.15*div)
    {
        ret.pop_back() ;
    }

    return ret ;

}

void Ellipse::sampleBoundingSurface (double linearDensity)
{
    if( getMinorAxis().norm() < POINT_TOLERANCE || getMajorAxis().norm() < POINT_TOLERANCE )
        return ;

    if(linearDensity*getPerimeter() < 3)
        return ;

    std::vector<Point> bound = getSamplingBoundingPoints(linearDensity) ;

    for(size_t i = 0 ; i <getBoundingPoints().size() ; i++)
        delete boundingPoints[i];

    getBoundingPoints().resize(bound.size()) ;

    for (size_t i = 0 ; i < bound.size() ; i++)
    {
        getBoundingPoints()[i] = new Point(bound[i]) ;
    }

}

void Ellipse::sampleSurface (double linearDensity)
{
    size_t num_points = round(linearDensity*getPerimeter()) ;
    if( getMinorAxis().norm() < POINT_TOLERANCE || getMajorAxis().norm() < POINT_TOLERANCE )
        return ;

    if(num_points < 3)
        num_points = 4 ;

    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

    inPoints.resize(1) ;
    inPoints[0] = new Point(center) ;

    sampleBoundingSurface(linearDensity*1.5) ;
    sampled = true ;
    if(getBoundingPoints().size() < 8)
        return ;

    double dist = perimeter/getBoundingPoints().size() ;

    double a = getMajorRadius() ;
    double b = getMinorRadius() ;
    std::vector<Point> toadd ;
    std::default_random_engine generator;
    std::uniform_real_distribution< double > distribution(0, dist*0.1);
    while(b > dist*1.5)
    {
        b -= dist*1.1 ;
        a -= dist*1.1 ;
        Ellipse elln(center,getMajorAxis()*a/getMajorRadius(), getMinorAxis()*b/getMinorRadius()) ;

            std::vector<Point> pn = elln.getSamplingBoundingPoints(linearDensity) ;
            if(pn.size() < 4)
                break ;
            for(size_t j = 0 ; j < pn.size() ; j++)
            {
//                pn[j] += Point( distribution(generator), distribution(generator) ) ;
/*                bool alone = in(pn[j]) ;
                if(!alone)
                {
                    pn[j].print() ;
                    std::cout << "not inside" << std::endl ;
                }

                size_t k = 0 ;
                while(alone && k < getBoundingPoints().size())
                {
                    if((pn[j]-getBoundingPoint(k)).norm() < dist)
                    {
                        alone = false ;
                    }
                    k++ ;
                }


                k = 0 ;
                while(alone && k < toadd.size())
                {
                    if((pn[j]-toadd[k]).norm() < dist)
                    {
                        alone = false ;
                    }
                    k++ ;
                }

                if(alone)*/
                    toadd.push_back(pn[j]) ;
//                else
//                    rm++ ;
            }
    }


    if(b > dist*0.6)
    {
        for(double x = dist ; x < a-dist*0.6 ; x += dist)
        {
            toadd.push_back( getCenter() + getMajorAxis()*x/getMajorRadius() ) ;
            toadd.push_back( getCenter() - getMajorAxis()*x/getMajorRadius() ) ;
        }
    }


        PointArray inTemp(inPoints) ;

        inPoints.resize(inPoints.size()+toadd.size()) ;

        for(size_t j = 0 ; j < inTemp.size() ; j++)
            inPoints[j] = inTemp[j] ;
        for(size_t j = 0 ; j < toadd.size() ; j++)
        {
            inPoints[j+inTemp.size()] = new Point(toadd[j]) ;
        }

/*    
    size_t ring = (num_points*3 / 5) - 1 ;

    if(getMinorRadius() / getMajorRadius() < 0.5)
        ring-- ;

    if(ring > 1)
    {
        std::vector<double> newa ;
        std::vector<double> newb ;
        std::vector<size_t> newn ;

        double minn = (double) num_points ;
        size_t inc = 1 ;
        if(ring > 4)
            inc++ ;

        for(size_t i = 0 ; i < ring ; i++)
        {
            newa.push_back(getMajorRadius() * (ring - i) / (ring + 1) ) ;
            newb.push_back(getMinorRadius() * (ring - i) / (ring + 1) ) ;
            newn.push_back(  minn + (size_t) ((double) num_points) * ((double) (ring-i)*(ring-i) / (double) ((ring+1)*(ring+1)))) ;
        }
        std::vector<Point> toadd ;
        int rm = 0 ;
        std::default_random_engine generator;
        std::uniform_real_distribution< double > distribution(0, dist*0.1);

        for(size_t i = 0 ; i < newa.size() ; i+=inc)
        {
            Ellipse elln(center,getMajorAxis()*newa[i]/getMajorRadius(),getMinorAxis()*newb[i]/getMinorRadius()) ;

//			int factor = 2 + i/2 ;
            size_t nell = newn[i] ;
            if(nell < minn) {
                nell = minn ;
            }

            std::vector<Point> pn = elln.getSamplingBoundingPoints(linearDensity) ;
            for(size_t j = 0 ; j < pn.size() ; j++)
            {
                pn[j] += Point( distribution(generator), distribution(generator) ) ;
                bool alone = in(pn[j]) ;
                if(!alone)
                {
                    pn[j].print() ;
                    std::cout << "not inside" << std::endl ;
                }

                size_t k = 0 ;
                while(alone && k < getBoundingPoints().size())
                {
                    if((pn[j]-getBoundingPoint(k)).norm() < dist)
                    {
                        alone = false ;
                    }
                    k++ ;
                }


                k = 0 ;
                while(alone && k < toadd.size())
                {
                    if((pn[j]-toadd[k]).norm() < dist)
                    {
                        alone = false ;
                    }
                    k++ ;
                }

                if(alone)
                    toadd.push_back(pn[j]) ;
                else
                    rm++ ;
            }
        }
        PointArray inTemp(inPoints) ;

        inPoints.resize(inPoints.size()+toadd.size()) ;

        for(size_t j = 0 ; j < inTemp.size() ; j++)
            inPoints[j] = inTemp[j] ;
        for(size_t j = 0 ; j < toadd.size() ; j++)
        {
            inPoints[j+inTemp.size()] = new Point(toadd[j]) ;
        }
    }
    else
    {
        if(ring == 1)
        {
            Ellipse elln(center, getMajorAxis() * 0.5, getMinorAxis() * 0.6666666) ;
            std::vector<Point> pn = elln.getSamplingBoundingPoints(linearDensityma) ;
            PointArray inTemp(inPoints) ;

            inPoints.resize(inPoints.size()+pn.size()) ;

            for(size_t j = 0 ; j < inTemp.size() ; j++)
                inPoints[j] = inTemp[j] ;
            for(size_t j = 0 ; j < pn.size() ; j++)
            {
                inPoints[j+inTemp.size()] = new Point(pn[j]) ;
            }
        }
    }*/
}

Ellipse Ellipse::getEllipseInLocalCoordinates() const
{
    Point c(center) ;
    Point a(getMajorRadius(), 0.) ;
    Point b(0.,getMinorRadius()) ;

    return Ellipse(c,a,b) ;

}

Point Ellipse::toLocalCoordinates(const Point & p) const
{
    double alpha = majorAxis.angle() ;
    Point prot((p-center).getX()*cos(-alpha)-(p-center).getY()*sin(-alpha),
               +(p-center).getX()*sin(-alpha)+(p-center).getY()*cos(-alpha)) ;
    return prot ;
}




std::vector<Point> Ellipse::getBoundingBox() const
{
    std::vector<Point> bbox(4) ;

    Point A = center + majorAxis + minorAxis ;
    Point B = center + majorAxis - minorAxis ;
    Point C = center - majorAxis - minorAxis ;
    Point D = center - majorAxis + minorAxis ;

    double minx = std::min(A.getX(), std::min(B.getX(), std::min(C.getX(), D.getX()))) ;
    double maxx = std::max(A.getX(), std::max(B.getX(), std::max(C.getX(), D.getX()))) ;
    double miny = std::min(A.getY(), std::min(B.getY(), std::min(C.getY(), D.getY()))) ;
    double maxy = std::max(A.getY(), std::max(B.getY(), std::max(C.getY(), D.getY()))) ;

    bbox[0] = Point(minx, maxy) ;
    bbox[1] = Point(maxx, maxy) ;
    bbox[2] = Point(maxx, miny) ;
    bbox[3] = Point(minx, miny) ;

    return  bbox ;



    double a = getMajorRadius() ;
    bbox[0] = center + Point(a,a) ;
    bbox[1] = center + Point(a,-a) ;
    bbox[2] = center + Point(-a,-a) ;
    bbox[3] = center + Point(-a,a) ;

    return bbox ;

//	double a = getMajorRadius() ;
    double b = getMinorRadius() ;
    double alpha = getMajorAxis().angle() ;

    double thetax1 = std::acos(std::sqrt((a*a*std::cos(alpha)*std::cos(alpha)) / (a*a*std::cos(alpha)*std::cos(alpha) + b*b*std::sin(alpha)*std::sin(alpha)))) ;
    double thetax2 = std::acos(- std::sqrt((a*a*std::cos(alpha)*std::cos(alpha)) / (a*a*std::cos(alpha)*std::cos(alpha) + b*b*std::sin(alpha)*std::sin(alpha)))) ;
    double thetay1 = std::asin(std::sqrt((b*b*std::cos(alpha)*std::cos(alpha)) / (b*b*std::cos(alpha)*std::cos(alpha) + a*a*std::sin(alpha)*std::sin(alpha)))) ;
    double thetay2 = std::asin(- std::sqrt((b*b*std::cos(alpha)*std::cos(alpha)) / (b*b*std::cos(alpha)*std::cos(alpha) + a*a*std::sin(alpha)*std::sin(alpha)))) ;

    double x1 = a*std::cos(alpha)*std::cos(thetax1) - b*std::sin(alpha)*std::sin(thetax1) ;
    double x2 = a*std::cos(alpha)*std::cos(thetax2) - b*std::sin(alpha)*std::sin(thetax2) ;

    double xmin = std::min(x1,x2) ;
    double xmax = std::max(x1,x2) ;

    double y1 = a*std::sin(alpha)*std::cos(thetay1) + b*std::cos(alpha)*std::sin(thetay1) ;
    double y2 = a*std::sin(alpha)*std::cos(thetay2) + b*std::cos(alpha)*std::sin(thetay2) ;

    double ymin = std::min(y1,y2) ;
    double ymax = std::max(y1,y2) ;

    bbox[0] = center + Point(xmax,ymax) ;
    bbox[1] = center + Point(xmax,ymin) ;
    bbox[2] = center + Point(xmin,ymin) ;
    bbox[3] = center + Point(xmin,ymax) ;

    return bbox ;

}


void Polygon::computeCenter()
{

    for(size_t i = 0 ; i < originalPoints.size() ; i++)
    {
        center += originalPoints[i] ;
    }
    center /= originalPoints.size() ;
}


Polygon::Polygon(const std::valarray<Point *> & points) : NonConvexGeometry(0), originalPoints(points.size())
{
    gType = POLYGON ;

    boundingPoints.resize(points.size()) ;
    for(size_t i = 0 ; i < points.size() ; i++)
    {
        boundingPoints[i] = points[i] ;
        originalPoints[i] = *(points[i]) ;
    }

    computeCenter();
}

Polygon::Polygon(const std::valarray<Point> & points) : NonConvexGeometry(0), originalPoints(points.size())
{
    gType = POLYGON ;

//    boundingPoints.resize(points.size()) ;
    for(size_t i = 0 ; i < points.size() ; i++)
    {
//        boundingPoints[i] = new Point(points[i]) ;
        originalPoints[i] = (points[i]) ;
    }

    computeCenter();
}

Polygon::~Polygon() { }

void Polygon::sampleBoundingSurface(double linearDensity)
{
    std::vector<Point> newPoints = getSamplingBoundingPoints(linearDensity) ;

    boundingPoints.resize(newPoints.size());

    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        boundingPoints[i] = new Point(newPoints[i]) ;

}

std::vector<Point> Polygon::getSamplingBoundingPoints(double linearDensity) const
{
    std::vector<Point> ret ;

    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        double fraction = dist(originalPoints[i], originalPoints[inext]) ;
        int numPointsOnSegment = std::max(round(fraction*linearDensity)+1, (fraction > linearDensity?3.:2.) ) ;

/*        if( numPointsOnSegment == 2 && dist(originalPoints[i], originalPoints[inext]) > perimeter*0.25/(num_points) )
            numPointsOnSegment++ ;*/

//        if(numPointsOnSegment > 2 || dist(originalPoints[i], originalPoints[inext]) > linearDensity)
//            ret.push_back(originalPoints[i]);        
        for(int j = 0 ; j < numPointsOnSegment-1 ; j++)
        {
            double f = ((double)j/(numPointsOnSegment-1)) ;
            ret.push_back(originalPoints[i]*(1.-f) + (originalPoints[inext]*f) );
        }
    }
    return ret ;
}

double Polygon::getPerimeter() const
{
    double perimeter = 0 ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        perimeter += dist(originalPoints[i], originalPoints[inext]) ;
    }
    return perimeter ;
}

std::vector<Polygon> Polygon::getInscribedPolygons( double delta ) const
{
    // first, get tentative segments inside p
    std::valarray<Point> realvertex = getOriginalPoints() ;
    std::vector<Point> vertex ;
    for(size_t i = 0 ; i < realvertex.size() ; i++)
    {
       size_t i_next = (i+1)%realvertex.size() ;
       size_t i_prev = (i==0 ? realvertex.size()-1 : i-1) ;
       if(dist( realvertex[i], realvertex[i_next]) > delta*0.5 && dist( realvertex[i], realvertex[i_prev]) > delta*0.5)
           vertex.push_back(realvertex[i]) ;
       else if(dist(realvertex[i], realvertex[i_next]) < delta*0.5 && dist( realvertex[i], realvertex[i_prev]) > delta*0.5 )
	{
           vertex.push_back( (realvertex[i]+realvertex[i_next])*0.5 ) ;
	}
    }

    std::vector<Segment> segments ;
    for(size_t i = 0 ; i < vertex.size() ; i++)
    {

       size_t i_next = (i+1)%vertex.size() ;
       size_t i_next_next = (i+2)%vertex.size() ;
       size_t i_prev = (i==0 ? vertex.size()-1 : i-1) ;
       Segment edge_prev( vertex[i_prev], vertex[i] ) ;
       Segment edge( vertex[i], vertex[i_next] ) ;
       Segment edge_next( vertex[i_next], vertex[i_next_next] ) ;

       Point n = edge.normal() ;
       if(n.norm() < POINT_TOLERANCE)
           continue ; // original edge is too small
       n *= 0.75*delta/n.norm() ;       
       if(!in(edge.midPoint() + n*0.15))
       {
          n *= -1 ;
          if(!in(edge.midPoint() + n*0.15))
             continue ; // original polygon is too thin
       }
       Line tentative( edge.midPoint()+n, vertex[i_next]-vertex[i] ) ;

       n = edge_prev.normal() ;
       if(n.norm() < POINT_TOLERANCE)
           continue ; // original edge is too small
       n *= 0.75*delta/n.norm() ;       
       if(!in(edge_prev.midPoint() + n*0.15))
       {
          n *= -1 ;
          if(!in(edge_prev.midPoint() + n*0.15))
             continue ; // original polygon is too thin
       }
       Line tentative_prev( edge_prev.midPoint()+n, vertex[i]-vertex[i_prev] ) ;

       n = edge_next.normal() ;
       if(n.norm() < POINT_TOLERANCE)
           continue ; // original edge is too small
       n *= 0.75*delta/n.norm() ;       
       if(!in(edge_next.midPoint() + n*0.15))
       {
          n *= -1 ;
          if(!in(edge_next.midPoint() + n*0.15))
             continue ; // original polygon is too thin
       }
       Line tentative_next( edge_next.midPoint()+n, vertex[i_next_next]-vertex[i_next] ) ;

/*

       Line bisectrix_prev( vertex[i], vertex[i_prev]+vertex[i_next]-(vertex[i]*2.) ) ;
       Line bisectrix_next( vertex[i_next], vertex[i]+vertex[i_next_next]-(vertex[i_next]*2.) ) ;*/

       if( (!tentative.intersects( tentative_prev )) || (!tentative.intersects( tentative_next ))) 
           continue ; // unable to find points for next tentative edge

       Point first = tentative.intersection( tentative_prev ) ;
       Point second = tentative.intersection( tentative_next ) ;
       if(!in( first ) || !in( second) )
       { 
           Segment s(first, second) ;
           if(in(first) && s.intersects(this))
           {
               Point test = s.intersection(this)[0]-first ;
               second = first + test*(test.norm()-delta)/test.norm() ;
           }
           else if(in(second) && s.intersects(this))
           {
               Point test = s.intersection(this)[0]-second ;
               first = second + test*(test.norm()-delta)/test.norm() ;
           }
       }

       if(!in( first ) || !in( second) )
           continue ;

       segments.push_back( Segment( first, second ) ) ;
    }

    if(segments.size() < 3 && vertex.size() > 2)
    {
        Point c = vertex[0] ;
        for(size_t i = 1 ; i < vertex.size() ; i++)
            c += vertex[i] ;
        c /= vertex.size() ;

       std::valarray<Point> nextvertex(vertex.size()) ;
       bool valid = true ;
       for(size_t i = 0 ; i < nextvertex.size() ; i++)
       {
           Point r = vertex[i]-c ;
           if(r.norm() > delta)
               nextvertex[i] = c + r*(r.norm()-delta)/r.norm() ;
           else
               valid = false ;
       }
       std::vector<Polygon> poly ;
       if(valid)
           poly.push_back( Polygon(nextvertex) ) ;
       return poly ;

    }


    // second, create list of connectivity
    std::vector<std::pair< Point, std::vector<int> > > nodes ;
    std::vector<int> dummy ;
    for(size_t i = 0 ; i < segments.size() ; i++)
    {
        if(nodes.empty())
        {
            nodes.push_back( std::make_pair( segments[i].first(), dummy ) ) ;
            nodes.push_back( std::make_pair( segments[i].second(), dummy ) ) ;
            nodes[0].second.push_back(1) ;
            nodes[1].second.push_back(0) ;
            continue ;
        }

        int found_first_index = -1 ;
        int found_second_index = -1 ;
        for(size_t j = 0 ; j < nodes.size() ; j++)
        {
            if(dist(segments[i].first(), nodes[j].first) < POINT_TOLERANCE)
            {
                found_first_index = j ;
                break ;
            }            
        }
        for(size_t j = 0 ; j < nodes.size() ; j++)
        {
            if(dist(segments[i].second(), nodes[j].first) < POINT_TOLERANCE)
            {
                found_second_index = j ;
                break ;
            }            
        }
        if(found_first_index < 0)
        {
            found_first_index = nodes.size() ;
            nodes.push_back( std::make_pair( segments[i].first(), dummy ) ) ;
        }
        if(found_second_index < 0)
        {
            found_second_index = nodes.size() ;
            nodes.push_back( std::make_pair( segments[i].second(), dummy ) ) ;
        }
        nodes[found_first_index].second.push_back(found_second_index) ;
        nodes[found_second_index].second.push_back(found_first_index) ;
    }

    // third, merge points too close
    std::valarray<bool> done( nodes.size() ) ;
    std::valarray<bool> kept( nodes.size() ) ;
    done = false ;
    kept = false ;
    for(size_t i = 0 ; i < nodes.size() ; i++)
    {
        if(done[i])
            continue ;
        std::vector<int> nodes_close ;
        for(size_t j = 0 ; j < nodes.size() ; j++)
        {
            if(i == j)
                continue ;

            if(dist( nodes[i].first, nodes[j].first ) < 0.5*delta)
                nodes_close.push_back(j) ;            
        }

        Point c = nodes[i].first ;
        for(size_t j = 0 ; j < nodes_close.size() ; j++)
            c += nodes[ nodes_close[j] ].first ;
        c /= 1+nodes_close.size() ; 
        nodes[i].first = c ;
        for(size_t j = 0 ; j < nodes_close.size() ; j++)
        {
            done[ nodes_close[j] ] = true ;
            kept[ nodes_close[j] ] = false ;
            for(size_t k = 0 ; k < nodes[ nodes_close[j] ].second.size() ; k++)
                nodes[i].second.push_back( nodes[ nodes_close[j] ].second[ k ] ) ;
            for(size_t k = 0 ; k < nodes.size() ; k++)
            {
                for(size_t l = 0 ; l < nodes[k].second.size() ; l++)
                {
                    if(nodes[k].second[l] == nodes_close[j])
                        nodes[k].second[l] = i ;
                }
            }
        }
        std::sort( nodes[i].second.begin(), nodes[i].second.end()) ;
        nodes[i].second.erase( std::unique( nodes[i].second.begin(), nodes[i].second.end()), nodes[i].second.end()) ;

        done[i] = true ;
        kept[i] = true ;
    }
    int keptNodes = 0 ;
    for(size_t i = 0 ; i < kept.size() ; i++) { keptNodes += (int) kept[i] ; }

    // fourth, recreate tentative segments
    segments.clear() ;
    done = false ;
    for(size_t i = 0 ; i < nodes.size() ; i++)
    {
        if(!kept[i])
            continue ; // node has been merged with another

        for(size_t j = 0 ; j < nodes[i].second.size() ; j++)
        {
            if((size_t) nodes[i].second[j] == i)
                continue ; // node connects with itself

            if(done[ nodes[i].second[j] ])
                continue ; // segment has already been added

            segments.push_back( Segment( nodes[i].first, nodes[ nodes[i].second[j] ].first ) ) ;
        }
        done[i] = true ;
    }

    // fifth, remove points too close to other edge 
    for(size_t i = 0 ; i < segments.size() ; i++)
    {
        for(size_t j = i+1 ; j < segments.size() ; j++)
        {
            if( dist( segments[i].first(), segments[j].second() ) < POINT_TOLERANCE || 
                dist( segments[i].first(), segments[j].first() ) < POINT_TOLERANCE ||
                dist( segments[i].second(), segments[j].first() ) < POINT_TOLERANCE ||
                dist( segments[i].second(), segments[j].second() ) < POINT_TOLERANCE)
                continue ; // segments are connected
            
            Point p = segments[j].project( segments[i].first() ) ;
            if(dist(p, segments[i].first()) < 0.5*delta) // first point is too close to next segment
            {
                // merge points
                for(size_t k = 0 ; k < segments.size() ; k++)
                {
                    if(k == i)
                       continue ;
                    if( dist( segments[i].first(), segments[k].first() ) < POINT_TOLERANCE)
                        segments[k].setFirst( p ) ;
                    if( dist( segments[i].first(), segments[k].second() ) < POINT_TOLERANCE)
                        segments[k].setSecond( p ) ;
                }
                segments[i].setFirst(p) ;
                // split segment
                Point s = segments[j].second() ;
                segments[j].setSecond( p ) ;
                segments.push_back( Segment( p, s ) ) ;
            }


            p = segments[j].project( segments[i].second() ) ;
            if(dist(p, segments[i].second()) < 0.25*delta) // second point is too close to next segment
            {
                for(size_t k = 0 ; k < segments.size() ; k++)
                {
                    if(k == i)
                       continue ;
                    if( dist( segments[i].second(), segments[k].first() ) < POINT_TOLERANCE)
                        segments[k].setFirst( p ) ;
                    if( dist( segments[i].second(), segments[k].second() ) < POINT_TOLERANCE)
                        segments[k].setSecond( p ) ;
                }
                segments[i].setSecond(p) ;
                // split segment
                Point s = segments[j].second() ;
                segments[j].setSecond( p ) ;
                segments.push_back( Segment( p, s ) ) ;
            }
        }
    }

    // last: reconstruct ordered polygons
    done.resize(segments.size()) ;
    done = false ;
    std::vector<Polygon> nextPolygons ;
    for(size_t i = 0 ; i < segments.size() ; i++)
    {
        if(done[i])
            continue ;

        bool close = false ;
        bool allchecked = false ;
        std::vector<Point> pts ;
        size_t j = i+1 ;
        pts.push_back( segments[i].first() ) ;
        pts.push_back( segments[i].second() ) ;
        done[i] = true ;
        size_t tries = 0 ;
        while(tries < segments.size() && !(allchecked || close) )
        {
             if(j == segments.size()) { j = 0 ; }
             if(j == i || done[j]) // move on
             {
                 j++ ;
                 tries++ ;
                 continue ;
             }
             Point last = pts[ pts.size()-1 ] ;
             if(dist( last, segments[j].second() ) < POINT_TOLERANCE ) // second point is last point added, so add first
             {
                 pts.push_back( segments[j].first()) ;
                 tries = 0 ;
                 done[j] = true ;
             }
             else if(dist( last, segments[j].first() ) < POINT_TOLERANCE ) // vice-versa
             {
                 pts.push_back( segments[j].second()) ;
                 tries = 0 ;
                 done[j] = true ;
             }

             j++ ;
             tries++ ;
             close = (dist(pts[0], pts[ pts.size()-1 ]) < POINT_TOLERANCE) ;
             allchecked = true ;
             for(size_t k = 0 ; allchecked && k < done.size() ; k++)
                 allchecked &= done[k] ;
        }

        if(pts.size() > 2)
        {
            std::valarray<Point> nextPoly( pts.size() - (int) close ) ;
            for(size_t k = 0 ; k < pts.size()-(int) close ; k++)
                nextPoly[k] = pts[k] ;
            nextPolygons.push_back( Polygon(nextPoly) ) ;
        }
    }

    return nextPolygons ;
}


void Polygon::sampleSurface(double linearDensity)
{
    for(size_t i = 0 ; i < boundingPoints.size() ; i++)
        delete boundingPoints[i] ;
    for(size_t i = 0 ; i < inPoints.size() ; i++)
        delete inPoints[i] ;

//    num_points *= 2 ;
//    double real_num = num_points*2. ;//std::sqrt((getRadius()*2.*M_PI/getPerimeter())) ; //*M_PI*(getRadius())/std::sqrt(area()) ;

    sampleBoundingSurface( linearDensity*3.);

    std::vector<Polygon> clusters ;
    clusters.push_back( Polygon(originalPoints) ) ;
    double perimeter = clusters[0].getPerimeter() ;
    double delta = 2.*perimeter/getBoundingPoints().size() ;
    std::vector<Point> newPoints ;
//    int t = 0 ; // for debugging purpose

    while(clusters.size() > 0)
    {
        std::vector<Polygon> nextCluster ;
        std::vector<Point> nextPoints ;
        for(size_t i = 0 ; i < clusters.size() ; i++)
        {
            std::vector<Polygon> inscribed = clusters[i].getInscribedPolygons( delta ) ;
            if( inscribed.size() == 0)// && clusters[i].getOriginalPoints().size() > 2)
                newPoints.push_back( clusters[i].getCenter() ) ;
            for(size_t j = 0 ; j < inscribed.size() ; j++)
            {
                std::vector<Point> next = inscribed[j].getSamplingBoundingPoints(linearDensity*2.) ;
                if(next.size() > 4)
                    nextPoints.insert( nextPoints.end(), next.begin(), next.end()) ;
//                else
//                    nextPoints.push_back( inscribed[j].getCenter() ) ;
//                std::valarray<Point> pts = inscribed[j].getOriginalPoints() ;
/*                for(size_t k = 0 ; k < pts.size() ; k++)
                {
                    bool alone = true ;
                    for(size_t l = 0 ; alone && l < newPoints.size() ; l++)
                       alone &= dist( pts[k], newPoints[l] ) > delta*0.5 ;
                    if(alone)
                       nextPoints.push_back( pts[k] ) ;
                    
                    if(k==1 && pts.size() == 2)
                        break ; // polygon is a line, do not repeat mesh points

                    size_t k_next = (k+1)%pts.size() ;
                    int numPointsOnSegment = ceil(dist(pts[k], pts[k_next])/delta) ;

                    for(int n = 0 ; n < numPointsOnSegment-1 ; n++)
                    {
                       alone = true ;
                       Point tentative = pts[k]*(double)n/(numPointsOnSegment-1) + pts[k_next]*(1.-(double)n/(numPointsOnSegment-1)) ;
                       for(size_t l = 0 ; alone && l < newPoints.size() ; l++)
                           alone &= dist( tentative, newPoints[l] ) > delta*0.5 ;
                       if(alone)
                           newPoints.push_back(tentative);
                    }
                }*/
                nextCluster.push_back( inscribed[j] ) ;
            }
        }
        for(size_t i = 0 ; i < nextPoints.size() ; i++)
             newPoints.push_back( nextPoints[i] ) ;

        clusters.clear() ;
        for(size_t i = 0 ; i < nextCluster.size() ; i++)
            clusters.push_back( nextCluster[i] ) ;

    }

    inPoints.resize(newPoints.size());
    for(size_t i = 0 ; i < newPoints.size() ; i++ )
        inPoints[i] = new Point(newPoints[i]) ;
}

Point Polygon::getOrientation() const
{
    Polygon convex = getConvexPolygon() ;
    std::valarray<Point> pts = convex.getOriginalPoints() ;

    size_t k = 0 ;
    double minwidth = 0. ;

    for(size_t i = 0 ; i < pts.size() ; i++)
    {
       size_t inext = (i+1)%pts.size() ;
       Line edge( pts[i], pts[inext]-pts[i] ) ;
       double minx = 0 ;
       double maxx = 0 ;
       double miny = 0 ;
       double maxy = 0 ;
       for(size_t j = 0 ; j < pts.size() ; j++)
       {
           Point p = edge.projection( pts[j] ) ;
           double x = dist( p, pts[i] ) ;
           double y = dist( p, pts[j] ) ;
           if(x < minx)
               minx = x ;
           if(x > maxx)
               maxx = x ;
           if(y < miny)
               miny = y ;
           if(y > maxy)
               maxy = y ;
       }
//       double height = maxx-minx ;
       double width = maxy-miny ;
       if( k == 0 && i == 0 )
           minwidth = width ;
       else
       {
           if(width < minwidth)
           {
               k = i ;
               minwidth = width ;
           }
       }

    }
    
    Point dir = pts[(k+1)%pts.size()] - pts[k] ;
    if(dir.norm() < POINT_TOLERANCE)
       dir.setX(1) ;
    else
       dir /= dir.norm() ;

    return dir ;
}

double Polygon::getAspectRatio() const
{
    Point dir = getOrientation() ;
    Line axis(getCenter(), dir) ;
    Point norm( dir.getY(), -dir.getX() ) ;
    double minx = 0 ;
    double maxx = 0 ;
    double miny = 0 ;
    double maxy = 0 ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++)
    {
        Point p = axis.projection( originalPoints[i] ) ;
        double x = (p-getCenter())*dir ;
        double y = (p-originalPoints[i])*norm ;
        if(i == 0)
        {
            minx = x ;
            maxx = x ;
            miny = y ;
            maxy = y ;
        }
        else
        {
           if(x < minx)
               minx = x ;
           if(x > maxx)
               maxx = x ;
           if(y < miny)
               miny = y ;
           if(y > maxy)
               maxy = y ;
        }
    }
    return (maxy-miny)/(maxx-minx) ;
}

Polygon Polygon::getConvexPolygon() const
{
    std::vector<Point> pts ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++)
    {
        size_t inext = (i+1)%originalPoints.size() ;
        size_t iprev = ( i == 0 ? originalPoints.size()-1 : i-1 ) ;
        Segment s( originalPoints[iprev], originalPoints[inext] ) ;
        Point test = originalPoints[i]-s.midPoint() ;
        if(test.norm() > POINT_TOLERANCE)
           test *= 1e-4 ;
        test = originalPoints[i]-test ;
        if(in(test))
           pts.push_back( originalPoints[i] ) ;
    }
    std::valarray<Point> next(pts.size()) ;
    for(size_t i = 0 ; i < pts.size() ; i++)
       next[i] = pts[i] ;
    return Polygon(next) ;
}


bool Polygon::in(const Point & v) const
{
    Point out = center + Point(0, getRadius()*2.) ;
    Segment s(out, v+Point(POINT_TOLERANCE,-POINT_TOLERANCE)) ;
    Segment s0(out, v+Point(-POINT_TOLERANCE,-POINT_TOLERANCE)) ;
    Segment s1(out, v+Point(-POINT_TOLERANCE,POINT_TOLERANCE)) ;
    Segment s2(out, v+Point(POINT_TOLERANCE,POINT_TOLERANCE)) ;

    int interCount = 0 ;
    int interheadcount = 0 ;
    int interCount0 = 0 ;
    int interheadcount0 = 0 ;
    int interCount1 = 0 ;
    int interheadcount1 = 0 ;
    int interCount2 = 0 ;
    int interheadcount2 = 0 ;

    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        Segment seg(originalPoints[i], originalPoints[inext]) ;
        bool inter = seg.intersects(s);
        interCount += inter ;
        if(inter)
        {
            Point intersec = seg.intersection(s) ;
            if(intersec == seg.first() || intersec == seg.second())
                interheadcount++ ;
        }
        inter = seg.intersects(s0);
        interCount0 += inter ;
        if(inter)
        {
            Point intersec = seg.intersection(s0) ;
            if(intersec == seg.first() || intersec == seg.second())
                interheadcount0++ ;
        }
        inter = seg.intersects(s1);
        interCount1 += inter ;
        if(inter)
        {
            Point intersec = seg.intersection(s1) ;
            if(intersec == seg.first() || intersec == seg.second())
                interheadcount1++ ;
        }
        inter = seg.intersects(s2);
        interCount2 += inter ;
        if(inter)
        {
            Point intersec = seg.intersection(s2) ;
            if(intersec == seg.first() || intersec == seg.second())
                interheadcount2++ ;
        }

    }
    return (interCount-interheadcount/2)%2 && (interCount0-interheadcount0/2)%2 && (interCount1-interheadcount1/2)%2 && (interCount2-interheadcount2/2)%2;
}

double Polygon::area() const
{
    std::vector<Point> bb = getBoundingBox() ;
    double maxx = bb[1].getX() ;
    double minx = bb[0].getX() ;
    double miny = bb[1].getY() ;
    double maxy = bb[2].getY() ;
    std::default_random_engine generator;
    std::uniform_real_distribution< double > distributionx(minx, maxx);
    std::uniform_real_distribution< double > distributiony(miny, maxy);

    double incount = 0 ;
    for(size_t i = 0 ; i < 1024  ; i++)
    {
        Point test(distributionx(generator), distributiony(generator)) ;
        if(in(test))
            incount++ ;
    }
    return (maxx-minx)*(maxy-miny)*(incount/1024.) ;

}

double Polygon::volume() const {
    return 0 ;
}

void Polygon::project(Point * init) const
{
    std::map<double, Point> potential ;
    for(size_t i = 0 ; i < originalPoints.size() ; i++ )
    {
        int inext = (i+1)%originalPoints.size() ;
        Segment s(originalPoints[i], originalPoints[inext]);
        Point p = s.project(*init);
        potential[squareDist2D(p, *init)] = p ;
    }
    init->set(potential.begin()->second) ;
}

double Polygon::getRadius() const
{
    double r = squareDist2D(center, originalPoints[0]) ;
    for(const auto & p : originalPoints)
        r = std::max(r, squareDist2D(p, center)) ;

    return sqrt(r) ;
}

SpaceDimensionality Polygon::spaceDimensions() const
{
    return SPACE_TWO_DIMENSIONAL ;
}

std::vector<Point> Polygon::getBoundingBox() const
{
    double maxx = originalPoints[0].getX() ;
    double minx = originalPoints[0].getX() ;
    double maxy = originalPoints[0].getY() ;
    double miny = originalPoints[0].getY() ;
    for(size_t i = 1 ; i < originalPoints.size() ; i++ )
    {
        maxx = std::max(maxx,originalPoints[i].getX()) ;
        minx = std::min(minx,originalPoints[i].getX()) ;
        maxy = std::max(maxy,originalPoints[i].getY()) ;
        miny = std::min(miny,originalPoints[i].getY()) ;
    }

    return {Point(minx, miny), Point(maxx, miny), Point(maxx, maxy), Point(minx,maxy)} ;

}

void Polygon::setCenter(const Point & newCenter)
{
    Point delta = newCenter-getCenter();
    Geometry::setCenter(newCenter);
    for(size_t i = 0 ; i < originalPoints.size() ; i++)
        originalPoints[i] += delta ;
}

}
