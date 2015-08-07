
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
// Author: Ruzena Chamrova <ruzena.chamrova@epfl.ch>, (C) 2007
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2009-2011
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_base.h"
//#include "../utilities/xml.h"
#include <limits>
#include <iomanip>

#include "../mesher/delaunay.h"
#include "space_time_geometry_2D.h"
#include "../polynomial/vm_function_base.h"

namespace Amie
{



#ifdef HAVE_SSE3
Point::Point(const Point & p) : vecxy(p.vecxy), veczt(p.veczt), id(p.getId())
{
}
#else
Point::Point(const Point & p)
{
    x = p.getX() ;
    y = p.getY() ;
    z = p.getZ() ;
    t = p.getT() ;
    id = p.getId() ;
}
#endif


    Point & operator*=( Point &p, const Matrix & m) {

        if(m.numCols() == 2)
        {
            Vector vec(2) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec = vec*m ;

            p.getX() = vec[0] ;
            p.getY() = vec[1] ;
        }
        else if(m.numCols() == 3)
        {
            Vector vec(3) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec[2] = p.getZ() ;
            vec = vec*m ;

            p.getX() = vec[0] ;
            p.getY() = vec[1] ;
            p.getZ() = vec[2] ;
        }
        else if(m.numCols() == 4)
        {
            Vector vec(4) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec[2] = p.getZ() ;
            vec[3] = p.getT() ;
            vec = vec*m ;

            p.getX() = vec[0] ;
            p.getY() = vec[1] ;
            p.getZ() = vec[2] ;
            p.getT() = vec[3] ;
        }
        
        return p ;
    }

   
   Point operator* (const Point & p, const Matrix & m) {
        if(m.numCols() == 2)
        {
            Vector vec(2) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec = vec*m ;
            return Point(vec[0], vec[1]) ;
        }
        if(m.numCols() == 3)
        {
            Vector vec(3) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec[2] = p.getZ() ;
            vec = vec*m ;
            return Point(vec[0], vec[1], vec[2]) ;
        }
        else if(m.numCols() == 4)
        {
            Vector vec(4) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec[2] = p.getZ() ;
            vec[3] = p.getT() ;
            vec = vec*m ;

            return Point(vec[0], vec[1], vec[2], vec[3]) ;
        }

        return Point() ;
    }
    
       Point operator* (const Matrix & m, const Point & p) {
        if(m.numCols() == 2)
        {
            Vector vec(2) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec = m*vec ;
            return Point(vec[0], vec[1]) ;
        }
        if(m.numCols() == 3)
        {
            Vector vec(3) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec[2] = p.getZ() ;
            vec = m*vec ;
            return Point(vec[0], vec[1], vec[2]) ;
        }
        else if(m.numCols() == 4)
        {
            Vector vec(4) ;
            vec[0] = p.getX() ;
            vec[1] = p.getY() ;
            vec[2] = p.getZ() ;
            vec[3] = p.getT() ;
            vec = m*vec ;

            return Point(vec[0], vec[1], vec[2], vec[3]) ;
        }

        return Point() ;
    }



Point & operator*=( Point & p, const double d) {
    p.getX() *= d ;
    p.getY() *= d ;
    p.getZ() *= d ;
    p.getT() *= d ;
    return p ;
}

Point & operator/=( Point & p, const double d) {
    double inv = 1./d ;
    p.getX() *= inv ;
    p.getY() *= inv ;
    p.getZ() *= inv ;
    p.getT() *= inv ;
    
    return p ;
}


Point& Point::operator = (const Point & p)
{
    getX() = p.getX() ;
    getY() = p.getY() ;
    getZ() = p.getZ() ;
    getT() = p.getT() ;
    getId() = p.getId() ;
    return *this ;
}


double Point::angle() const
{
    return atan2(y,x) ;
}

Point::Point(double x_, double y_, double z_, double t_): id(-1)
{
#ifdef HAVE_SSE3
    vecxy = _mm_setr_pd(x_, y_) ;
    veczt = _mm_setr_pd(z_, t_) ;
#else
    x= x_ ;
    y = y_ ;
    z = z_ ;
    t = t_ ;
#endif
}


void Point::print() const
{
    std::cout << " ( id = " << id << std::flush ;
    std::cout << " ; "<< x << std::flush ;
    std::cout << "; " << y << std::flush ;
    std::cout << "; " << z << std::flush ;
    std::cout << "; " << /*std::setprecision(16)<<*/ t << ") " << std::endl;
}

double Point::norm() const
{
// #ifdef HAVE_SSE4
// 	vecdouble r0 ;
// 	r0.vec = _mm_dp_pd(vecxy, vecxy, 61) ;
// 	r0.vec += _mm_dp_pd(veczt, veczt, 62) ;
// 	return sqrt(r0.val[0]+ r0.val[1]);
// #elif defined HAVE_SSE3
// 	vecdouble rzt ;
// 	rzt.vec = _mm_add_pd(_mm_mul_pd(veczt, veczt), _mm_mul_pd(vecxy, vecxy)) ;
// 	return sqrt(rzt.val[0]+ rzt.val[1]);
// #else
    return sqrt(x*x+y*y+z*z/*+t*t*/) ;
// #endif
}

double Point::sqNorm() const
{
// #ifdef HAVE_SSE4
// 	vecdouble r0 ;
// 	r0.vec = _mm_dp_pd(vecxy, vecxy, 61) ;
// 	r0.vec += _mm_dp_pd(veczt, veczt, 62) ;
// 	return r0.val[0]+ r0.val[1];
// #elif HAVE_SSE3
// 	vecdouble rzt ;
// 	rzt.vec = _mm_add_pd(_mm_mul_pd(veczt, veczt), _mm_mul_pd(vecxy, vecxy)) ;
// 	return rzt.val[0]+ rzt.val[1];
// #else
    return x*x+y*y+z*z/*+t*t*/ ;
// #endif
}

void Point::setX(double v)
{
    x = v ;
}

void Point::setY(double v)
{
    y = v ;
}

void Point::setZ(double v)
{
    z = v ;
}

void Point::setT(double v)
{
    t = v ;
}

void Point::set(const Point & p)
{
#ifdef HAVE_SSE3
    vecxy = p.vecxy ;
    veczt = p.veczt ;
#else
    x = p.getX() ;
    y = p.getY() ;
    z = p.getZ() ;
    t = p.getT() ;
#endif
}


void Point::set(double v, double vv, double vvv, double vvvv)
{
    x = v ;
    y = vv ;
    z = vvv ;
    t = vvvv ;
}

bool operator==(const Point & p_, const Point &p)
{


    if(std::abs(p.getX()-p_.getX()) > 2.*POINT_TOLERANCE)
        return false ;
    if(std::abs(p.getY()-p_.getY()) > 2.*POINT_TOLERANCE)
        return false ;
    if(std::abs(p.getZ()-p_.getZ()) > 2.*POINT_TOLERANCE)
        return false ;
    if(std::abs(p.getT()-p_.getT()) > 2.*POINT_TOLERANCE)
        return false ;

    return true ;

//     return dist(p, *this) < POINT_TOLERANCE ;
// 	Point mid = (p+ *this)*.5 ;
// 	double d =  std::max(p.norm(), norm()) ;
// 	if(d < POINT_TOLERANCE)
// 		return true ;
// 	return dist((p-mid)/d, (*this-mid)/d) <= POINT_TOLERANCE ;

// 	double delta = POINT_TOLERANCE ;
// 	Point a(p) ; a.getX() += delta ; a.getY() += delta ; a.getZ() += delta ;
// 	Point b(p) ; b.getX() += delta ; b.getY() += delta; b.getZ() -= delta ;
// 	Point c(p) ; c.getX() += delta ; c.getY() -= delta; c.getZ() += delta ;
// 	Point d(p) ; d.getX() += delta ; d.getY() -= delta; d.getZ() -= delta ;
// 	Point e(p) ; e.getX() -= delta ; e.getY() += delta; e.getZ() += delta ;
// 	Point f(p) ; f.getX() -= delta ; f.getY() += delta; f.getZ() -= delta ;
// 	Point g(p) ; g.getX() -= delta ; g.getY() -= delta; g.getZ() += delta ;
// 	Point h(p) ; h.getX() -= delta ; h.getY() -= delta; h.getZ() -= delta ;
//
// 	return squareDist( &p, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &a, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &b, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &c, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &d, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &e, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &f, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &g, this) < POINT_TOLERANCE*POINT_TOLERANCE
// 		|| squareDist( &h, this) < POINT_TOLERANCE*POINT_TOLERANCE ;
}

bool operator!=(const Point & p_ ,const Point & p)
{
    return !(p_ == p) ;
//     double pnorm = p.norm() ;
//     double tnorm = norm() ;
//     if(pnorm <= 4.*POINT_TOLERANCE && tnorm <= 4.*POINT_TOLERANCE)
//         return false ;
//     return dist(p, *this) >= 4.*POINT_TOLERANCE ;
}

Point operator-(const Point & p_, const Point &p)
{
    Point ret(p_) ;
#ifdef HAVE_SSE3
    ret.vecxy = _mm_sub_pd(ret.vecxy, p.vecxy) ;
    ret.veczt = _mm_sub_pd(ret.veczt, p.veczt) ;
    ret.setId( std::max(p_.getId(), p.getId())) ;
#else
    ret.getX() -= p.getX() ;
    ret.getY() -= p.getY() ;
    ret.getZ() -= p.getZ() ;
    ret.getT() -= p.getT() ;
#endif
    return ret ;
}

Point operator-(const Point & p_, const Vector &p)
{
    Point ret(p_) ;
    ret.getX() -= p[0] ;
    ret.getY() -= p[1] ;
    if(p.size() > 2)
        ret.getZ() -= p[2] ;
    if(p.size() > 3)
        ret.getT() -= p[3] ;
    ret.setId( p_.getId()) ;
    return ret ;
}
Point operator-(const Point & p_)
{
    Point ret(p_) ;
    ret.getX() = -p_.getX() ;
    ret.getY() = -p_.getY() ;
    
    ret.getZ() = -p_.getZ() ;
    ret.getT() = -p_.getT() ;
    ret.setId( p_.getId()) ;
    return ret ;
}

Point operator+(const Point & p_, const Point &p)
{
    Point ret(p_) ;
#ifdef HAVE_SSE3
    ret.vecxy = _mm_add_pd(ret.vecxy, p.vecxy) ;
    ret.veczt = _mm_add_pd(ret.veczt, p.veczt) ;
#else
    ret.getX() += p.getX() ;
    ret.getY() += p.getY() ;
    ret.getZ() += p.getZ() ;
    ret.getT() += p.getT() ;
#endif
    ret.getId() = std::max(p.getId(), p.getId()) ;
    return ret ;
}

Point & operator+=(Point &  p_, const Point &p)
{
#ifdef HAVE_SSE3
    p_.vecxy = _mm_add_pd(p_.vecxy, p.vecxy) ;
    p_.veczt = _mm_add_pd(p_.veczt, p.veczt) ;
#else
    p_.getX() += p.getX() ;
    p_.getY() += p.getY() ;
    p_.getZ() += p.getZ() ;
    p_.getT() += p.getT() ;
#endif
    return p_ ;
}

Point & operator-=(Point & p_,  const Point &p)
{
#ifdef HAVE_SSE3
    p_.vecxy = _mm_sub_pd(p_.vecxy, p.vecxy) ;
    p_.veczt = _mm_sub_pd(p_.veczt, p.veczt) ;
#else
    p_.getX() -= p.getX() ;
    p_.getY() -= p.getY() ;
    p_.getZ() -= p.getZ() ;
    p_.getT() -= p.getT() ;
#endif
    
    return p_ ;
}

Point operator+(const Point & p_, const Vector &p) 
{
    Point ret(p_) ;
    ret.getX() += p[0] ;
    ret.getY() += p[1] ;
    if(p.size() > 2)
        ret.getZ() += p[2] ;
    if(p.size() > 3)
        ret.getT() += p[3] ;
    ret.setId( p_.getId()) ;
    return ret ;
}

Point operator/(const Point & p_, const double p) 
{
    double inv = 1./p ;
    Point ret(p_.x*inv, p_.y*inv,p_.z*inv, p_.t*inv) ;
    ret.setId( p_.id) ;
    return ret ;
}

Geometry::~Geometry()
{

    for(size_t i = 0 ; i < inPoints.size() ; i++)
    {
        delete inPoints[i] ;
    }
}

bool operator <(const Point & p_, const Point &p)
{
// 	if(p == *this)
// 		return false ;
//
    if(p_.x < p.getX())
        return true ;
    else if(p_.x > p.getX())
        return false ;

    if(p_.y < p.getY())
        return true ;
    else if(p_.y > p.getY())
        return false ;

    if(p_.z < p.getZ())
        return true ;
    else if (p_.z > p.getZ())
        return false ;

    if(p_.t < p.getT())
        return true ;

    return false ;
}

bool operator >(const Point & p_, const Point &p)
{
    if(p == p_)
        return false ;

    double tol = POINT_TOLERANCE ;
    return (p_.y > p.getY() )
           || (( std::abs(p_.y - p.getY()) < tol)
               && (p_.x > p.getX()))
           || (( std::abs(p_.y - p.getY()) < tol)
               && ( std::abs(p_.x - p.getX()) < tol)
               && (p_.z> p.getZ()))
           ||(( std::abs(p_.y - p.getY()) < tol)
              && ( std::abs(p_.x - p.getX()) < tol)
              && ( std::abs(p_.z - p.getZ()) < tol)
              && (p_.t> p.getT()));
}

Point operator*(const Point & p_, const double p)
{
    Point ret(p_) ;
#ifdef HAVE_SSE3
    __m128d temp = _mm_load1_pd(&p) ;

    ret.vecxy = _mm_mul_pd(ret.vecxy, temp) ;
    ret.veczt = _mm_mul_pd(ret.veczt, temp) ;
#else
    ret.getX() *= p ;
    ret.getY() *= p ;
    ret.getZ() *= p ;
    ret.getT() *= p ;
#endif
    ret.setId( p_.id);
    return ret ;
}

Point operator*( const double p,const Point & p_)
{
    Point ret(p_) ;
#ifdef HAVE_SSE3
    __m128d temp = _mm_load1_pd(&p) ;

    ret.vecxy = _mm_mul_pd(ret.vecxy, temp) ;
    ret.veczt = _mm_mul_pd(ret.veczt, temp) ;
#else
    ret.getX() *= p ;
    ret.getY() *= p ;
    ret.getZ() *= p ;
    ret.getT() *= p ;
#endif
    ret.setId( p_.id);
    return ret ;
}

double operator*(const Point & p_, const Point &p)
{
// #ifdef HAVE_SSE4
//     vecdouble r ;
//     r.vec = _mm_dp_pd(p.vecxy, p_.vecxy, 61) ;
//     r.vec += _mm_dp_pd(p.veczt, p_.veczt, 62) ;
//     return r.val[0] + r.val[1];
// #elif defined HAVE_SSE3
//     vecdouble r ;
//     r.vec = _mm_add_pd(_mm_mul_pd(p.vecxy, vecxy), _mm_mul_pd(p.veczt, veczt)) ;
//     return r.val[0] + r.val[1];
// #endif
    return p.getX()*p_.x+p.getY()*p_.y+p.getZ()*p_.z/*+p.getT()*p_.t*/ ;

}

double operator*(const Point & p_, const Vector &p)
{
    double ret = p_.x*p[0] + p_.y*p[1] ;
    if(p.size() > 2)
        ret+=p_.z*p[2] ;
//     if(p.size() > 3)
//         ret+=p_.t*p[3] ;
    return ret ;
}


// Point Point::operator^(const Point &p) const
// {
// 	Point ret ;
//
// 	ret.getX() = y*p.getZ() - z*p.getY() ; //fma(y,p.getZ(),  -z*p.getY()) ;
// 	ret.getY() = z*p.getX() - x*p.getZ() ;//fma(z,p.getX() , -x*p.getZ()) ;
// 	ret.getZ() = x*p.getY() - y*p.getX() ; //fma(x,p.getY() , -y*p.getX()) ;
//
// 	ret.getId() = std::max(id, p.getId()) ;
// 	return ret ;
// }

// Point Point::operator^(const Vector &p) const
// {
// 	Point ret ;
// 	ret.getX() = y*p[2] - z*p[1] ;
// 	ret.getY() = z*p[0] - x*p[2] ;
// 	ret.getZ() = x*p[1] - y*p[0] ;
// 	ret.setId( id);
// 	return ret ;
// }

PtP operator^(const Point & p_, const Point &p)
{
    return PtP(p_,p) ;
}

PtV operator^(const Point & p_, const Vector &p)
{
    return PtV(p_,p) ;
}


PointSet::PointSet() : boundingPoints(0)
{
    this->chullEndPos = 0;
}

std::vector<Point> Geometry::getBoundingBox() const
{
    return std::vector<Point>() ;
}

Matrix rotationMatrix(double theta, size_t i)
{
    Matrix ret(4,4) ;
    ret[i][i] = 1 ;
    ret[3][3] = 1 ;

    int j = 0 ;
    int k = 1 ;
    if(i == 0)
    {
        j = 1 ;
        k = 2 ;
    }
    if(i == 1)
    {
        j = 2 ;
        k = 0 ;
    }

    ret[j][j] = cos(theta) ;
    ret[k][k] = cos(theta) ;
    ret[j][k] = -sin(theta) ;
    ret[k][j] = sin(theta) ;

    return ret ;
}

void transform(Geometry * g, GeometricTransformationType transformation, const Point& p)
{
    switch(transformation)
    {
    case SCALE:
        if( p.getX() < POINT_TOLERANCE || p.getY() < POINT_TOLERANCE || ( g->spaceDimensions() == SPACE_THREE_DIMENSIONAL && p.getZ() < POINT_TOLERANCE) )
        {
//            std::cout << "try to scale geometry with factor = 0 ... do nothing instead" << std::endl ;
            return ;
        }
        if(g->getGeometryType() == CIRCLE)
        {
            dynamic_cast<Circle *>(g)->setRadius( g->getRadius()*p.getX() );
            return ;
        }
        if(g->getGeometryType() == SPHERE)
        {
            dynamic_cast<Sphere *>(g)->setRadius( g->getRadius()*p.getX() );
            return ;
        }

        for(size_t i = 0 ; i < g->getInPoints().size() ; i++)
        {
            g->getInPoint(i).getX() = (g->getCenter().getX() + g->getInPoint(i).getX() - g->getCenter().getX()) * p.getX() ;
            g->getInPoint(i).getY() = (g->getCenter().getY() + g->getInPoint(i).getY() - g->getCenter().getY()) * p.getY() ;
            g->getInPoint(i).getZ() = (g->getCenter().getZ() + g->getInPoint(i).getZ() - g->getCenter().getZ()) * p.getZ() ;
        }
        for(size_t i = 0 ; i < g->getBoundingPoints().size() ; i++)
        {
            g->getBoundingPoint(i).getX() = g->getCenter().getX() + (g->getBoundingPoint(i).getX() - g->getCenter().getX()) * p.getX() ;
            g->getBoundingPoint(i).getY() = g->getCenter().getY() + (g->getBoundingPoint(i).getY() - g->getCenter().getY()) * p.getY() ;
            g->getBoundingPoint(i).getZ() = g->getCenter().getZ() + (g->getBoundingPoint(i).getZ() - g->getCenter().getZ()) * p.getZ() ;
        }

        if(g->getGeometryType() == ELLIPSE)
        {
            Point A = dynamic_cast<Ellipse *>(g)->getMajorAxis() ;
            Point B = dynamic_cast<Ellipse *>(g)->getMinorAxis() ;
            A.getX() *= p.getX() ;
            A.getY() *= p.getY() ;
            B.getX() *= p.getX() ;
            B.getY() *= p.getY() ;
            dynamic_cast<Ellipse *>(g)->setMajorAxis(A) ;
            dynamic_cast<Ellipse *>(g)->setMinorAxis(B) ;
        }

        if(g->getGeometryType() == POLYGON)
        {
	    std::valarray<Point> original = dynamic_cast<Polygon *>(g)->getOriginalPoints() ;
	    for(size_t i = 0 ; i < original.size() ; i++)
	    {
		original[i].setX( g->getCenter().getX() + (original[i].getX()-g->getCenter().getX()) * p.getX() ) ;
		original[i].setY( g->getCenter().getY() + (original[i].getY()-g->getCenter().getY()) * p.getY() ) ;
		original[i].setZ( g->getCenter().getZ() + (original[i].getZ()-g->getCenter().getZ()) * p.getZ() ) ;
	    }
	    dynamic_cast<Polygon *>(g)->setOriginalPoints(original, true) ;
        }

        break ;
    case ROTATE:
    {

        if(g->getGeometryType() == CIRCLE)
        {
            return ;
        }
        if(g->getGeometryType() == SPHERE)
        {
            return ;
        }

        Matrix rotation = rotationMatrix( p.getX(), 0) ;
        rotation *= rotationMatrix( p.getY(), 1 ) ;
        rotation *= rotationMatrix( p.getZ(), 2 ) ;

        Point c = g->getCenter() ;
        c *= -1 ;
        transform(g, TRANSLATE, c);

        for(size_t i = 0 ; i < g->getInPoints().size() ; i++)
            g->getInPoint(i) *= rotation ;
        for(size_t i = 0 ; i < g->getBoundingPoints().size() ; i++)
            g->getBoundingPoint(i) *= rotation ;

        if(g->getGeometryType() == POLYGON)
        {
            std::valarray<Point> pts = dynamic_cast<Polygon *>(g)->getOriginalPoints() ;
            for(size_t i = 0 ; i < pts.size() ; i++)
               pts[i] *= rotation ;
            dynamic_cast<Polygon *>(g)->setOriginalPoints(pts, true) ;
        }

        c *= -1 ;
        transform(g, TRANSLATE, c);


        if(g->getGeometryType() == ELLIPSE)
        {
            Point A = dynamic_cast<Ellipse *>(g)->getMajorAxis() ;
            Point B = dynamic_cast<Ellipse *>(g)->getMinorAxis() ;
            A *= rotation ;
            B *= rotation;
            dynamic_cast<Ellipse *>(g)->setMajorAxis(A) ;
            dynamic_cast<Ellipse *>(g)->setMinorAxis(B) ;
        }
        break ;
    }
    case TRANSLATE:
        g->setCenter( g->getCenter() + p ) ;
        break ;
    }
}

const std::valarray<Point *> & ConvexGeometry::getBoundingPoints() const
{
    return boundingPoints ;
}

std::valarray<Point *> & ConvexGeometry::getBoundingPoints()
{
    return boundingPoints ;
}

const Point & ConvexGeometry::getBoundingPoint(size_t i) const
{
    return *boundingPoints[i] ;
}

Point & ConvexGeometry::getBoundingPoint(size_t i)
{
    return *boundingPoints[i] ;
}

const Point & ConvexGeometry::getPoint(size_t i) const
{
    if (i < inPoints.size())
    {
        return *inPoints[i] ;
    }

    return *boundingPoints[i-inPoints.size()] ;
}

Point & ConvexGeometry::getPoint(size_t i)
{
    if (i < inPoints.size())
    {
        return *inPoints[i] ;
    }

    return *boundingPoints[i-inPoints.size()] ;
}

std::valarray<Point *> & NonConvexGeometry::getBoundingPoints()
{
    return boundingPoints ;
}

const std::valarray<Point *> & NonConvexGeometry::getBoundingPoints() const
{
    return boundingPoints ;
}

const Point & NonConvexGeometry::getBoundingPoint(size_t i) const
{
    return *boundingPoints[i] ;
}

Point & NonConvexGeometry::getBoundingPoint(size_t i)
{
    return *boundingPoints[i] ;
}

const Point & Geometry::getInPoint(size_t i) const
{
    return *inPoints[i] ;
}


void Geometry::setCenter(const Point & newCenter)
{
    Point delta = newCenter-getCenter();
    getCenter() = newCenter ;

    for(size_t  i = 0 ; i < getInPoints().size() ; i++)
        getInPoint(i) += delta ;
    for(size_t  i = 0 ; i < getBoundingPoints().size() ; i++)
        getBoundingPoint(i) += delta ;
}

Point & Geometry::getInPoint(size_t i)
{
    return *inPoints[i] ;
}

void NonConvexGeometry::setBoundingPoint(size_t i, Point * p)
{
    boundingPoints[i] = p ;
}

void NonConvexGeometry::setBoundingPoints(const PointArray & nb)
{
    boundingPoints.resize(nb.size()) ;
    boundingPoints = nb ;
}

void ConvexGeometry::setBoundingPoint(size_t i, Point * p)
{
    boundingPoints[i] = p ;
}

void ConvexGeometry::setBoundingPoints(const PointArray & nb)
{
    boundingPoints.resize(nb.size()) ;
    boundingPoints = nb ;
}

Plane::Plane(const Point & origin, const Point & vector) : p(origin), v(vector)
{
    v /= v.norm() ;
}

Plane::Plane(const Point & a, const Point & b, const Point & c) : p(a), v((b-a)^(c-a))
{
    v /= v.norm() ;
}

bool Plane::intersects(const Line &l) const
{
    double d = p*v ;
    double numerator = d - l.origin()*v ;
    double denominator = l.vector()*v ;

    // line is in the plane
    if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
        return true ;

    // the intersection exists and is unique
    if(std::abs(denominator) > POINT_TOLERANCE)
        return true ;

    return false ;
}

bool Plane::intersects(const Segment &l) const
{
    double d = p*v ;
    double numerator = d - l.second()*v ;
    double denominator = l.vector()*v ;

    // segment is in the plane
    if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
        return true ;

    // the intersection exists and is unique
    if(std::abs(denominator) > POINT_TOLERANCE)
    {
        double t = numerator/denominator ;

        return t > 0 && t < 1 ;
    }

    return false ;
}

bool Plane::intersects(const Geometry *g) const
{
    switch (g->getGeometryType())
    {
    case SPHERE:
    {
        double n = v.norm() ;
        double d = v*p ;
        d /= n ;

        return d >= g->getRadius() ;

    }
    case HEXAHEDRON:
    {
        std::vector<Point> bbox = g->getBoundingBox() ;
        double maxx =  bbox[0].getX() ;
        double minx =  bbox[0].getX() ;
        double maxy =  bbox[0].getY() ;
        double miny =  bbox[0].getY() ;
        double maxz =  bbox[0].getZ() ;
        double minz =  bbox[0].getZ() ;

        for(size_t i = 1 ; i < bbox.size() ; i++)
        {
            if(bbox[i].getX() > maxx)
                maxx = bbox[i].getX() ;
            if(bbox[i].getX() < minx)
                minx = bbox[i].getX() ;
            if(bbox[i].getY() > maxy)
                maxy = bbox[i].getY() ;
            if(bbox[i].getY() < miny)
                miny = bbox[i].getY() ;
            if(bbox[i].getZ() > maxz)
                maxz = bbox[i].getZ() ;
            if(bbox[i].getZ() < minz)
                minz = bbox[i].getZ() ;
        }

        Point corner1 (minx, miny, minz) ;
        Point corner2 (minx, miny, maxz) ;
        Point corner3 (minx, maxy, minz) ;
        Point corner4 (minx, maxy, maxz) ;
        Point corner5 (maxx, miny, minz) ;
        Point corner6 (maxx, miny, maxz) ;
        Point corner7 (maxx, maxy, minz) ;
        Point corner8 (maxx, maxy, maxz) ;

        Segment side1(corner1, corner2) ;
        if(intersects(side1))
            return true ;
        Segment side2(corner1, corner3) ;
        if(intersects(side2))
            return true ;
        Segment side3(corner1, corner5) ;
        if(intersects(side3))
            return true ;
        Segment side4(corner2, corner4) ;
        if(intersects(side4))
            return true ;
        Segment side5(corner2, corner6) ;
        if(intersects(side5))
            return true ;
        Segment side6(corner3, corner4) ;
        if(intersects(side6))
            return true ;
        Segment side7(corner3, corner7) ;
        if(intersects(side7))
            return true ;
        Segment side8(corner4, corner8) ;
        if(intersects(side8))
            return true ;
        Segment side9(corner5, corner6) ;
        if(intersects(side9))
            return true ;
        Segment side10(corner5, corner7) ;
        if(intersects(side10))
            return true ;
        Segment side11(corner6, corner8) ;
        if(intersects(side11))
            return true ;
        Segment side12(corner7, corner8) ;
        if(intersects(side12))
            return true ;

        return false ;
    }
    case TETRAHEDRON:
    {
        Point corner1 (g->getBoundingPoint(0)) ;
        Point corner2 (g->getBoundingPoint(1)) ;
        Point corner3 (g->getBoundingPoint(2)) ;
        Point corner4 (g->getBoundingPoint(3)) ;

        Segment side1(corner1, corner2) ;
        if(intersects(side1))
            return true ;
        Segment side2(corner1, corner3) ;
        if(intersects(side2))
            return true ;
        Segment side3(corner2, corner4) ;
        if(intersects(side3))
            return true ;
        Segment side4(corner3, corner4) ;
        if(intersects(side4))
            return true ;

        return false ;
    }
    default:
    {
        std::cout << "intersection type not implemented" << std::endl ;
        return false ;
    }
    }
}

bool Plane::intersects(const Plane & plane) const
{
    return v/v.norm() != plane.vector()/plane.vector().norm() ;
}

bool Plane::on(const Point &point) const
{
    double d = v*p ;
    double dtest = point*v ;
    return std::abs(d-dtest) < POINT_TOLERANCE ;
}

std::vector<Point> Plane::intersection(const Geometry * g) const
{
    switch (g->getGeometryType())
    {
    case SPHERE:
    {
        Point vec = projection(g->getCenter())-g->getCenter() ;
        double l = dist(g->getCenter(),projection(g->getCenter())) ;
        Point centerOfIntersection = g->getCenter() + vec/vec.norm()*l ;
        double radiusOfIntersection = sqrt(g->getRadius()*g->getRadius() - l*l) ;
        OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

        size_t num_points = std::max((int)round(23.*sqrt(g->getBoundingPoints().size())*radiusOfIntersection/g->getRadius()), 8 ) ;
        C.sampleSurface(num_points) ;

        return C.getSamplingBoundingPoints(num_points) ;

    }
    case HEXAHEDRON:
    {
        std::vector<Point> ret ;
        std::vector<Point> bbox = g->getBoundingBox() ;
        double maxx =  bbox[0].getX() ;
        double minx =  bbox[0].getX() ;
        double maxy =  bbox[0].getY() ;
        double miny =  bbox[0].getY() ;
        double maxz =  bbox[0].getZ() ;
        double minz =  bbox[0].getZ() ;

        for(size_t i = 1 ; i < bbox.size() ; i++)
        {
            if(bbox[i].getX() > maxx)
                maxx = bbox[i].getX() ;
            if(bbox[i].getX() < minx)
                minx = bbox[i].getX() ;
            if(bbox[i].getY() > maxy)
                maxy = bbox[i].getY() ;
            if(bbox[i].getY() < miny)
                miny = bbox[i].getY() ;
            if(bbox[i].getZ() > maxz)
                maxz = bbox[i].getZ() ;
            if(bbox[i].getZ() < minz)
                minz = bbox[i].getZ() ;
        }

        Point corner1 (minx, miny, minz) ;
        Point corner2 (minx, miny, maxz) ;
        Point corner3 (minx, maxy, minz) ;
        Point corner4 (minx, maxy, maxz) ;
        Point corner5 (maxx, miny, minz) ;
        Point corner6 (maxx, miny, maxz) ;
        Point corner7 (maxx, maxy, minz) ;
        Point corner8 (maxx, maxy, maxz) ;

        Segment side1(corner1, corner2) ;
        if(intersects(side1))
            ret.push_back(intersection(side1)) ;
        Segment side2(corner1, corner3) ;
        if(intersects(side2))
            ret.push_back(intersection(side2)) ;
        Segment side3(corner1, corner5) ;
        if(intersects(side3))
            ret.push_back(intersection(side3)) ;
        Segment side4(corner2, corner4) ;
        if(intersects(side4))
            ret.push_back(intersection(side4)) ;
        Segment side5(corner2, corner6) ;
        if(intersects(side5))
            ret.push_back(intersection(side5)) ;
        Segment side6(corner3, corner4) ;
        if(intersects(side6))
            ret.push_back(intersection(side6)) ;
        Segment side7(corner3, corner7) ;
        if(intersects(side7))
            ret.push_back(intersection(side7)) ;
        Segment side8(corner4, corner8) ;
        if(intersects(side8))
            ret.push_back(intersection(side8)) ;
        Segment side9(corner5, corner6) ;
        if(intersects(side9))
            ret.push_back(intersection(side9)) ;
        Segment side10(corner5, corner7) ;
        if(intersects(side10))
            ret.push_back(intersection(side10)) ;
        Segment side11(corner6, corner8) ;
        if(intersects(side11))
            ret.push_back(intersection(side11)) ;
        Segment side12(corner7, corner8) ;
        if(intersects(side12))
            ret.push_back(intersection(side12)) ;

        return ret ;
    }
    case TETRAHEDRON:
    {
        std::vector<Point> ret ;
        Point corner1 (g->getBoundingPoint(0)) ;
        Point corner2 (g->getBoundingPoint(1)) ;
        Point corner3 (g->getBoundingPoint(2)) ;
        Point corner4 (g->getBoundingPoint(3)) ;

        Segment side1(corner1, corner2) ;
        if(intersects(side1))
            ret.push_back(intersection(side1)) ;
        Segment side2(corner1, corner3) ;
        if(intersects(side2))
            ret.push_back(intersection(side2)) ;
        Segment side3(corner2, corner4) ;
        if(intersects(side3))
            ret.push_back(intersection(side3)) ;
        Segment side4(corner3, corner4) ;
        if(intersects(side4))
            ret.push_back(intersection(side4)) ;

        return ret ;
    }
    default:
    {
        std::cout << "intersection type not implemented" << std::endl ;
        return std::vector<Point>() ;
    }
    }
}

Point Plane::intersection(const Line &l) const
{
    double lambda = 0  ;

    if(std::abs(v*(l.origin() + l.vector() - p )) < POINT_TOLERANCE)
    {
        return l.origin() ;
    }
    else if(std::abs(v*l.vector()) >= POINT_TOLERANCE)
    {
        lambda = -(v*(l.origin()-p))/(v*l.vector()) ;
        return l.origin() + l.vector()*lambda ;
    }
    return Point() ;

// 	if(std::abs(l.origin().getX()*l.vector().getX()) > POINT_TOLERANCE)
// 	{
// 		lambda = ( v.getX() * p.getX() - v.getX()*l.origin().getX() ) / (l.origin().getX()*l.vector().getX()) ;
// 		return l.origin() + l.vector()*lambda ;
// 	}
// 	else if (std::abs(l.origin().getY()*l.vector().getY()) > POINT_TOLERANCE)
// 	{
// 		lambda = ( v.getY() * p.getY() - v.getY()*l.origin().getY() ) / (l.origin().getY()*l.vector().getY()) ;
// 		return l.origin() + l.vector()*lambda ;
// 	}
// 	else
// 		return Point() ;
//
// 	double d = p*v ;
// 	double numerator = d - l.origin()*v ;
// 	double denominator = l.vector()*v ;
//
// 	// segment is in the plane
// 	if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
// 		return p ;
//
// 	//there is no intersection, but we need to return something
// 	return p ;
//
// 	// the intersection exists and is unique
// 	double t = numerator/denominator ;
// 	return l.origin()+l.vector()*t ;

}

Point Plane::intersection(const Segment &l) const
{
    double d = p*v ;
    double numerator = d - l.second()*v ;
    double denominator = l.vector()*v ;

    // segment is in the plane
    if(std::abs(numerator) < POINT_TOLERANCE && std::abs(denominator) < POINT_TOLERANCE)
        return l.midPoint() ;

    // the intersection exists and is unique
    if(std::abs(denominator) > POINT_TOLERANCE)
    {
        double t = numerator/denominator ;

        return l.second()+l.vector()*t ;
    }

    return p ;
}

Line Plane::intersection(const Plane & plane) const
{
    Point vec = v^plane.vector() ;
    double d1 = p*v ;
    double d2 = plane.vector()*plane.origin() ;
    if(!((v.getX() != 0  && plane.vector().getY() != 0) || (v.getY() != 0  && plane.vector().getX() != 0)))
    {
        Matrix m(2, 2) ;
        m[0][0] = v.getX() ;
        m[0][1] = v.getY() ;
        m[1][0] = plane.vector().getX() ;
        m[1][1] = plane.vector().getY() ;

        invert2x2Matrix(m) ;

        Vector d (2) ;
        d[0] = d1 ;
        d[1] = d2 ;

        Vector xy = m*d ;

        return Line(Point(xy[0], xy[1], 0), vec) ;
    }
    else if(!((v.getX() != 0  && plane.vector().getZ() != 0) || (v.getZ() != 0  && plane.vector().getX() != 0)))
    {
        Matrix m(2, 2) ;
        m[0][0] = v.getX() ;
        m[0][1] = v.getZ() ;
        m[1][0] = plane.vector().getX() ;
        m[1][1] = plane.vector().getZ() ;

        invert2x2Matrix(m) ;

        Vector d (2) ;
        d[0] = d1 ;
        d[1] = d2 ;

        Vector xz = m*d ;

        return Line(Point(xz[0], 0, xz[1]), vec) ;
    }
    else if(!((v.getY() != 0  && plane.vector().getZ() != 0) || (v.getZ() != 0  && plane.vector().getY() != 0)))
    {
        Matrix m(2, 2) ;
        m[0][0] = v.getY() ;
        m[0][1] = v.getZ() ;
        m[1][0] = plane.vector().getX() ;
        m[1][1] = plane.vector().getZ() ;

        invert2x2Matrix(m) ;

        Vector d (2) ;
        d[0] = d1 ;
        d[1] = d2 ;

        Vector yz = m*d ;

        return Line(Point(0, yz[0], yz[1]), vec) ;
    }

    return Line(Point(0, 0, 0), vec) ;
}

const Point & Plane::vector() const
{
    return v ;
}

const Point & Plane::origin() const
{
    return p ;
}

Point Plane::projection(const Point &toProject ) const
{
    return intersection(Line(toProject, v)) ;
}

double Geometry::overlapFraction(const Geometry * target) const
{
    Point c = getCenter() ;
    if(dynamic_cast<const Triangle *>(this))
        c = dynamic_cast<const Triangle *>(this)->getCircumCenter() ;
    else if(dynamic_cast<const Tetrahedron *>(this))
        c = dynamic_cast<const Tetrahedron *>(this)->getCircumCenter() ;

    Point tc = target->getCenter() ;
    if(dynamic_cast<const Triangle *>(target))
        tc = dynamic_cast<const Triangle *>(target)->getCircumCenter() ;
    else if(dynamic_cast<const Tetrahedron *>(target))
        tc = dynamic_cast<const Tetrahedron *>(target)->getCircumCenter() ;


    if(dist(c, tc)  > getRadius() + target->getRadius())
        return 0 ;

    if(dist(c, tc) + getRadius() <= target->getRadius())
        return 1 ;

    if(spaceDimensions() == SPACE_TWO_DIMENSIONAL)
    {
        double selfin = 0 ;
        double targin = 0 ;
        for(double i = c.getX()-getRadius() ; i < c.getX()+getRadius() ; i += getRadius()/8.)
        {
            for(double j = c.getY()-getRadius() ; j < c.getY()+getRadius() ; j += getRadius()/8.)
            {
                Point p(i,j) ;
                if(in(p))
                {
                    selfin++ ;
                    if(target->in(p))
                        targin++ ;
                }
            }
        }

        return targin/selfin ;
    }
    else
    {
        double selfin = 0 ;
        double targin = 0 ;
        for(double i = c.getX()-getRadius() ; i < c.getX()+getRadius() ; i += getRadius()/8.)
        {
            for(double j = c.getY()-getRadius() ; j < c.getY()+getRadius() ; j += getRadius()/8.)
            {
                for(double k = c.getZ()-getRadius() ; k < c.getZ()+getRadius() ; k += getRadius()/8.)
                {
                    Point p(i,j,k) ;
                    if(in(p))
                    {
                        selfin++ ;
                        if(target->in(p))
                            targin++ ;
                    }
                }
            }
        }

        return targin/selfin ;
    }
}

const std::valarray<Point* > & Geometry::getInPoints() const
{
    return inPoints ;
}

std::valarray<Point* > & Geometry::getInPoints()
{
    return inPoints ;
}

bool Geometry::intersects(const Geometry *g) const
{
    switch(getGeometryType())
    {
    case TRIANGLE :
    {
        int numpoints = getBoundingPoints().size() ;
        Segment s0(getBoundingPoint(0), getBoundingPoint(numpoints/3)) ;
        Segment s1(getBoundingPoint(numpoints/3), getBoundingPoint(2*numpoints/3)) ;
        Segment s2(getBoundingPoint(0),getBoundingPoint(2*numpoints/3)) ;

        if(g->getGeometryType() == TRIANGLE)
        {
            int gnumpoints = g->getBoundingPoints().size() ;
            Segment t0(g->getBoundingPoint(0), g->getBoundingPoint(gnumpoints/3)) ;
            Segment t1(g->getBoundingPoint(gnumpoints/3), g->getBoundingPoint(2*gnumpoints/3)) ;
            Segment t2(g->getBoundingPoint(0), g->getBoundingPoint(2*gnumpoints/3)) ;

            return s0.intersects(t0) || s0.intersects(t1) || s0.intersects(t2)
                   || s1.intersects(t0) || s1.intersects(t1) || s1.intersects(t2)
                   || s2.intersects(t0) || s2.intersects(t1) || s2.intersects(t2) ;

        }

        return s0.intersects(g) || s1.intersects(g)  || s2.intersects(g) ;
    }
    case RECTANGLE:
    {
        std::vector<Point> box = getBoundingBox() ;
        Segment s0(box[0], box[1]) ;
        Segment s1(box[1], box[2]) ;
        Segment s2(box[2], box[3]) ;
        Segment s3(box[3], box[0]) ;

        return s0.intersects(g) || s1.intersects(g)  || s2.intersects(g)  || s3.intersects(g);
    }
    case POLYGON:
    {
	std::valarray<Point> original = dynamic_cast<const Polygon *>(this)->getOriginalPoints() ;
	int inext = 0 ;
	for(size_t i = 0 ; i < original.size() ; i++)
	{
		inext = (i+1)%original.size() ;
		Segment s(original[i], original[inext]) ;
		if(s.intersects(g))
			return true ;
	}
	return false ;
    }
    case PARALLELOGRAMME:
    {
        size_t n = getBoundingPoints().size() ;

        std::vector<Point> box ;
        box.push_back(getBoundingPoint(0)) ;
        box.push_back(getBoundingPoint(n*1/4)) ;
        box.push_back(getBoundingPoint(n*2/4)) ;
        box.push_back(getBoundingPoint(n*3/4)) ;
        Segment s0(box[0], box[1]) ;
        Segment s1(box[1], box[2]) ;
        Segment s2(box[2], box[3]) ;
        Segment s3(box[3], box[0]) ;

        return s0.intersects(g) || s1.intersects(g)  || s2.intersects(g)  || s3.intersects(g);
    }
    case SEGMENTED_LINE:
    {
        std::vector<Segment> segs ;
        for(size_t i = 0 ; i < getBoundingPoints().size()-1 ; i++)
        {
            segs.push_back(Segment(getBoundingPoint(i), getBoundingPoint(i+1))) ;
        }
        bool intersects = false ;
        for(size_t i = 0 ; i < segs.size() ; i++)
        {
            intersects = intersects || segs[i].intersects(g) ;
        }
        return intersects ;
    }
    case CIRCLE:
    {
        if(g->getGeometryType() == CIRCLE)
        {
            double birad = getRadius()+g->getRadius() ;
            if(getCenter().getX() + birad < g->getCenter().getX())
                return false ;
            if(getCenter().getX() - birad > g->getCenter().getX())
                return false ;
            if(getCenter().getY() + birad < g->getCenter().getY())
                return false ;
            if(getCenter().getY() - birad > g->getCenter().getY())
                return false ;

            return squareDist2D( getCenter(), g->getCenter()) < birad*birad ;
        }

        std::vector<Segment> segs ;

        if(g->getGeometryType() == ELLIPSE)
	{
	    Point p(getCenter().getX(), getCenter().getY());
            Ellipse falseCircle( p , Point(this->getRadius(),0.), Point(0.,this->getRadius() ) ) ;
            return g->intersects(&falseCircle) ;
	}

        if(g->getGeometryType() == TRIANGLE)
        {
            segs.push_back(Segment(g->getBoundingPoint(0),
                                   g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3),
                                   g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(0),
                                   g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
        }
        if(g->getGeometryType() == RECTANGLE)
        {
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(3*g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(3*g->getBoundingPoints().size()/4))) ;
        }
        if(g->getGeometryType() == SEGMENTED_LINE)
        {
            for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
            {
                segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
            }
        }
        if(g->getGeometryType() == POLYGON)
        {
            std::valarray<Point> original = dynamic_cast<const Polygon *>(g)->getOriginalPoints() ;
            for(size_t i = 0 ; i < original.size()-1 ; i++)
            {
                segs.push_back(Segment(original[i], original[i+1])) ;
            }
            segs.push_back(Segment(original[original.size()-1], original[0])) ;
        }

        bool intersects = false ;
        for(size_t i = 0 ; i < segs.size() ; i++)
        {
            intersects = intersects || segs[i].intersects(this) ;
        }

       return intersects ;
    }
    case TIME_DEPENDENT_CIRCLE:
    {
        if(g->timePlanes() < 2)
            return false ;
        if(g->getGeometryType() == TRIANGLE)
        {

            bool in = false ;
            bool out = false ;

            for(size_t i = 0 ; i < g->getBoundingPoints().size() ; i++)
            {
                if(this->in(g->getBoundingPoint(i)))
                    in = true ;
                else
                    out = false ;

                if(in && out)
                    return true ;
            }


            size_t pointsPerPlane = g->getBoundingPoints().size() / g->timePlanes();
            size_t pointsPerEdge = pointsPerPlane/3 ;

            for(size_t i = 0 ; i < g->timePlanes() ; i++)
            {

                Segment s0( g->getBoundingPoint(i*pointsPerPlane), g->getBoundingPoint(i*pointsPerPlane+ pointsPerEdge)) ;
                Segment s1( g->getBoundingPoint(i*pointsPerPlane+pointsPerEdge), g->getBoundingPoint(i*pointsPerPlane+ pointsPerEdge*2)) ;
                Segment s2( g->getBoundingPoint(i*pointsPerPlane+pointsPerEdge*2), g->getBoundingPoint(i*pointsPerPlane)) ;
                if(s0.intersects(this) || s1.intersects(this) || s2.intersects(this))
                {
                    return true ;
                }
            }
        }

        return false ;
    }
    case ELLIPSE:
    {

        if(g->getRadius() < this->getRadius() /*&& g->getGeometryType() != POLYGON*/)
        {
            return g->intersects(this) ;
        }

        if(g->getGeometryType() == CIRCLE)
        {
	    Point p(g->getCenter().getX(), g->getCenter().getY()) ;
            Ellipse falseCircle(p, Point(g->getRadius(), 0.), Point(0.,g->getRadius())) ;
            return falseCircle.intersects(this) ;
        }

        if(g->getGeometryType() == ELLIPSE)
        {
            if(dist(g->getCenter(),getCenter()) > (getRadius()+g->getRadius())+POINT_TOLERANCE)
	    {
                return false ;
	    }

            if(dist(g->getCenter(), getCenter()) < dynamic_cast<const Ellipse *>(this)->getMinorRadius() || dist(g->getCenter(), getCenter()) < dynamic_cast<const Ellipse *>(g)->getMinorRadius())
                return true ;

            Point gcenter(g->getCenter().getX(), g->getCenter().getY()) ;
            Point thiscenter(getCenter().getX(), getCenter().getY()) ;
            Point ga(dynamic_cast<const Ellipse *>(g)->getMajorAxis().getX(), dynamic_cast<const Ellipse *>(g)->getMajorAxis().getY()) ;
            Point gb(dynamic_cast<const Ellipse *>(g)->getMinorAxis().getX(), dynamic_cast<const Ellipse *>(g)->getMinorAxis().getY()) ;
            Point thisa(dynamic_cast<const Ellipse *>(this)->getMajorAxis().getX(), dynamic_cast<const Ellipse *>(this)->getMajorAxis().getY()) ;
            Point thisb(dynamic_cast<const Ellipse *>(this)->getMinorAxis().getX(), dynamic_cast<const Ellipse *>(this)->getMinorAxis().getY()) ;

            Ellipse gcopy(gcenter, ga, gb) ;
            gcopy.sampleBoundingSurface(64) ;

            Ellipse thiscopy(thiscenter, thisa, thisb) ;
            thiscopy.sampleBoundingSurface(64) ;

		if(thiscopy.in(gcenter))
		{
			return true ;
		}

            for(size_t i = 0 ; i < gcopy.getBoundingPoints().size() ; i++)
            {
                if(thiscopy.in(gcopy.getBoundingPoint(i)))
		{
			return true ;
		}
            }


            for(size_t i = 0 ; i < thiscopy.getBoundingPoints().size() ; i++)
            {
                if(gcopy.in(thiscopy.getBoundingPoint(i)))
		{
			return true ;
		}
            }


            return false ;

        }

        std::vector<Segment> segs ;
        bool isInSegments = false ;

        if(g->getGeometryType() == TRIANGLE)
        {
            segs.push_back(Segment(g->getBoundingPoint(0),
                                   g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3),
                                   g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(0),
                                   g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
            isInSegments = true ;
        }
        if(g->getGeometryType() == RECTANGLE)
        {
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(3*g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(3*g->getBoundingPoints().size()/4))) ;
            isInSegments = true ;
        }
        if(g->getGeometryType() == SEGMENTED_LINE)
        {
            for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
            {
                segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
            }
            isInSegments = true ;
        }
        if(g->getGeometryType() == POLYGON)
        {
            std::valarray<Point> vertex = dynamic_cast<const Polygon*>(g)->getOriginalPoints() ;
            for(size_t i = 0 ; i < vertex.size() ; i++)
            {
                segs.push_back(Segment(vertex[i], vertex[(i+1)%vertex.size()])) ;
            }
            isInSegments = true ;
        }

        if(isInSegments)
        {
            for(size_t i = 0 ; i < segs.size() ; i++)
            {
                if(segs[i].intersects(this))
                    return true ;
            }
            return false ;
        }

        if(g->getGeometryType() == ELLIPSE && (g->in(getCenter()) || in(g->getCenter())))
            return true ;

        std::vector<Point> box ;
        for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
            box.push_back(getBoundingPoint(i)) ;
        if(box.size() == 0)
            box = this->getBoundingBox() ;
        box.push_back(box[0]) ;
        bool intersects = false ;
        for(size_t i = 0 ; i < box.size()-1 ; i++) {
            Segment seg(box[i],box[i+1]) ;
            intersects = intersects || seg.intersects(g) ;
        }

        return intersects ;

    }
    case SPHERE:
    {
        if(g->getGeometryType() == SPHERE)
        {
            if(in(g->getCenter()) || g->in(getCenter()))
                return true ;

            double birad = getRadius()+g->getRadius() ;
// 				if(((getCenter().getX() + birad < g->getCenter().getX()) ||
// 					return false ;
// 				if(getCenter().getX() - birad > g->getCenter().getX()))
// 					return false ;
// 				if((getCenter().getY() + birad < g->getCenter().getY()) ||
// 					return false ;
// 				if(getCenter().getY() - birad > g->getCenter().getY()))
// 					return false ;
// 				if((getCenter().getZ() + birad < g->getCenter().getZ()) ||
// 					return false ;
// 				if(getCenter().getZ() - birad > g->getCenter().getZ())))
// 					return false ;

            return dist(getCenter(), g->getCenter()) < birad ;

        }
        if(g->getGeometryType() == HEXAHEDRON)
        {
            if((getCenter().getX() < g->getCenter().getX()-static_cast<const Hexahedron *>(g)->getXSize()*.5 +getRadius()&&
               getCenter().getX() > g->getCenter().getX()-static_cast<const Hexahedron *>(g)->getXSize()*.5 -getRadius())||
               (getCenter().getX() < g->getCenter().getX()+static_cast<const Hexahedron *>(g)->getXSize()*.5 +getRadius()&&
               getCenter().getX() > g->getCenter().getX()+static_cast<const Hexahedron *>(g)->getXSize()*.5 -getRadius())||
               (getCenter().getY() < g->getCenter().getY()-static_cast<const Hexahedron *>(g)->getYSize()*.5 +getRadius()&&
               getCenter().getY() > g->getCenter().getY()-static_cast<const Hexahedron *>(g)->getYSize()*.5 -getRadius())||
               (getCenter().getY() < g->getCenter().getY()+static_cast<const Hexahedron *>(g)->getYSize()*.5 +getRadius()&&
               getCenter().getY() > g->getCenter().getY()+static_cast<const Hexahedron *>(g)->getYSize()*.5 -getRadius())||
               (getCenter().getZ() < g->getCenter().getZ()-static_cast<const Hexahedron *>(g)->getZSize()*.5 +getRadius()&&
               getCenter().getZ() > g->getCenter().getZ()-static_cast<const Hexahedron *>(g)->getZSize()*.5 -getRadius())|| 
               (getCenter().getZ() < g->getCenter().getZ()+static_cast<const Hexahedron *>(g)->getZSize()*.5 +getRadius()&&
               getCenter().getZ() > g->getCenter().getZ()+static_cast<const Hexahedron *>(g)->getZSize()*.5 -getRadius())
            )
                return true ;
            
            return false ;
            
            Point proj0(g->getCenter()) ;
            project(&proj0) ;
            Point proj1(getCenter()) ;
            g->project(&proj1) ;
            return g->in(proj0) || in(proj1) ;
        }
        if(g->getGeometryType() == TETRAHEDRON)
        {
            if(squareDist3D( static_cast<const Tetrahedron*>(g)->getCircumCenter(), getCenter()) > (g->getRadius()+getRadius())*(g->getRadius()+getRadius()))
                return false ;

            int incount = 0 ;
            int outcount = 0 ;
            for(size_t i = 0 ; i < g->getBoundingPoints().size() ; i++)
            {
                if(in(g->getBoundingPoint(i)))
                    incount++ ;
                else
                    outcount++ ;
            }

            if(incount && outcount)
                return true ;

            std::multimap<double, const Point *> pts ;
            for(size_t i = 0 ; i < g->getBoundingPoints().size() ; i++)
            {
                pts.insert(std::make_pair(std::abs(g->getRadius()*g->getRadius()-squareDist3D(static_cast<const Tetrahedron*>(g)->getCircumCenter(), g->getBoundingPoint(i))), & g->getBoundingPoint(i))) ;
            }
            auto p = pts.begin() ;
            const Point * a = p->second ;
            ++p ;
            const Point * b = p->second ;
            ++p ;
            const Point * c = p->second ;
            ++p ;
            const Point * d = p->second ;

            Segment s(g->getCenter(), getCenter()) ;
            TriPoint t0(a, b, c) ;
            TriPoint t1(a, b, d) ;
            TriPoint t2(a, c, d) ;
            TriPoint t3(b, c, d) ;
            if(s.intersects(t0))
                return in(s.intersection(t0)[0]) ;
            if(s.intersects(t1))
                return in(s.intersection(t1)[0]) ;
            if(s.intersects(t2))
                return in(s.intersection(t2)[0]) ;
            if(s.intersects(t3))
                return in(s.intersection(t3)[0]) ;

            return false ;

        }
        return false ;

    }
    case HEXAHEDRON:
    {
        if(g->getGeometryType() == SPHERE)
        {
            return g->intersects(this) ;
        }
        if(g->getGeometryType() == HEXAHEDRON)
        {
            size_t count = 0 ;
            size_t count_alt = 0 ;

            std::vector<Point> myBB = getBoundingBox() ;
            std::vector<Point> gBB = g->getBoundingBox() ;

            for(size_t i = 0 ; i < 8 ; i++)
            {
                if(in(gBB[i]))
                    count++ ;
                if(g->in(myBB[i]))
                    count_alt++ ;
            }

            return count != 0 && count !=8 && count_alt != 0 && count_alt !=8;

        }
        if(g->getGeometryType() == TETRAHEDRON)
            return g->intersects(this) ;
        return false ;
    }
    case TETRAHEDRON:
    {
        if(g->getGeometryType() == SPHERE)
        {
		std::cout << "really???" << std::endl ;
            return g->intersects(this) ;
        }

        if(g->getGeometryType() == HEXAHEDRON)
        {
            if(g->in(getCenter()) || in(g->getCenter()))
                return true ;

            Segment s(getCenter(), g->getCenter()) ;

            std::vector<Point> inter = s.intersection(this) ;
            for(size_t i = 0 ; i < inter.size() ; i++)
                if(g->in(inter[i]))
                    return true ;
            inter = s.intersection(g) ;
            for(size_t i = 0 ; i < inter.size() ; i++)
                if(in(inter[i]))
                    return true ;
        }
    }
    default:
    {
        return false ;
    }
    }
}

std::vector<Point> Geometry::intersection(const Geometry * g) const
{
    std::vector<Point> ret ;
    switch(getGeometryType())
    {
    case POLYGON:
    {
        std::valarray<Point> original = dynamic_cast<const Polygon *>(this)->getOriginalPoints() ;
        int inext = 0 ;
        for(size_t i = 0 ; i < original.size() ; i++)
        {
            inext = (i+1)%original.size() ;
            Segment s(original[i], original[inext]) ;
            if(s.intersects(g))
            {
                std::vector<Point> inter = s.intersection(g) ;
                for(size_t j = 0 ; j < inter.size() ; j++)
                {
                    bool alone = true ;
                    for(size_t k = 0 ; alone && k < ret.size() ; k++)
                    {
                        if(dist(ret[k],inter[j]) < POINT_TOLERANCE)
                            alone = false ;
                    }
                    if(alone)
                        ret.push_back(inter[j]) ;
                }
            }
        }
        return ret ;
    }
    case TRIANGLE :
    {
        Segment s0(getBoundingPoint(0),
                   getBoundingPoint(getBoundingPoints().size()/3)) ;
        Segment s1(getBoundingPoint(getBoundingPoints().size()/3),
                   getBoundingPoint(2*getBoundingPoints().size()/3)) ;
        Segment s2(getBoundingPoint(0),
                   getBoundingPoint(2*getBoundingPoints().size()/3)) ;
        std::vector<Point> intersection = s0.intersection(g) ;
        ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
        intersection = s1.intersection(g) ;
        ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
        intersection = s2.intersection(g) ;
        ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
        return ret ;
    }
    case RECTANGLE:
    {
        std::vector<Point> box = getBoundingBox() ;
        Segment s0(box[0], box[1]) ;
        Segment s1(box[1], box[2]) ;
        Segment s2(box[2], box[3]) ;
        Segment s3(box[3], box[0]) ;
        if(true)//g->getGeometryType() == RECTANGLE)
        {
            std::vector<Point> intersection = s0.intersection(g) ;
            std::vector<Point> it = s1.intersection(g) ;
            intersection.insert(intersection.end(), it.begin(), it.end()) ;
            it = s2.intersection(g) ;
            intersection.insert(intersection.end(), it.begin(), it.end()) ;
            it = s3.intersection(g) ;
            intersection.insert(intersection.end(), it.begin(), it.end()) ;

            bool haveDuplicates = true ;
            while(haveDuplicates)
            {
                haveDuplicates = false ;
                for(size_t i  = 0 ; i < intersection.size() ; i++)
                {
                    for(size_t j  = i+1 ; j < intersection.size() ; j++)
                    {
                        if(squareDist3D(intersection[i], intersection[j])< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                        {
                            haveDuplicates = true ;
                            intersection.erase(intersection.begin()+j) ;
                            break ;
                        }
                    }
			if(intersection.size() > 1)
				ret.push_back(intersection[intersection.size()-1]) ;

                    if(haveDuplicates)
                        break ;
                }
            }
            return intersection ;
        }

        if(g->getGeometryType() == POLYGON)
        {
            std::vector<Point> intersection = s0.intersection(g) ;
            if(s0.intersects(g) && s1.intersects(g) && g->in(box[1]))
                intersection.push_back(box[1]) ;
            std::vector<Point> it = s1.intersection(g) ;
            intersection.insert(intersection.end(), it.begin(), it.end()) ;
            if(s1.intersects(g) && s2.intersects(g) && g->in(box[2]))
                intersection.push_back(box[2]) ;
            it = s2.intersection(g) ;
            intersection.insert(intersection.end(), it.begin(), it.end()) ;
            if(s2.intersects(g) && s3.intersects(g) && g->in(box[3]))
                intersection.push_back(box[3]) ;
            if(s3.intersects(g) && s0.intersects(g) && g->in(box[0]))
                intersection.push_back(box[0]) ;
            it = s3.intersection(g) ;
            intersection.insert(intersection.end(), it.begin(), it.end()) ;

            if(intersection.size() > 0)
            {
                for(size_t i = 0 ; i < box.size() ; i++)
                {
                   if( g->in(box[i]) )
                       intersection.push_back(box[i]) ;
                }
            }

            bool haveDuplicates = true ;
            while(haveDuplicates)
            {
                haveDuplicates = false ;
                for(size_t i  = 0 ; i < intersection.size() ; i++)
                {
                    for(size_t j  = i+1 ; j < intersection.size() ; j++)
                    {
                        if(squareDist3D(intersection[i], intersection[j])< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                        {
                            haveDuplicates = true ;
                            intersection.erase(intersection.begin()+j) ;
                            break ;
                        }
                    }
			if(intersection.size() > 1)
				ret.push_back(intersection[intersection.size()-1]) ;

                    if(haveDuplicates)
                        break ;
                }
            }
            return intersection ;
        }

	std::vector<Point> intersection = s0.intersection(g) ;
	if(intersection.size() == 1)
	{
	    if(g->in(s0.second()))
	        intersection.push_back(s0.second()) ;
	    else if(g->in(s0.first()))
	        intersection.push_back(s0.first()) ;
	}
	double perimetre = s0.norm()+s1.norm()+s2.norm()+s3.norm() ;
	for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
	{
		ret.push_back(intersection[i]) ;
		double l = dist(intersection[i], intersection[i+1]) ;
		double jmax = round(1.5*getBoundingPoints().size()*(l/perimetre)) ;
		for(double j = 1 ; j < jmax+1 ; j++ )
		{
			double frac = j/(jmax+1) ;
			ret.push_back(intersection[i]*(1.-frac)+intersection[i+1]*frac) ;
		}
	}
	if(intersection.size() > 1)
		ret.push_back(intersection[intersection.size()-1]) ;


        intersection = s1.intersection(g) ;
        if(intersection.size() == 1)
        {
            if(g->in(s1.second()))
                intersection.push_back(s1.second()) ;
            else if(g->in(s1.first()))
                intersection.push_back(s1.first()) ;
        }
        for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
        {
            ret.push_back(intersection[i]) ;
            double l = dist(intersection[i], intersection[i+1]) ;
            double jmax = round(1.5*getBoundingPoints().size()*(l/perimetre)) ;
            for(double j = 1 ; j < jmax+1 ; j++ )
            {
                double frac = j/(jmax+1) ;
                ret.push_back(intersection[i]*(1.-frac)+intersection[i+1]*frac) ;
            }

        }
	if(intersection.size() > 1)
		ret.push_back(intersection[intersection.size()-1]) ;

        intersection = s2.intersection(g) ;
        if(intersection.size() == 1)
        {
            if(g->in(s2.second()))
                intersection.push_back(s2.second()) ;
            else if(g->in(s2.first()))
                intersection.push_back(s2.first()) ;
        }

        for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
        {
            ret.push_back(intersection[i]) ;
            double l = dist(intersection[i], intersection[i+1]) ;
            double jmax = round(1.5*getBoundingPoints().size()*(l/perimetre)) ;
            for(double j = 1 ; j < jmax+1 ; j++ )
            {
                double frac = j/(jmax+1) ;
                ret.push_back(intersection[i]*(1.-frac)+intersection[i+1]*frac) ;
            }
        }
	if(intersection.size() > 1)
		ret.push_back(intersection[intersection.size()-1]) ;

        intersection = s3.intersection(g) ;
        if(intersection.size() == 1)
        {
            if(g->in(s3.second()))
                intersection.push_back(s3.second()) ;
            else if(g->in(s3.first()))
                intersection.push_back(s3.first()) ;
        }
        for(int i = 0 ; i < (int)intersection.size()-1 ; i++)
        {
            ret.push_back(intersection[i]) ;
            double l = dist(intersection[i], intersection[i+1]) ;
            double jmax = round(1.5*getBoundingPoints().size()*(l/perimetre)) ;
            for(double j = 1 ; j < jmax+1 ; j++ )
            {
                double frac = j/(jmax+1) ;
                ret.push_back(intersection[i]*(1.-frac)+intersection[i+1]*frac) ;
            }
        }
	if(intersection.size() > 1)
		ret.push_back(intersection[intersection.size()-1]) ;


        bool haveDuplicates = true ;
        while(haveDuplicates)
        {
            haveDuplicates = false ;
            for(size_t i  = 0 ; i < ret.size() ; i++)
            {
                for(size_t j  = i+1 ; j < ret.size() ; j++)
                {
                    if(squareDist3D(ret[i], ret[j])< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                    {
                        haveDuplicates = true ;
                        ret.erase(ret.begin()+j) ;
                        break ;
                    }
                }

                if(haveDuplicates)
                    break ;
            }
        }
        return ret ;
    }
    case SEGMENTED_LINE:
    {
        std::vector<Segment> segs ;
        for(size_t i = 0 ; i < getBoundingPoints().size()-1 ; i++)
        {
            segs.push_back(Segment(getBoundingPoint(i), getBoundingPoint(i+1))) ;
        }
        for(size_t i = 0 ; i < segs.size() ; i++)
        {
            std::vector<Point> intersection = segs[i].intersection(g) ;
            ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
        }
        return ret ;
    }
    case CIRCLE:
    {
        if(g->getGeometryType() == CIRCLE)
        {
            double r = getRadius() ;
            double R = g->getRadius() ;


            //first quadratic equation
            double a = 1. ;
            double b = 1. - 2.*R*R ;
            double c = R*R*R*R - r*r ;

            double delta = b*b - 4.*a*c ;

            if(delta < 0)
                return ret ;

            double y_squared_0 = (-b + sqrt(delta))/(2.*a) ;
            if(delta == 0)
            {
                if(y_squared_0 < 0)
                    return ret ;



                if((r*r-y_squared_0) >= 0)
                {
                    double y = sqrt(y_squared_0) ;
                    double x = sqrt(r*r-y_squared_0) ;
                    Point A(x,y) ;
                    if(std::abs(squareDist(A, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(A, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(A) ;

                    Point B(-x,y) ;
                    if(std::abs(squareDist(B, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(B, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(B) ;

                    Point C(-x,-y) ;
                    if(std::abs(squareDist(C, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(C, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(C) ;

                    Point D(x,-y) ;
                    if(std::abs(squareDist(D, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(D, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(D) ;


                }

                return ret ;
            }

            double y_squared_1 = (-b - sqrt(delta))/(2.*a) ;

            if(y_squared_0 >= 0)
            {
                if((r*r-y_squared_0) >= 0)
                {
                    double x = sqrt(r*r-y_squared_0) ;
                    double y = sqrt(y_squared_0) ;

                    Point A(x,y) ;
                    if(std::abs(squareDist(A, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(A, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(A) ;

                    Point B(-x,y) ;
                    if(std::abs(squareDist(B, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(B, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(B) ;

                    Point C(-x,-y) ;
                    if(std::abs(squareDist(C, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(C, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(C) ;

                    Point D(x,-y) ;
                    if(std::abs(squareDist(D, getCenter()) -r*r) < POINT_TOLERANCE*POINT_TOLERANCE &&
                            std::abs(squareDist(D, g->getCenter()) -R*R) < POINT_TOLERANCE*POINT_TOLERANCE
                      )
                        ret.push_back(D) ;
                }
            }
            if(y_squared_1 >= 0)
            {
                if((r*r-y_squared_1) >= 0)
                {

                    double x = sqrt(r*r-y_squared_1) ;
                    double y = sqrt(y_squared_1) ;

                    Point A(x,y) ;
                    if(std::abs(squareDist(A, getCenter()) -r*r) < POINT_TOLERANCE &&
                            std::abs(squareDist(A, g->getCenter()) -R*R) < POINT_TOLERANCE
                      )
                        ret.push_back(A) ;

                    Point B(-x,y) ;
                    if(std::abs(squareDist(B, getCenter()) -r*r) < POINT_TOLERANCE &&
                            std::abs(squareDist(B, g->getCenter()) -R*R) < POINT_TOLERANCE
                      )
                        ret.push_back(B) ;

                    Point C(-x,-y) ;
                    if(std::abs(squareDist(C, getCenter()) -r*r) < POINT_TOLERANCE &&
                            std::abs(squareDist(C, g->getCenter()) -R*R) < POINT_TOLERANCE
                      )
                        ret.push_back(C) ;

                    Point D(x,-y) ;
                    if(std::abs(squareDist(D, getCenter()) -r*r) < POINT_TOLERANCE &&
                            std::abs(squareDist(D, g->getCenter()) -R*R) < POINT_TOLERANCE
                      )
                        ret.push_back(D) ;
                }
            }

            return ret ;


        }

        std::vector<Segment> segs ;

        if(g->getGeometryType() == TRIANGLE)
        {
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3), g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
        }
        if(g->getGeometryType() == RECTANGLE || g->getGeometryType() == POLYGON)
        {
            std::vector<Point> inter = g->intersection(this) ;
//            if(inter.size() < 2 || getBoundingPoints().size() < 2)
                return inter ;
/*            for(size_t i = 0 ; i < inter.size()-1 ; i++)
            {
                double num = getBoundingPoints().size() ;
                Point start = inter[i]-getCenter() ;
                Point end = inter[i+1]-getCenter() ;
                double a = start.angle() ;
                double b = end.angle() ;
                if(b > a)
                    b -= M_PI*2. ;
                double alpha = a-b ;
                num *= std::abs(alpha)/(2.*M_PI) ;
                std::vector<Point> it = dynamic_cast<const Circle *>(this)->getSamplingBoundingPointsOnArc( round(num), inter[i], inter[i+1] ) ;
                if(it.size() < 2)
                    return inter ;
                if( g->in(it[it.size()/2]) )
                {
                    ret.insert(ret.end(), it.begin(), it.end()) ;
                    if(i == inter.size()-2)
                       ret.push_back(inter[i+1]) ;
                }
                else
                {
                    it.clear() ;
                    it = dynamic_cast<const Circle *>(this)->getSamplingBoundingPointsOnArc( getBoundingPoints().size()-round(num), inter[i], inter[i+1], true ) ;
                    if( g->in(it[it.size()/2]) )
                    {
			std::cout << "REVERSE" << std::endl ;
                        ret.insert(ret.end(), it.begin(), it.end()) ;
                        if(i == inter.size()-2)
                           ret.push_back(inter[i+1]) ;
                     }
                     else
                     {
			std::cout << "WRONG REVERSE" << std::endl ;
                        it.clear() ;
                        it = dynamic_cast<const Circle *>(this)->getSamplingBoundingPointsOnArc( getBoundingPoints().size()-round(num), inter[i+1], inter[i], true ) ;
                        if( !g->in(it[it.size()/2]) )
				std::cout << "WRONG REVERSE (LINE 2)" << std::endl ;
                        ret.insert(ret.end(), it.begin(), it.end()) ;
                        if(i == inter.size()-2)
                           ret.push_back(inter[i]) ;
                     }
                }
            }
            return ret ;*/
        }
        if(g->getGeometryType() == SEGMENTED_LINE)
        {
            for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
            {
                segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
            }
        }

        for(size_t i = 0 ; i < segs.size() ; i++)
        {
            std::vector<Point>intersection = segs[i].intersection(this) ;
            ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
        }
        return ret ;
    }
    case TIME_DEPENDENT_CIRCLE:
    {
        if(g->getGeometryType() == TRIANGLE && g->timePlanes() > 1)
        {
            size_t pointsPerPlane = g->getBoundingPoints().size() / g->timePlanes() ;
            size_t pointsPerEdge = pointsPerPlane / 3 ;
            std::vector<Point> inter ;
            for(size_t i = 0 ; i < g->timePlanes() ; i++)
            {
                Segment s0(g->getBoundingPoint(pointsPerPlane*i), g->getBoundingPoint(pointsPerPlane*i+pointsPerEdge)) ;
                inter.clear() ;
                inter = s0.intersection(this) ;
                ret.insert(ret.end(), inter.begin(), inter.end() ) ;
                Segment s1(g->getBoundingPoint(pointsPerPlane*i+pointsPerEdge), g->getBoundingPoint(pointsPerPlane*i+pointsPerEdge*2)) ;
                inter.clear() ;
                inter = s1.intersection(this) ;
                ret.insert(ret.end(), inter.begin(), inter.end() ) ;
                Segment s2(g->getBoundingPoint(pointsPerPlane*i+pointsPerEdge*2), g->getBoundingPoint(pointsPerPlane*i)) ;
                inter.clear() ;
                inter = s2.intersection(this) ;
                ret.insert(ret.end(), inter.begin(), inter.end() ) ;
            }

        }
        return ret ;
    }
    case ELLIPSE:
    {

        std::vector<Segment> segs ;
        bool isInSegments = false ;

        if(g->getGeometryType() == TRIANGLE)
        {
            segs.push_back(Segment(g->getBoundingPoint(0),
                                   g->getBoundingPoint(g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/3),
                                   g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
            segs.push_back(Segment(g->getBoundingPoint(0),
                                   g->getBoundingPoint(2*g->getBoundingPoints().size()/3))) ;
            isInSegments = true ;
        }
        if(g->getGeometryType() == RECTANGLE)
        {
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(3*g->getBoundingPoints().size()/4), g->getBoundingPoint(2*g->getBoundingPoints().size()/4))) ;
            segs.push_back(Segment(g->getBoundingPoint(0), g->getBoundingPoint(3*g->getBoundingPoints().size()/4))) ;
            isInSegments = true ;
        }
        if(g->getGeometryType() == SEGMENTED_LINE)
        {
            for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
            {
                segs.push_back(Segment(g->getBoundingPoint(i), g->getBoundingPoint(i+1))) ;
            }
            isInSegments = true ;
        }
        if(g->getGeometryType() == POLYGON)
        {
            std::valarray<Point> vertex = dynamic_cast<const Polygon *>(g)->getOriginalPoints() ;
            for(size_t i = 0 ; i < vertex.size() ; i++)
            {
                segs.push_back(Segment(vertex[i], vertex[(i+1)%vertex.size()])) ;
            }
            isInSegments = true ;
        }

        if(isInSegments)
        {
            for(size_t i = 0 ; i < segs.size() ; i++)
            {
                std::vector<Point> vtemp = segs[i].intersection(this) ;
                for(size_t j = 0 ; j < vtemp.size() ; j++)
                {
                    ret.push_back(vtemp[j]) ;
                }
            }
            return ret ;
        }

        for(size_t i = 0 ; i < getBoundingPoints().size() - 1 ; i++)
        {
            segs.push_back(Segment(getBoundingPoint(i),getBoundingPoint(i+1))) ;
        }
        segs.push_back(Segment(getBoundingPoint(getBoundingPoints().size()-1),getBoundingPoint(0))) ;

        for(size_t i = 0 ; i < segs.size() ; i++)
        {
            std::vector<Point>intersection = segs[i].intersection(g) ;
            for(size_t j = 0 ; j < intersection.size() ; j++)
            {
//					intersection[j].print() ;
                ret.push_back(intersection[j]) ;
            }
//				ret.insert(ret.end(), intersection.begin(), intersection.end()) ;
        }
        return ret ;

    }
    case SPHERE:
    {

        if(g->getGeometryType() == SPHERE)
        {
            double dc = dist(getCenter(), g->getCenter()) ;
            double r0 = getRadius() ;
            double r1 = g->getRadius() ;

            if (dc < POINT_TOLERANCE)
                return ret ;

            double cosAngleInSelf = ((r0*r0+dc*dc-r1*r1)/(2.*r0*dc)) ;

            double d = cosAngleInSelf*r0/*0.5*getRadius() +0.5*g->getRadius() - dc*/;

            Point v = g->getCenter()-getCenter() ;

            Point centerOfIntersection = getCenter() + v*(d/dc) ;
            double radiusOfIntersection = sqrt(r0*r0 - d*d) ;

            OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

            size_t num_points = (size_t)(2*sqrt(g->getBoundingPoints().size())) ;

            if(num_points < 3 )
                return ret ;

            C.sampleBoundingSurface(num_points) ;
            for(size_t i = 0 ;  i < num_points ; i++)
            {
                ret.push_back(C.getBoundingPoint(i)) ;
            }

            return ret ;

        }
        if(g->getGeometryType() == HEXAHEDRON)
        {
            std::vector<Point> bbox = g->getBoundingBox() ;
            double maxx =  bbox[0].getX() ;
            double minx =  bbox[0].getX() ;
            double maxy =  bbox[0].getY() ;
            double miny =  bbox[0].getY() ;
            double maxz =  bbox[0].getZ() ;
            double minz =  bbox[0].getZ() ;

            for(size_t i = 1 ; i < bbox.size() ; i++)
            {
                if(bbox[i].getX() > maxx)
                    maxx = bbox[i].getX() ;
                if(bbox[i].getX() < minx)
                    minx = bbox[i].getX() ;
                if(bbox[i].getY() > maxy)
                    maxy = bbox[i].getY() ;
                if(bbox[i].getY() < miny)
                    miny = bbox[i].getY() ;
                if(bbox[i].getZ() > maxz)
                    maxz = bbox[i].getZ() ;
                if(bbox[i].getZ() < minz)
                    minz = bbox[i].getZ() ;
            }



            if(std::abs(center.getX() - minx) < getRadius())
            {
                double d = std::abs(center.getX() - minx) ;
                Point v(-1,0,0) ;
                Point centerOfIntersection(minx, center.getY(), center.getZ()) ;
                double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
                OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

                size_t num_points = std::max(round(0.5*radiusOfIntersection/getBoundingPoints().size()/(getRadius()*getRadius())), 3.) ;

                C.sampleSurface(num_points) ;

                Circle planeCircle(radiusOfIntersection, Point( center.getY(), center.getZ())) ;
                planeCircle.sampleSurface(num_points) ;
                Rectangle planeRect(maxy-miny, maxz-minz, g->getCenter().getY(), g->getCenter().getZ()) ;
                planeRect.sampleSurface(g->getBoundingPoints().size()/2) ;
                std::vector<Point> planeIntersection = planeRect.intersection(&planeCircle) ;
                for(size_t i = 0 ;  i < planeIntersection.size() ; i++)
                {
                    Point candidate(minx, planeIntersection[i].getX(), planeIntersection[i].getY()) ;
                    if(g->in(candidate) && in(candidate))
                        ret.push_back(candidate) ;
                }


                for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
                {
                    if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
                        ret.push_back(C.getBoundingPoint(i)) ;

                }
                for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
                {
                    if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
                        ret.push_back(C.getInPoint(i)) ;
                }

            }
            if(std::abs(center.getX() - maxx) < getRadius())
            {
                double d = std::abs(center.getX() - maxx) ;
                Point v(1,0,0) ;
                Point centerOfIntersection(maxx, center.getY(), center.getZ()) ;
                double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
                OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

                size_t num_points = std::max(round(0.5*radiusOfIntersection/getBoundingPoints().size()/(getRadius()*getRadius())), 3.) ;
                C.sampleSurface(num_points) ;

                Circle planeCircle(radiusOfIntersection, Point( center.getY(), center.getZ())) ;
                planeCircle.sampleSurface(num_points) ;
                Rectangle planeRect(maxy-miny, maxz-minz, g->getCenter().getY(), g->getCenter().getZ()) ;
                planeRect.sampleSurface(g->getBoundingPoints().size()/2) ;
                std::vector<Point> planeIntersection = planeRect.intersection(&planeCircle) ;
                for(size_t i = 0 ;  i < planeIntersection.size() ; i++)
                {

                    Point candidate(maxx, planeIntersection[i].getX(), planeIntersection[i].getY()) ;
                    if(g->in(candidate) && in(candidate))
                        ret.push_back(candidate) ;
                }

                for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
                {
                    if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
                        ret.push_back(C.getBoundingPoint(i)) ;

                }
                for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
                {
                    if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
                        ret.push_back(C.getInPoint(i)) ;
                }
            }
            if(std::abs(center.getY() - miny) < getRadius())
            {
                double d = std::abs(center.getY() - miny) ;
                Point v(0,-1,0) ;
                Point centerOfIntersection(center.getX(), miny, center.getZ()) ;
                double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
                OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

                size_t num_points = std::max(round(0.5*radiusOfIntersection/getBoundingPoints().size()/(getRadius()*getRadius())), 3.) ;
                C.sampleSurface(num_points) ;

                Circle planeCircle(radiusOfIntersection, Point( center.getX(), center.getZ())) ;
                planeCircle.sampleSurface(num_points) ;
                Rectangle planeRect(maxx-minx, maxz-minz, g->getCenter().getX(), g->getCenter().getZ()) ;
                planeRect.sampleSurface(g->getBoundingPoints().size()/2) ;
                std::vector<Point> planeIntersection = planeRect.intersection(&planeCircle) ;
                for(size_t i = 0 ;  i < planeIntersection.size() ; i++)
                {
                    Point candidate(planeIntersection[i].getX(), miny, planeIntersection[i].getY()) ;
                    if(g->in(candidate) && in(candidate))
                        ret.push_back(candidate) ;
                }

                for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
                {
                    if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
                        ret.push_back(C.getBoundingPoint(i)) ;

                }
                for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
                {
                    if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
                        ret.push_back(C.getInPoint(i)) ;
                }
            }
            if(std::abs(center.getY() - maxy) < getRadius())
            {
                double d = std::abs(center.getY() - maxy) ;
                Point v(0,1,0) ;
                Point centerOfIntersection(center.getX(), maxy, center.getZ()) ;
                double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
                OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

                size_t num_points = std::max(round(0.5*radiusOfIntersection/getBoundingPoints().size()/(getRadius()*getRadius())), 3.) ;
                C.sampleSurface(num_points) ;

                Circle planeCircle(radiusOfIntersection, Point( center.getX(), center.getZ())) ;
                planeCircle.sampleSurface(num_points) ;
                Rectangle planeRect(maxx-minx, maxz-minz, g->getCenter().getX(), g->getCenter().getZ()) ;
                planeRect.sampleSurface(g->getBoundingPoints().size()/2) ;
                std::vector<Point> planeIntersection = planeRect.intersection(&planeCircle) ;
                for(size_t i = 0 ;  i < planeIntersection.size() ; i++)
                {
                    Point candidate(planeIntersection[i].getX(), maxy, planeIntersection[i].getY()) ;
                    if(g->in(candidate) && in(candidate))
                        ret.push_back(candidate) ;
                }

                for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
                {
                    if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
                        ret.push_back(C.getBoundingPoint(i)) ;

                }
                for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
                {
                    if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
                        ret.push_back(C.getInPoint(i)) ;
                }
            }
            if(std::abs(center.getZ() - minz) < getRadius())
            {
                double d = std::abs(center.getZ() - minz) ;
                Point v(0,0,-1) ;
                Point centerOfIntersection(center.getX(), center.getY(), minz) ;
                double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
                OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

                size_t num_points = std::max(round(0.5*radiusOfIntersection/getBoundingPoints().size()/(getRadius()*getRadius())), 3.) ;
                C.sampleSurface(num_points) ;

                Circle planeCircle(radiusOfIntersection, Point( center.getX(), center.getY())) ;
                planeCircle.sampleSurface(num_points) ;
                Rectangle planeRect(maxx-minx, maxy-miny, g->getCenter().getX(), g->getCenter().getY()) ;
                planeRect.sampleSurface(g->getBoundingPoints().size()/2) ;
                std::vector<Point> planeIntersection = planeRect.intersection(&planeCircle) ;
                for(size_t i = 0 ;  i < planeIntersection.size() ; i++)
                {
                    Point candidate(planeIntersection[i].getX(), planeIntersection[i].getY(), minz) ;
                    if(g->in(candidate) && in(candidate))
                        ret.push_back(candidate) ;
                }

                for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
                {
                    if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
                        ret.push_back(C.getBoundingPoint(i)) ;

                }
                for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
                {
                    if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
                        ret.push_back(C.getInPoint(i)) ;
                }
            }
            if(std::abs(center.getZ() - maxz) < getRadius())
            {
                double d = std::abs(center.getZ() - maxz) ;
                Point v(0,0,1) ;
                Point centerOfIntersection(center.getX(), center.getY(), maxz) ;
                double radiusOfIntersection = sqrt(getRadius()*getRadius() - d*d) ;
                OrientableCircle C(radiusOfIntersection, centerOfIntersection, v) ;

                size_t num_points = std::max(round(0.5*radiusOfIntersection/getBoundingPoints().size()/(getRadius()*getRadius())), 3.) ;
                C.sampleSurface(num_points) ;

                Circle planeCircle(radiusOfIntersection, Point( center.getX(), center.getY())) ;
                planeCircle.sampleSurface(num_points) ;
                Rectangle planeRect(maxx-minx, maxy-miny, g->getCenter().getX(), g->getCenter().getY()) ;
                planeRect.sampleSurface(g->getBoundingPoints().size()/2) ;
                std::vector<Point> planeIntersection = planeRect.intersection(&planeCircle) ;
                for(size_t i = 0 ;  i < planeIntersection.size() ; i++)
                {
                    Point candidate(planeIntersection[i].getX(), planeIntersection[i].getY(), maxz) ;
                    if(g->in(candidate) && in(candidate))
                        ret.push_back(candidate) ;
                }

                for(size_t i = 0 ;  i < C.getBoundingPoints().size() ; i++)
                {
                    if(g->in(C.getBoundingPoint(i)) && in(C.getBoundingPoint(i)))
                        ret.push_back(C.getBoundingPoint(i)) ;

                }
                for(size_t i = 0 ;  i < C.getInPoints().size() ; i++)
                {
                    if(g->in(C.getInPoint(i)) && in(C.getInPoint(i)))
                        ret.push_back(C.getInPoint(i)) ;
                }
            }

            bool haveDuplicates = true ;
            while(haveDuplicates)
            {
                haveDuplicates = false ;
                for(size_t i  = 0 ; i < ret.size() ; i++)
                {
                    for(size_t j  = i+1 ; j < ret.size() ; j++)
                    {
                        if(squareDist3D(ret[i], ret[j])< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                        {
                            haveDuplicates = true ;
                            ret.erase(ret.begin()+j) ;
                            break ;
                        }
                    }

                    if(haveDuplicates)
                        break ;
                }
            }

        }
        if(g->getGeometryType() == TETRAHEDRON)
        {
            int mult = 1 ;
            if(g->getBoundingPoints().size() == 10)
                mult = 2 ;
            Segment s0(g->getBoundingPoint(0), g->getBoundingPoint(1*mult)) ;
            Segment s1(g->getBoundingPoint(0), g->getBoundingPoint(2*mult)) ;
            Segment s2(g->getBoundingPoint(0), g->getBoundingPoint(3*mult)) ;
            Segment s3(g->getBoundingPoint(1*mult), g->getBoundingPoint(2*mult)) ;
            Segment s4(g->getBoundingPoint(1*mult), g->getBoundingPoint(3*mult)) ;
            Segment s5(g->getBoundingPoint(2*mult), g->getBoundingPoint(3*mult)) ;
            std::vector<Segment *> intersectingSegments ;
            if(s0.intersects(this))
            {
                std::vector<Point> inter = s0.intersection(this) ;
                if(inter.size() == 1)
                    ret.insert(ret.end(), inter.begin(), inter.end()) ;
            }
            if(s1.intersects(this))
            {
                std::vector<Point> inter = s1.intersection(this) ;
                if(inter.size() == 1)
                    ret.insert(ret.end(), inter.begin(), inter.end()) ;
            }
            if(s2.intersects(this))
            {
                std::vector<Point> inter = s2.intersection(this) ;
                if(inter.size() == 1)
                    ret.insert(ret.end(), inter.begin(), inter.end()) ;
            }
            if(s3.intersects(this))
            {
                std::vector<Point> inter = s3.intersection(this) ;
                if(inter.size() == 1)
                    ret.insert(ret.end(), inter.begin(), inter.end()) ;
            }
            if(s4.intersects(this))
            {
                std::vector<Point> inter = s4.intersection(this) ;
                if(inter.size() == 1)
                    ret.insert(ret.end(), inter.begin(), inter.end()) ;
            }
            if(s5.intersects(this))
            {
                std::vector<Point> inter = s5.intersection(this) ;
                if(inter.size() == 1)
                    ret.insert(ret.end(), inter.begin(), inter.end()) ;
            }

            bool haveDuplicates = true ;
            while(haveDuplicates)
            {
                haveDuplicates = false ;
                for(size_t i  = 0 ; i < ret.size() ; i++)
                {
                    for(size_t j  = i+1 ; j < ret.size() ; j++)
                    {
                        if(squareDist3D(ret[i], ret[j])< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                        {
                            haveDuplicates = true ;
                            ret.erase(ret.begin()+j) ;
                            break ;
                        }
                    }

                    if(haveDuplicates)
                        break ;
                }
            }

// 				TriPoint t0(&g->getBoundingPoint(0), &g->getBoundingPoint(1*mult), &g->getBoundingPoint(2*mult)) ;
// 				TriPoint t1(&g->getBoundingPoint(0), &g->getBoundingPoint(1*mult), &g->getBoundingPoint(3*mult)) ;
// 				TriPoint t2(&g->getBoundingPoint(0), &g->getBoundingPoint(2*mult), &g->getBoundingPoint(3*mult)) ;
// 				TriPoint t3(&g->getBoundingPoint(1*mult), &g->getBoundingPoint(2*mult), &g->getBoundingPoint(3*mult)) ;
//
// 				Plane p0(*t0.point[0], t0.normal) ;
// 				Plane p1(*t1.point[0], t1.normal) ;
// 				Plane p2(*t2.point[0], t2.normal) ;
// 				Plane p3(*t3.point[0], t3.normal) ;

        }

        return ret ;
    }
    case TETRAHEDRON:
    {
        return g->intersection(this) ;
    }
    case HEXAHEDRON:
    {

        if(g->getGeometryType() == SPHERE)
            return g->intersection(this) ;
    }
    case LOFTED_POLYGON:
    {
        if(g->getGeometryType() == LOFTED_POLYGON)
        {
            // first, find an intersection point, assuming two lofted cylinders
//             double selfBaseDiameter = dynamic_cast<const LoftedPolygonPrism *>(this)->base.getRadius() ;
//             double gBaseDiameter = dynamic_cast<const LoftedPolygonPrism *>(g)->base.getRadius() ;
//             Point start = (dynamic_cast<const LoftedPolygonPrism *>(this)->interpolationPoints[0]+dynamic_cast<const LoftedPolygonPrism *>(g)->interpolationPoints[0])*0.5 ;
        }
    }
    default:
    {
        return ret ;
    }
    }

}


ConvexGeometry::ConvexGeometry(): ConvexPolygon((size_t)0)
{
}

ConvexGeometry::ConvexGeometry(size_t s) : ConvexPolygon(s)
{
}

const Point & Geometry::getCenter() const
{
    return this->center;
}

size_t & Geometry::timePlanes()
{
    return this->time_planes ;
}

size_t Geometry::timePlanes() const
{
    return this->time_planes ;
}

Point & Geometry::getCenter()
{
    return this->center;
}

PointSet::PointSet(size_t npoints)
{
    this->boundingPoints.resize(npoints) ;
    for(size_t i = 0 ; i < npoints ; i++)
        boundingPoints[i] = nullptr ;
    this->chullEndPos = 0;
} 

size_t ConvexGeometry::size() const
{
    return 	boundingPoints.size() + inPoints.size() ;
}

size_t NonConvexGeometry::size() const
{
    return 	boundingPoints.size() + inPoints.size() ;
}

double PointSet::x(size_t i)
{
    return boundingPoints[i]->getX() ;
}

double PointSet::y(size_t i)
{
    return boundingPoints[i]->getY() ;
}

double PointSet::z(size_t i)
{
    return boundingPoints[i]->getZ() ;
}

void PointSet::setX(size_t i, double vv)
{
    boundingPoints[i]->Point::setX(vv) ;
}

void PointSet::setY(size_t i, double vv)
{
    boundingPoints[i]->Point::setY(vv) ;
}

void PointSet::setZ(size_t i, double vv)
{
    boundingPoints[i]->Point::setZ(vv) ;
}

void PointSet::set(size_t i, Point * p)
{
    delete boundingPoints[i] ;
    boundingPoints[i] = p;
}

void PointSet::set(size_t i, double x, double y)
{
    int id = -1 ;
    if(boundingPoints[i] != nullptr)
        id = boundingPoints[i]->getId() ;

    delete boundingPoints[i] ;

    boundingPoints[i] = new Point(x,y) ;
    boundingPoints[i]->setId( id ) ;
}

void PointSet::set(size_t i, double x, double y, double z)
{
    int id = -1 ;
    if(boundingPoints[i] != nullptr)
        id = boundingPoints[i]->getId() ;

    delete boundingPoints[i] ;

    boundingPoints[i] = new Point(x,y,z) ;
    boundingPoints[i]->setId( id ) ;
}

Point * PointSet::getPoint(size_t i) const
{
    return boundingPoints[i] ;
}

Point * PointSet::getPoint(size_t i)
{
    return boundingPoints[i] ;
}

PointSet::const_iterator PointSet::begin() const
{
    return &boundingPoints[0] ;
}

PointSet::const_iterator PointSet::end() const
{
    return  &boundingPoints[boundingPoints.size()] ;
}

PointSet::iterator PointSet::begin()
{
    return &boundingPoints[0] ;
}

PointSet::iterator PointSet::end()
{
    return  &boundingPoints[boundingPoints.size()] ;
}

ConvexPolygon * PointSet::convexHull() const
{

    if(this->chullEndPos != 0 )
    {
        ConvexPolygon *hull = new ConvexPolygon(chullEndPos) ;
        std::copy(begin(), end(), hull->begin()) ;
        return hull ;
    }

    ConvexPolygon *ret = new ConvexPolygon(this) ;

    return ret ;

}

size_t PointSet::size() const
{
    return  boundingPoints.size() ;
}

Point PointSet::computeCenter() const
{
    Point ret;
    for (size_t i = 0 ; i < boundingPoints.size() ; i++)
    {
        ret.getX() += boundingPoints[i]->getX()/boundingPoints.size() ;
        ret.getY() += boundingPoints[i]->getY()/boundingPoints.size() ;
    }
    return ret ;
}

void PointSet::removePoint(size_t index)
{
    std::valarray<Point *> n(size()-1) ;
    for (size_t i = 0 ; i < index ; i++)
    {
        n[i] = boundingPoints[i] ;
    }
    for (size_t i = index+1 ; i < size() ; i++)
    {
        n[i] = boundingPoints[i] ;
    }

    boundingPoints.resize(n.size()) ;
    boundingPoints = n ;
}



Nurb::Nurb(std::vector<double> knot, size_t degree)
{
    this->knot = knot ;
    this->degree = degree;
}

Nurb::Nurb(std::vector<Point> controlPoint, std::vector<double> weight,std::vector<double> knot, size_t degree)
{
    this->controlPoint = controlPoint ;
    this->weight = weight ;
    this->knot = knot ;
    this->degree = degree ;
}

Nurb::Nurb(std::vector<Point> controlPoint, std::vector<double> weight, std::vector<double> knot)
{
    this->controlPoint = controlPoint ;
    this->weight = weight ;
    this->knot = knot ;

}

int Nurb::computeDegree(const std::vector<double> &vec,const std::vector<Point> &Pv)
{
    return vec.size() - Pv.size() - 1;
}

std::vector<double> Nurb::transVector(std::vector<double> &vec)
{
    double diff = vec.back() - vec.front();
    for(size_t i = 0 ; i != vec.size(); ++i)
    {
        vec[i] = ( vec[i] - vec.front() ) / diff ;
    }
    return vec;
}

double Nurb::getBasis(double u, int i, int k, const std::vector<double> &vec)
{
    if (k == 0)
    {
        if (vec[i] <= u && u <= vec[i+1])
        {
            return 1;
        }
        return 0;
    }
    double denom1 = vec[i+k]-vec[i];
    double denom2 = vec[i+k+1]-vec[i+1];
    double part1 = 0;
    double part2 = 0;
    if (denom1 > 0)
    {
        part1 = ((u-vec[i])/denom1)*getBasis(u,i,k-1,vec);
    }
    if (denom2 > 0)
    {
        part2 = ((vec[i+k+1]-u)/denom2)*getBasis(u,i+1,k-1,vec);
    }
    return part1 + part2;
}

Point Nurb::getNurbPoint(double u, double r)
{
    sort(knot.begin(),knot.end());
    transVector(knot);
    int k = computeDegree(knot,controlPoint);
    Point outPointNom;
    double outDenom = 0;
    for(size_t i = 0 ; i != controlPoint.size() ; ++i)
    {
        double koef = getBasis(u,i,k,knot);
        if(koef > r)
        {
            outPointNom += controlPoint[i]*weight[i]*koef ;
            outDenom += weight[i]*koef ;
        }
    }
    (outPointNom / outDenom).print();
    return outPointNom / outDenom;
}

std::vector<Point> Nurb::sampleNurb(size_t n)
{
    std::vector<Point>  nurbPoint;
    double step = 1./n;
    for(double u = 0; u != 1; u += step)
    {
        nurbPoint.push_back(getNurbPoint(u));
    }
    return nurbPoint;
}

void WeightedPoint::print() const
{
    std::cout << " ( id = " << id << " ; "<< x << ", " << y << ", "<< z << ", "<< w <<") " << std::endl ;
}

WeightedPoint WeightedPoint::operator*(const double & p)  const
{
    return WeightedPoint(x*p, y*p, z*p, w*p);
}

WeightedPoint WeightedPoint::operator+(const WeightedPoint &wp)  const
{
    return WeightedPoint(x + wp.getX(), y + wp.getY(),z + wp.getZ(), w + wp.w);
}

WeightedPoint WeightedPoint::operator/(const double p) const
{
    return WeightedPoint(x/p, y/p, z/p, w/p);
}

Point Nurb::pointOnNurb(double u)
{
    size_t i = 0;

    while (u >= knot[i])
        ++i;

    size_t leftBound = i - 1;

    size_t numberOfCycles, knotMultiplicity  = 0;

    if ( u == knot[leftBound] )
    {
        knotMultiplicity = count(knot.begin(),knot.end(),u);
        numberOfCycles = degree - knotMultiplicity;
    }
    else numberOfCycles = degree;

    size_t r = leftBound - knotMultiplicity + 1;

    std::valarray<WeightedPoint> originalWPoint(controlPoint.size());
    for (size_t i = 0; i != originalWPoint.size(); ++i)
    {
        originalWPoint[i] = WeightedPoint(controlPoint[i]*weight[i],weight[i]);
    }

    double coeff = 0;
    std::valarray<WeightedPoint> generatedWPoint(originalWPoint.size());
    for (size_t j = 1; j <= numberOfCycles; ++j)
    {
        for (size_t k = leftBound - degree + j; k != r ; ++k)
        {
            coeff = ( u - knot[k] ) / ( knot[k+1+degree-j] - knot[k] );
            generatedWPoint[k] = originalWPoint[k-1] * ( 1 - coeff) + originalWPoint[k] *coeff;
        }


        for (size_t k = 0; k != originalWPoint.size(); ++k)
        {
            originalWPoint[k]=generatedWPoint[k];
        }
    }

    Point resultPoint = originalWPoint[r - 1] / originalWPoint[r - 1].w;
    resultPoint.print();
    return resultPoint;
}

Geometry::Geometry(): inPoints(0),gType(nullptr_GEOMETRY)
{
    sampled = false ;
    time_planes = 1 ;
}

Geometry::Geometry(size_t numPoints):inPoints(numPoints), gType(nullptr_GEOMETRY)
{
    time_planes = 1 ;
}

GeometryType Geometry::getGeometryType() const
{
    return gType ;
}

ConvexPolygon::ConvexPolygon(size_t npoints) : PointSet(npoints)
{
}

ConvexPolygon::ConvexPolygon(const PointSet * po) : PointSet(po->size())
{
    std::vector<Point*> *points = new std::vector<Point*>(po->size()) ;
    std::copy(po->begin(), po->end(), points->begin()) ;

    Point *pivot( (*points)[0] ) ;
    for(size_t i = 1 ; i < points->size() ; i++)
    {
        if (pivot->getY() < (*points)[i]->getY())
            pivot = (*points)[i] ;
    }

    //! we build a map ordered by the angle to the pivot.

    std::map< double, Point *>  pointSet ;

    for(size_t i = 0 ; i < points->size() ; i++)
    {
        if ((*(*points)[i]) != (*pivot))
        {
            double angle = atan2(pivot->getY()-(*points)[i]->getY(), pivot->getX()-(*points)[i]->getX()) ;

            if(pointSet.find(angle) != pointSet.end())
            {
                if(squareDist((*pointSet[angle]), (*pivot)) < squareDist((*(*points)[i]), (*pivot)))
                {
                    pointSet[angle] = (*points)[i] ;
                }
            }
            else
                pointSet[angle] = (*points)[i] ;
        }
    }

    pointSet[-10] = pivot ;

    //!we build the convex hull : we are on the convex hull if and only if we only turned left.

    std::vector<Point *> temphull ;
    std::vector<Point *> orderedset ;

    for(auto i = pointSet.begin() ; i!= pointSet.end() ; ++i)
    {
        orderedset.push_back(i->second) ;
    }

    temphull.push_back(orderedset[0] ) ;
    temphull.push_back(orderedset[1] ) ;

    for(auto i = orderedset.begin()+2 ; i != orderedset.end() ; ++i)
    {
        for(size_t j = 0 ; j < temphull.size() ; j++)

            //! this is a usual cross product of two vectors...
            if(  ((*i)->getY() - (*temphull.rbegin())->getY())*((*temphull.rbegin())->getX() - temphull[temphull.size()-2]->getX() ) -
                 ((*i)->getX() - (*temphull.rbegin())->getX())*((*temphull.rbegin())->getY() - temphull[temphull.size()-2]->getY() ) >
                    POINT_TOLERANCE )
            {
                temphull.push_back(*i) ;
            }
            else
            {
                while( !(((*i)->getY() - (*temphull.rbegin())->getY())*((*temphull.rbegin())->getX() - temphull[temphull.size()-2]->getX() ) -
                         ((*i)->getX() - (*temphull.rbegin())->getX())*((*temphull.rbegin())->getY() - temphull[temphull.size()-2]->getY() ) >
                         POINT_TOLERANCE))
                {
                    temphull.pop_back();
                }
                temphull.push_back(*i) ;
            }


    }

    delete points ;
    std::copy(temphull.begin(), temphull.end(),this->begin()) ;
}

Segment::Segment(const Segment & l)
{
    f = l.first() ;
    s = l.second() ;
    vec = l.vector() ;
    mid = l.midPoint() ;
}

Segment& Segment::operator=(const Segment& l)
{
    f = l.first() ;
    s = l.second() ;
    vec = l.vector() ;
    mid = l.midPoint() ;
    return *this ;
}

bool ConvexPolygon::isTrigoOriented()  const
{

    for (size_t i = 0 ;  i <  boundingPoints.size()-1; i++)
    {
        Point v_0 = *boundingPoints[(i+1)%boundingPoints.size()]
                    - *boundingPoints[(i)%boundingPoints.size()] ;
        Point v_1 = *boundingPoints[(i+2)%boundingPoints.size()]
                    - *boundingPoints[(i+1)%boundingPoints.size()] ;
        if((v_0^v_1).getZ() < 0 )
            return false ;
    }

    return true ;

}

NonConvexGeometry::NonConvexGeometry() : PointSet(1)
{
    orderedSet.resize(1) ;
    orderedSet[0] = inPoints[0] ;
}

NonConvexGeometry::NonConvexGeometry(size_t numPoints) : PointSet(numPoints)
{
    orderedSet.resize(numPoints) ;
    for(size_t i = 0 ; i < numPoints ; i++)
        orderedSet[i] = boundingPoints[i] ;
}

NonConvexGeometry::NonConvexGeometry(const PointArray & p)
{
    boundingPoints.resize(p.size()) ;
    std::copy(&p[0], &p[p.size()], begin()) ;//&boundingPoints[0]) ;
}

const Point & NonConvexGeometry::getPoint(size_t i) const
{
    if (i < inPoints.size())
        return *inPoints[i] ;
    return *boundingPoints[i-inPoints.size()] ;
}

Point & NonConvexGeometry::getPoint(size_t i)
{
    if (i < inPoints.size())
        return *inPoints[i] ;
    return *boundingPoints[i-inPoints.size()] ;
}

Line::Line(const Point & origin, const Point &vector)
{
    p = origin ;
    v = vector ;
}

Line::Line(const Segment & base)
{
    p = base.first() ;
    v = base.vector() ;
}

Line::Line(const Line & l, const Point & through)
{
    p = through ;
    v = l.vector() ;
}

Line::Line()
{
    p = Point(0,0) ;
    v = Point(1,1) ;
}

bool Line::intersects(const Line &l) const
{
//    std::cout << v.getX() * l.vector().getY() - v.getY() * l.vector().getX() << std::endl ;
    return (std::abs(v.getX() * l.vector().getY() - v.getY() * l.vector().getX()) > POINT_TOLERANCE) ;
}

bool Line::intersects(const Segment &s) const
{

    if (std::abs(-v.getY() * s.vector().getX() + v.getX() * s.vector().getY()) <= std::numeric_limits<double>::epsilon())
        return false ;

    Matrix m(2,2) ;
    Vector vv(2) ;

    m[0][0] = s.vector().getX() ;
    m[0][1] = -v.getX() ;
    m[1][0] = s.vector().getY() ;
    m[1][1] = -v.getY() ;

    vv[0] = p.getX()-s.second().getX()  ;
    vv[1] = p.getY()-s.second().getY()   ;

    invert2x2Matrix(m) ;

    Vector fac = m * vv ;

    return fac[0] < 1 && fac[0] > 0 ;
}

bool Line::on(const Point &m) const
{
    if( m == p)
        return true ;
    return isAligned(m, p, (p+v)) ;
}

const Point & Line::vector() const
{
    return v ;
}

const Point & Line::origin() const
{
    return p ;
}

Point Line::intersection(const Line &l) const
{
    Matrix mxy(2,2) ;
    Matrix mxz(2,2) ;
    Matrix myz(2,2) ;
    Vector vec(2) ;

    mxy[0][0] = v.getX() ;
    mxy[0][1] = -l.vector().getX() ;
    mxy[1][0] = v.getY() ;
    mxy[1][1] = -l.vector().getY() ;

    if(std::abs(det(mxy)) > POINT_TOLERANCE)
    {
        vec[0] = l.origin().getX() - p.getX() ;
        vec[1] = l.origin().getY() - p.getY() ;
        vec[0] = l.origin().getX() - p.getX() ;
        vec[1] = l.origin().getY() - p.getY() ;

        invert2x2Matrix(mxy) ;

        Vector fac = mxy * vec ;
        return p + v*fac[0];
    }

    mxz[0][0] = v.getX() ;
    mxz[0][1] = -l.vector().getX() ;
    mxz[1][0] = v.getZ() ;
    mxz[1][1] = -l.vector().getZ() ;

    if(std::abs(det(mxz)) > POINT_TOLERANCE)
    {
        vec[0] = l.origin().getX() - p.getX() ;
        vec[1] = l.origin().getZ() - p.getZ() ;
        vec[0] = l.origin().getX() - p.getX() ;
        vec[1] = l.origin().getZ() - p.getZ() ;

        invert2x2Matrix(mxz) ;

        Vector fac = mxz * vec ;
        return p + v*fac[0];
    }

    myz[0][0] = v.getY() ;
    myz[0][1] = -l.vector().getY() ;
    myz[1][0] = v.getZ() ;
    myz[1][1] = -l.vector().getZ() ;

    if(std::abs(det(myz)) > POINT_TOLERANCE)
    {
        vec[0] = l.origin().getY() - p.getY() ;
        vec[1] = l.origin().getZ() - p.getZ() ;
        vec[0] = l.origin().getY() - p.getY() ;
        vec[1] = l.origin().getZ() - p.getZ() ;

        invert2x2Matrix(myz) ;

        Vector fac = myz * vec ;
        return p + v*fac[0];
    }


    return Point() ;

// 	vec[0] = l.origin().getX() - p.getX() ; vec[1] = l.origin().getY() - p.getY() ;
// 	m.print();
// 	invert3x3Matrix(m) ;
// 	m.print();
// 	Vector fac = m * vec ;
// 	return p + v*fac[0];

}

Point Line::intersection(const Segment &s) const
{
    Matrix m(2,2) ;
    Vector vec(2) ;

    m[0][0] = v.getX() ;
    m[0][1] = -s.vector().getX() ;
    m[1][0] = v.getY() ;
    m[1][1] = -s.vector().getY() ;

    vec[0] = s.first().getX() - p.getX() ;
    vec[1] = s.first().getY() - p.getY() ;

    invert2x2Matrix(m) ;

    Vector fac = m * vec ;

    return p + v*fac[0];
}


bool Line::intersects(const TriPoint &g) const
{
    if(isCoplanar(&g.point[1],&g.point[0], &g.point[2],&p))
        return true ;

    Matrix mat(3,3) ;

    mat[0][0] = v.getX();
    mat[0][1] = g.point[1].getX()-g.point[0].getX();
    mat[0][2] = g.point[2].getX()-g.point[0].getX();
    mat[1][0] = v.getY();
    mat[1][1] = g.point[1].getY()-g.point[0].getY();
    mat[1][2] = g.point[2].getY()-g.point[0].getY();
    mat[2][0] = v.getZ();
    mat[2][1] = g.point[1].getZ()-g.point[0].getZ();
    mat[2][2] = g.point[2].getZ()-g.point[0].getZ();
    return abs(det(mat)) > POINT_TOLERANCE ;


}

std::vector<Point> Line::intersection(const TriPoint &s) const
{
    Point u = s.point[1]-s.point[0] ;
    Point vec = s.point[2]-s.point[0] ;
    Point wp = p-s.point[0] ;
    double a = s.normal*wp ;
    double b = s.normal*v ;
    
    if(std::abs(b) < POINT_TOLERANCE)
    {
        if(abs(a) > POINT_TOLERANCE)
            return std::vector<Point>(0) ;
        
        Triangle t(s.point[0], s.point[1], s.point[2]) ;
        return intersection(&t) ;
    }
    
    double r = a/b ;
    
    std::vector<Point> ret ;
    Point inter0 = p + r*v ;
    Point inter1 = p - r*v ;
    double uu, uv, vv, wu0, wv0, wu1, wv1, D;
    uu = u*u;
    uv = u*vec;
    vv = vec*vec;
    Point w0 = inter0 -s.point[0];
    wu0 = w0*u;
    wv0 = w0*vec;
    
    Point w1 = inter1 -s.point[0];
    wu1 = w1*u;
    wv1 = w1*vec;
    
    D = uv * uv - uu * vv;
    double ss0, ss1, t0, t1 ;
    ss0 = (uv * wv0 - vv * wu0) / D;
    bool in0 = true ;
    if (ss0 < 0.0 || ss0 > 1.0)
         in0 = false ;
    t0 = (uv * wu0 - uu * wv0) / D;
    if (t0 < 0.0 || (ss0 + t0) > 1.0)
         in0 = false ;
    
    ss1 = (uv * wv1 - vv * wu1) / D;
    bool in1 = true ;
    if (ss1 < -POINT_TOLERANCE  || ss1 > 1.0-POINT_TOLERANCE)
         in1 = false ;
    t1 = (uv * wu1 - uu * wv1) / D;
    if (t1 < -POINT_TOLERANCE || (ss1 + t1) > 1.0-POINT_TOLERANCE)
         in1 = false ;
    
    if(in0)
        ret.push_back(inter0);
    if(in1)
        ret.push_back(inter1);
    
    return ret ;
    
//     Point vec = s.point[2]-s.point[0] ;
//     double a = -(s.normal*(p-s.point[0])) ;
//     double b = vec*s.normal ;
// 
//     if(std::abs(b) < POINT_TOLERANCE)
//     {
//         if(abs(a) > POINT_TOLERANCE)
//             return std::vector<Point>(0) ;
// 
//         Triangle t(s.point[0], s.point[1], s.point[2]) ;
//         return intersection(&t) ;
//     }
// 
//     double r = a/b ;
// 
//     std::vector<Point> ret ;
//     Point p0 = p + v*r ;
//     Point p1 = p - v*r ;
//     if(s.in(p0))
//         ret.push_back( p0) ;
//     if(s.in(p1))
//         ret.push_back( p1) ;
//     return ret ;
}

bool Line::intersects(const Geometry *g) const
{
    switch(g->getGeometryType())
    {
    case CIRCLE:
    {
        return squareDist(projection(g->getCenter()), g->getCenter()) < g->getRadius()*g->getRadius() ;
    }
    case ELLIPSE:
    {
        // see Geometry->intersects() [same method as ellipse-ellipse intersection]
        Point onEllipse(g->getCenter()) ;
        Point onLine = this->projection(onEllipse) ;
        if(dynamic_cast<const Ellipse*>(g)->in(onLine))
            return true ;
        if((onEllipse-onLine).norm() > dynamic_cast<const Ellipse*>(g)->getMajorRadius()*1.1)
            return false ;
        Point onEllipse2(onEllipse) ;
        onEllipse = dynamic_cast<const Ellipse*>(g)->project(onLine) ;
        int nnn = 0 ;
        while((onEllipse-onEllipse2).norm() > POINT_TOLERANCE)
        {
            onEllipse2 = onEllipse ;
            onLine = this->projection(onEllipse) ;
            onEllipse = dynamic_cast<const Ellipse*>(g)->project(onLine) ;
            nnn++ ;
        }

        return (onLine-onEllipse).norm() < POINT_TOLERANCE*100. ;
    }
    case TRIANGLE:
    {
        bool ret = false ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            ret = ret || s.intersects(*this) ;
        }

        return ret ;
    }
    case RECTANGLE:
    {

        bool ret = false ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            ret = ret || s.intersects(*this) ;
        }

        return ret ;
    }
    case CONVEX_POLYGON:
    {

        bool ret = false ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            ret = ret || s.intersects(*this) ;
        }

        return ret ;
    }
    case SEGMENTED_LINE:
    {

        bool ret = false ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size()-1 ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
            ret = ret || s.intersects(*this) ;
        }

        return ret ;
    }
    case SPHERE:
    {
        double a = v.sqNorm() ;
        double b = ((p.getX()-g->getCenter().getX())*v.getX() + (p.getY()-g->getCenter().getY())*v.getY() + (p.getZ()-g->getCenter().getZ())*v.getZ())*2. ;
        double c = (p.getX()-g->getCenter().getX())*(p.getX()-g->getCenter().getX())
                   +(p.getY()-g->getCenter().getY())*(p.getY()-g->getCenter().getY())
                   +(p.getZ()-g->getCenter().getZ())*(p.getZ()-g->getCenter().getZ())
                   -g->getRadius()*g->getRadius() ;
        double delta = b*b - 4.*a*c ;

        return delta >= 0 ;
    }
    default:
    {
        std::cout << "case not solved : " << g->getGeometryType() << std::endl ;
        return false ;
    }
    }
}


std::vector<Point> Line::intersection(const Geometry * g) const
{
    switch(g->getGeometryType())
    {
    case CIRCLE:
    {

        double a = v.sqNorm() ;
        double b = ((p.getX()-g->getCenter().getX())*v.getX() + (p.getY()-g->getCenter().getY())*v.getY())*2. ;
        double c = (p.getX()-g->getCenter().getX())*(p.getX()-g->getCenter().getX())
                   +(p.getY()-g->getCenter().getY())*(p.getY()-g->getCenter().getY())
                   -g->getRadius()*g->getRadius() ;
        double delta = b*b - 4.*a*c ;

        if (delta > POINT_TOLERANCE)
        {
            std::vector<Point> ret ;
            ret.push_back(p+v*((-b + sqrt(delta))/(2.*a))) ;
            ret.push_back(p+v*((-b - sqrt(delta))/(2.*a))) ;
            return ret ;
        }
        else if (delta > 0)
        {

            std::vector<Point> ret ;
            ret.push_back(p+v*(-b/(2.*a))) ;
            return ret ;
        }


        return std::vector<Point>(0) ;

    }
    case ELLIPSE:
    {
        std::vector<Point> ret ;
        if(!this->intersects(g))
            return ret ;

        // get first point by iterative projections
        Point onEllipse(g->getCenter()) ;
        Point onLine = this->projection(onEllipse) ;
        Point onEllipse2(onEllipse) ;
        onEllipse = dynamic_cast<const Ellipse*>(g)->project(onLine) ;
        int nnn = 0 ;
        while((onEllipse-onEllipse2).norm() > POINT_TOLERANCE)
        {
            onEllipse2 = onEllipse ;
            onLine = this->projection(onEllipse) ;
            onEllipse = dynamic_cast<const Ellipse*>(g)->project(onLine) ;
            nnn++ ;
        }

        ret.push_back(onEllipse) ;

        // get direction; where to look next
        // the goal is to get a point on the line inside the ellipse
        double xfirst = (onEllipse - this->origin())*this->vector() ;
        double dx = 0.1 ;
        double xhere = xfirst ;
        double xprev = xhere -dx ;
        double xnext = xhere+dx ;
        bool dir = false ;
        while(!dir)
        {
            xprev = xhere-dx ;
            Point prev = this->origin() + this->vector()*xprev ;
            if(dynamic_cast<const Ellipse*>(g)->in(prev))
            {
                dir = true ;
                dx = -dx ;
            }
            xnext = xhere+dx ;
            Point next = this->origin() + this->vector()*xnext ;
            if(dynamic_cast<const Ellipse*>(g)->in(next))
            {
                dir = true ;
            }
            dx = dx/10. ;
            if(dx < POINT_TOLERANCE * 100.) {
                // line is tangent
                return ret ;
            }
        }

        // then ye get a point on a line outside the ellipse
        xnext = xhere+dx ;
        Point next = this->origin() + this->vector()*xnext ;
        while(dynamic_cast<const Ellipse*>(g)->in(next))
        {
            next += this->vector()*dx ;
        }

        // another run of iterative projections to get the second intersection
        onEllipse = (g->getCenter()) ;
        onLine = next ;
        onEllipse2 = (onEllipse) ;
        onEllipse = dynamic_cast<const Ellipse*>(g)->project(onLine) ;
        nnn = 0 ;
        while((onEllipse-onEllipse2).norm() > POINT_TOLERANCE)
        {
            onEllipse2 = onEllipse ;
            onLine = this->projection(onEllipse) ;
            onEllipse = dynamic_cast<const Ellipse*>(g)->project(onLine) ;
            nnn++ ;
        }

        ret.push_back(onEllipse) ;
        return ret ;

    }
    case TRIANGLE:
    {

        std::vector<Point> ret ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            if(s.intersects(*this))
            {
                Point inter =  s.intersection(*this) ;
                ret.push_back(inter) ;
            }
        }

        return ret ;
    }
    case RECTANGLE:
    {

        std::vector<Point> ret ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            if(s.intersects(*this))
                ret.push_back(intersection(s)) ;
        }
        return ret ;
    }
    case CONVEX_POLYGON:
    {

        std::vector<Point> ret ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            if(s.intersects(*this))
                ret.push_back(intersection(s)) ;
        }

        return ret ;
    }
    case SEGMENTED_LINE:
    {

        std::vector<Point> ret ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size()-1 ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
            if(s.intersects(*this))
                ret.push_back(intersection(s)) ;
        }

        return ret ;
    }
    case SPHERE:
    {
        double a = v.sqNorm() ;
        double b = ((p.getX()-g->getCenter().getX())*v.getX() + (p.getY()-g->getCenter().getY())*v.getY() + (p.getZ()-g->getCenter().getZ())*v.getZ())*2. ;
        double c = (p.getX()-g->getCenter().getX())*(p.getX()-g->getCenter().getX())
                   +(p.getY()-g->getCenter().getY())*(p.getY()-g->getCenter().getY())
                   +(p.getZ()-g->getCenter().getZ())*(p.getZ()-g->getCenter().getZ())
                   -g->getRadius()*g->getRadius() ;
        double delta = b*b - 4*a*c ;

        if(std::abs(delta) < POINT_TOLERANCE)
        {
            std::vector<Point> ret ;
            ret.push_back(p+v*(-b/(2.*a))) ;
            return ret ;
        }
        else if (delta >= POINT_TOLERANCE)
        {
            std::vector<Point> ret ;
            ret.push_back(p+v*((-b + sqrt(delta))/(2.*a))) ;
            ret.push_back(p+v*((-b - sqrt(delta))/(2.*a))) ;
            return ret ;
        }

        return std::vector<Point>(0) ;

    }
    default:
    {
        std::cout << "line-geometry type " << g->getGeometryType() << " not implemented" << std::endl ;
        return std::vector<Point>() ;
    }
    }
}

Point Line::projection(const Point &m ) const
{
    if(on(m))
        return m ;

    Point n0 = v^(m-p) ;
    n0 /= n0.norm() ;
    Point n = v^n0 ;
    n /= n.norm() ;
    Line line(m, n) ;

    return line.intersection(*this) ;
}

Point Segment::normal() const
{
    double n = norm();
    if(n > POINT_TOLERANCE)
        return Point(-vec.getY(), vec.getX())/n ;

    return Point(0, 0) ;
}

Point Segment::normal(const Point & inside) const
{
    bool sameSide = (mid-inside)*normal() > 0 ;
    double sign = 1. ;
    if(!sameSide)
        sign = -1 ;

    double n = norm();
    if(n > POINT_TOLERANCE)
        return Point(-vec.getY()*sign, vec.getX()*sign)/n ;

    return Point(0, 0) ;

}

Vector Segment::normalv(const Point & p) const
{
    bool sameSide = (mid-p)*normal() > 0 ;
    double sign = 1. ;
    if(!sameSide)
        sign = -1 ;
    double n = vec.norm();
    if(n > POINT_TOLERANCE)
    {
        Vector ret(2) ;
        ret[0] = -vec.getY()/n*sign ;
        ret[1] = vec.getX()/n*sign ;
        return ret ;
    }

    return Vector(0., 2) ;
}

Segment::Segment(const Point & p0, const Point & p1)
{

    f = p0 ;
    s = p1 ;
    mid = (p0+p1)*0.5;
    vec = f-s ;
}

double Segment::norm() const
{
    return vec.norm() ;
}

std::valarray<std::pair<Point, double> > Segment::getGaussPoints(bool timeDependent) const
{
    std::valarray< std::pair<Point, double> > gp(2+2*timeDependent) ;
    Point a( f.getX()*0.788675134594813+ s.getX()*(1.-0.788675134594813),f.getY()*0.788675134594813+ s.getY()*(1.-0.788675134594813)) ;
    Point b( s.getX()*0.788675134594813+ f.getX()*(1.-0.788675134594813),s.getY()*0.788675134594813+ f.getY()*(1.-0.788675134594813)) ;
//	double n = norm() ;
    if(!timeDependent)
    {
        gp[0] = std::pair<Point, double>(a, 1) ;
        gp[1] = std::pair<Point, double>(b, 1) ;
    }
    else
    {
        gp[0] = std::pair<Point, double>(Point(a.getX(), a.getY(), a.getZ(), -0.577350269189626), 1) ;
        gp[1] = std::pair<Point, double>(Point(b.getX(), b.getY(), b.getZ(), -0.577350269189626), 1) ;
        gp[2] = std::pair<Point, double>(Point(a.getX(), a.getY(), a.getZ(), 0.577350269189626), 1) ;
        gp[3] = std::pair<Point, double>(Point(b.getX(), b.getY(), b.getZ(), 0.577350269189626), 1) ;
    }
    return gp ;
}

Segment::Segment()
{
    f = Point(0,0) ;
    s = Point(0,0) ;
    mid = f*0.5 + s*0.5;
    vec = f-s ;
}
Segment::~Segment()
{
}


void Segment::print() const
{
    std::cout << "[ (" << f.getX() << ", " << f.getY() << ") ; (" << s.getX() << ", " << s.getY() << ") ]" << std::endl ;
}

bool Segment::intersects(const Line & l) const
{
    if (std::abs(-vec.getX()*l.vector().getY()+l.vector().getX()*vec.getY()) <= std::numeric_limits<double>::epsilon())
        return false ;

    Matrix m(2,2) ;
    Vector vv(2) ;

    m[0][0] = vec.getX() ;
    m[0][1] = -l.vector().getX() ;
    m[1][0] = vec.getY() ;
    m[1][1] = -l.vector().getY() ;

    vv[0] = l.origin().getX() - s.getX() ;
    vv[1] = l.origin().getY() - s.getY() ;

    invert2x2Matrix(m) ;

    Vector fac = m * vv ;
    return fac[0] < 1 && fac[0] > 0 ;

}

bool Segment::intersects(const Plane & l) const
{
    return l.intersects(*this) ;

}

bool Segment::intersects(const Geometry *g) const
{
    switch(g->getGeometryType())
    {
    case TIME_DEPENDENT_CIRCLE:
    {
// 			s.print() ;
// 			f.print() ;



        if(g->in(s) && g->in(f))
            return false ;
        if(g->in(s) || g->in(f))
            return true ;

        Point center(g->getCenter()) ;
        center.getT() = s.getT() ;
        Point proj = project(center) ;
        if(g->in(proj)) {
            return true ;
        }
        center = g->getCenter() ;
        center.getT() = f.getT() ;
        proj = project(center) ;
        if(g->in(proj)) {
            return true ;
        }
        return false ;
    }
    case CIRCLE:
    {
        if(g->in(s) && g->in(f))
            return false ;
        if(g->in(s) || g->in(f))
            return true ;

        Point center(g->getCenter()) ;
        Point proj = project(center) ;
        return g->in(proj) ;
    }
    case ELLIPSE:
    {
/*	Point pfirst = dynamic_cast<const Ellipse *>(g)->toLocalCoordinates(f) ;
	Point psecond = dynamic_cast<const Ellipse *>(g)->toLocalCoordinates(s) ;
	Segment inLocalCoordinates(pfirst, psecond) ;
	Circle c(0.,0.,1.) ;
	return inLocalCoordinates.intersects(&c) ;*/

        Ellipse ell(g->getCenter(), dynamic_cast<const Ellipse *>(g)->getMajorAxis(), dynamic_cast<const Ellipse *>(g)->getMinorAxis()) ;

			ell.sampleBoundingSurface(128) ;

        for(size_t i = 0 ; i < ell.getBoundingPoints().size()-1 ; i++)
        {
            Segment test(ell.getBoundingPoint(i), ell.getBoundingPoint(i+1)) ;
            if(this->intersects(test))
                return true ;
        }
        Segment test(ell.getBoundingPoint(ell.getBoundingPoints().size()-1), ell.getBoundingPoint(0)) ;
        return this->intersects(test) ;


        Line l(f,s-f) ;


        std::vector<Point> in = l.intersection(g) ;
        if(in.size() == 0)
            return false ;
        for(size_t i = 0 ; i < in.size() ; i++)
        {
            if(this->on(in[i]))
                return true ;
        }
        return false ;
    }
    case TRIANGLE:
    {
        std::vector<Point> pts ;
        std::multimap<double, Point> pt ;
        size_t planes = 1 ;
        if(g->timePlanes() > 1)
            planes = g->timePlanes() ;
        if( g->getBoundingPoints().size() != 3*planes)
        {
            for(size_t i = 0 ; i < g->getBoundingPoints().size() ;  i++)
            {
                pt.insert(std::make_pair(
                              std::abs(
                                  squareDist2D(dynamic_cast<const Triangle *>(g)->getCircumCenter(), g->getBoundingPoint(i))-g->getRadius()*g->getRadius()), g->getBoundingPoint(i)));
            }
            int count  = 0 ;
            for(auto i = pt.begin() ; count < 3 ; i++ )
            {
                count++ ;
                pts.push_back(i->second);
            }
        }
        else
        {
            pts.push_back(g->getBoundingPoint(0));
            pts.push_back(g->getBoundingPoint(1));
            pts.push_back(g->getBoundingPoint(2));
        }


        if(this->on(pts[0]) || this->on(pts[1]) || this->on(pts[2]))
            return true ;

        Segment sa(pts[0],pts[1]) ;
        Segment sb(pts[1],pts[2]) ;
        Segment sc(pts[2],pts[0]) ;


        return sa.intersects(*this) || sb.intersects(*this) || sc.intersects(*this) ;

    }
    case RECTANGLE:
    {
        std::vector<Point> box = g->getBoundingBox() ;
        Segment s0(box[0], box[1]) ;
        Segment s1(box[1], box[2]) ;
        Segment s2(box[2], box[3]) ;
        Segment s3(box[3], box[0]) ;

        return intersects(s0) || intersects(s1) || intersects(s2) || intersects(s3) ;
    }
    case PARALLELOGRAMME:
    {
        size_t n = g->getBoundingPoints().size() ;

        std::vector<Point> box ;
        box.push_back(g->getBoundingPoint(0)) ;
        box.push_back(g->getBoundingPoint(n*1/4)) ;
        box.push_back(g->getBoundingPoint(n*2/4)) ;
        box.push_back(g->getBoundingPoint(n*3/4)) ;
        Segment s0(box[0], box[1]) ;
        Segment s1(box[1], box[2]) ;
        Segment s2(box[2], box[3]) ;
        Segment s3(box[3], box[0]) ;

        return intersects(s0) || intersects(s1) || intersects(s2) || intersects(s3) ;
    }
    case SEGMENTED_LINE :
    {
        bool ret = false ;

        for(size_t i = 0 ; i <  g->getBoundingPoints().size()-1 ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
            ret = ret || s.intersects(*this) ;
        }

        return ret ;
    }
    case POLYGON:
    {

        bool ret = false ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            ret = ret || s.intersects(*this) ;
        }

        return ret ;
    }
    case CONVEX_POLYGON:
    {

        bool ret = false ;
        for(size_t i = 0 ; i <  g->getBoundingPoints().size() ;  i++)
        {
            Segment s(g->getBoundingPoint(i), g->getBoundingPoint((i+1)%g->getBoundingPoints().size())) ;
            ret = ret || s.intersects(*this) ;
        }

        return ret ;
    }
    case SPHERE:
    {
        if(g->in(s) && g->in(f))
            return false ;
        if(g->in(s) || g->in(f))
            return true ;
        Point center(g->getCenter()) ;
        Point v = s-f ;
        double uc = v*g->getCenter() ;
        double delta = uc*uc - v.sqNorm()*g->getCenter().sqNorm() - g->getRadius()*g->getRadius() ;
        return delta >= 0 ;
    }
    case POLYGON_PRISM:
    {
            if(g->in(s) == !g->in(f))
            return true ;
        
        
        Point dir = first() - second();
        int countin =0 ;
        int pointcount = 0 ;
        for(double i = 0 ; i < 1 ; i+=.2 )
        {
            Point toTest = second()+i*dir ;
            pointcount++;
            countin += g->in(toTest) ;
        }
        
        if(countin && countin != pointcount)
            return true ;
        return false ;   
    }
    case LOFTED_POLYGON:
    {
        if(g->in(s) == !g->in(f))
            return true ;
        
        Point dir = first() - second();
        int countin =0 ;
        int pointcount = 0 ;
        for(double i = 0 ; i < 1 ; i+=.2 )
        {
            Point toTest = second()+i*dir ;
            pointcount++;
            countin += g->in(toTest) ;
        }
        
        if(countin && countin != pointcount)
            return true ;
        return false ;
    }
    default:
        return false ;
    }
}

bool Segment::intersects(const TriPoint *g) const
{

    Point v(g->point[1] - g->point[0]) ;
    Point u(g->point[2] - g->point[0]) ;

    Point dir = first() - second();
    Point w0 = second() - g->point[0];
    Point n = u ^ v ;
    double a = -(n*w0) ;
    double b = n*dir;
    if (std::abs(b) < POINT_TOLERANCE)
    {
        if (std::abs(a) < POINT_TOLERANCE)
            return true;
        return false ;
    }

    double r = a / b ;

    if(r < 0 || r > 1)
        return false ;

    Point intersect(second()+dir*r) ;
    Point w(intersect-g->point[0]) ;

    double uv = u*v;
    double wu = u*w;
    double wv = v*w;
    double uu = u*u;
    double vv = v*v;
    double d = uv*uv - uu*vv ;
    double s = (uv*wv - vv*wu)/d ;
    if (s < -POINT_TOLERANCE || s > 1.+POINT_TOLERANCE)        // I is outside T
        return false;
    double t = (uv*wu - uu*wv)/d ;
    if (t < -POINT_TOLERANCE || (s + t) > 1.+POINT_TOLERANCE)
        return false ;

    return true ;

// 	if(isCoplanar(f, *g->point[0],*g->point[1],*g->point[2]))
// 		return g->in(f) ;
// 	if(isCoplanar(s, *g->point[0],*g->point[1],*g->point[2]))
// 		return g->in(s) ;
// 	Vector vec(3) ;
// 	vec[0] = f.getX() - g->point[0]->getX() ;
// 	vec[1] = f.getY() - g->point[0]->getY() ;
// 	vec[2] = f.getZ() - g->point[0]->getZ() ;
//
// 	Matrix mat(3,3) ;
// 	mat[0][0] = f.getX()-s.getX(); mat[0][1] = g->point[1]->getX()-g->point[0]->getX(); mat[0][2] = g->point[2]->getX()-g->point[0]->getX();
// 	mat[1][0] = f.getY()-s.getY(); mat[1][1] = g->point[1]->getY()-g->point[0]->getY(); mat[1][2] = g->point[2]->getY()-g->point[0]->getY();
// 	mat[2][0] = f.getZ()-s.getZ(); mat[2][1] = g->point[1]->getZ()-g->point[0]->getZ(); mat[2][2] = g->point[2]->getZ()-g->point[0]->getZ();
// 	Vector tuv = inverse3x3Matrix(mat)*vec ;
//
// 	return tuv.max() <=1 && tuv.min() >= 0 &&tuv[1] +tuv[2] <= 1  ;

}

std::vector<Point> Segment::intersection(const Geometry *g) const
{
    switch(g->getGeometryType())
    {
    case TRIANGLE:
    {
        std::vector<Point> ret ;

        std::vector<Point> pts ;
        const Triangle * gt = dynamic_cast<const Triangle *>(g) ;

        for(size_t i = 0 ; i <  gt->getBoundingPoints().size()/gt->timePlanes() ;  i++)
        {
            if(std::abs(dist(gt->getCircumCenter(), gt->getBoundingPoint(i))-gt->getRadius()) < POINT_TOLERANCE)
                pts.push_back(gt->getBoundingPoint(i));
        }

        for(size_t i = 0 ; i <  pts.size() ;  i++)
        {
            Segment s(pts[i], pts[(i+1)%pts.size()]) ;
            if(s.intersects(*this))
            {
                if(on(pts[i]))
                    ret.push_back(pts[i]) ;
                else
                    ret.push_back( s.intersection(*this)) ;
            }
        }
        std::sort(ret.begin(), ret.end()) ;
        auto e = std::unique(ret.begin(), ret.end()) ;
        ret.erase(e, ret.end()) ;
        return ret ;
    }
    case RECTANGLE:
    {
        std::vector<Point> ret ;
        std::vector<Point> bbox = g->getBoundingBox() ;
        Segment s0(bbox[0], bbox[1]) ;
        if(s0.intersects(*this))
        {
            if(isAligned(f, s0.first(), s0.second()) && isAligned(s, s0.first(), s0.second()))
            {
                if(s0.on(f)) {
                    ret.push_back(f) ;
                }
                if(s0.on(s)) {
                    ret.push_back(s) ;
                }
                if(on(s0.first())) {
                    ret.push_back(s0.first()) ;
                }
                if(on(s0.second())) {
                    ret.push_back(s0.second()) ;
                }
            }
            else {
                ret.push_back( s0.intersection(*this)) ;
            }

        }
        Segment s1(bbox[1], bbox[2]) ;

        if(s1.intersects(*this))
        {
            if(isAligned(f, s1.first(), s1.second()) && isAligned(s, s1.first(), s1.second()))
            {
                if(s1.on(f)) {
                    ret.push_back(f) ;
                }
                if(s1.on(s)) {
                    ret.push_back(s) ;
                }
                if(on(s1.first())) {
                    ret.push_back(s1.first()) ;
                }
                if(on(s1.second())) {
                    ret.push_back(s1.second()) ;
                }
            }
            else {
                ret.push_back( s1.intersection(*this)) ;
            }


        }

        Segment s2( bbox[2],  bbox[3]) ;
        if(s2.intersects(*this))
        {

            if(isAligned(f, s2.first(), s2.second()) && isAligned(s, s2.first(), s2.second()))
            {

                if(s2.on(f)) {
                    ret.push_back(f) ;
                }
                if(s2.on(s)) {
                    ret.push_back(s) ;
                }
                if(on(s2.first())) {
                    ret.push_back(s2.first()) ;
                }
                if(on(s2.second())) {
                    ret.push_back(s2.second()) ;
                }
            }
            else {
                ret.push_back( s2.intersection(*this)) ;
            }

        }

        Segment s3(bbox[3],bbox[0]) ;
        if(s3.intersects(*this))
        {

            if(isAligned(f, s3.first(), s3.second()) && isAligned(s, s3.first(), s3.second()))
            {

                if(s3.on(f)) {
                    ret.push_back(f) ;
                }
                if(s3.on(s)) {
                    ret.push_back(s) ;
                }
                if(on(s3.first())) {
                    ret.push_back(s3.first()) ;
                }
                if(on(s3.second())) {
                    ret.push_back(s3.second()) ;
                }
            }
            else {
                ret.push_back( s3.intersection(*this)) ;
            }

        }
        bool haveDuplicates = true ;
        while(haveDuplicates)
        {
            haveDuplicates = false ;
            for(size_t i  = 0 ; i < ret.size() ; i++)
            {
                for(size_t j  = i+1 ; j < ret.size() ; j++)
                {
                    if(squareDist3D(ret[i], ret[j])< 128*POINT_TOLERANCE*POINT_TOLERANCE)
                    {
                        haveDuplicates = true ;
                        ret.erase(ret.begin()+j) ;
                        break ;
                    }
                }

                if(haveDuplicates)
                    break ;
            }
        }
        return ret ;
    }
    case TIME_DEPENDENT_CIRCLE:
    {
        Circle c( dynamic_cast<const TimeDependentCircle *>(g)->radiusAtTime(s), g->getCenter().getX(), g->getCenter().getY()) ;
        std::vector<Point> tmp = this->intersection(&c) ;
        for(size_t i = 0 ; i < tmp.size() ; i++)
            tmp[i].getT() = s.getT() ;
        return tmp ;
    }
    case CIRCLE:
    {
        double a = vec.getX()*vec.getX() + vec.getY()*vec.getY() ;
        if(a < POINT_TOLERANCE)
            return std::vector<Point>(0) ;
        double dx = s.getX()-g->getCenter().getX() ;
        double dy = s.getY()-g->getCenter().getY() ;
        double b = 2.*(dx*vec.getX() + dy*vec.getY() );
        double c = dx*dx + dy*dy -g->getRadius()*g->getRadius() ;
        double delta = b*b - 4.*a*c ;

        if(std::abs(delta) < POINT_TOLERANCE)
        {
            std::vector<Point> ret ;
            Point A(s+vec*(-b/(2.*a))) ;
            if(on(A))
                ret.push_back(A) ;
            return ret ;
        }
        else if (delta >= POINT_TOLERANCE)
        {
            std::vector<Point> ret ;
            Point A(s+vec*(-b + sqrt(delta))/(2.*a)) ;
            if(on(A))
                ret.push_back(A) ;
            Point B(s+vec*(-b - sqrt(delta))/(2.*a)) ;
            if(on(B))
                ret.push_back(B) ;
            return ret ;
        }
        else
        {
            return std::vector<Point>(0) ;
        }
    }
    case ELLIPSE:
    {
        std::vector<Point> ret ;
//			Ellipse ell(g->getCenter(), dynamic_cast<const Ellipse *>(g)->getMajorAxis(), dynamic_cast<const Ellipse *>(g)->getMinorAxis()) ;

//			ell.sampleBoundingSurface(128) ;

        for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
        {
            Segment test(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
            if(this->intersects(test))
            {
                ret.push_back(this->intersection(test)) ;
            }
        }
        Segment test(g->getBoundingPoint(g->getBoundingPoints().size()-1), g->getBoundingPoint(0)) ;
        if(this->intersects(test))
        {
            ret.push_back(this->intersection(test)) ;
        }
        return ret ;
//			Segment test(ell.getBoundingPoint(ell.getBoundingPoints().size()-1), ell.getBoundingPoint(0)) ;
//			return this->intersects(test) ;


        /*			Line l(f,s-f) ;
        			std::vector<Point> ret = l.intersection(g) ;
        			for(size_t i = ret.size() ; i > 0 ; i--)
        			{
        				if(!(this->on(ret[i-1])))
        					ret.erase(ret.begin()+i-1) ;
        			}
        			return ret ;*/

    }
    case POLYGON:
    {
        std::vector<Point> ret ;
	std::valarray<Point> corners = dynamic_cast<const Polygon *>(g)->getOriginalPoints() ;
        for(size_t i = 0 ; i < corners.size()-1 ; i++)
        {
            Segment test(corners[i], corners[i+1]) ;
            if(test.intersects(*this))
            {
		if(dist( this->midPoint(), test.intersection(*this) ) > POINT_TOLERANCE)
                    ret.push_back(test.intersection(*this)) ;
                else
                    ret.push_back(test.first()) ;
            }
	   
        }
        Segment test(corners[corners.size()-1], corners[0]) ;
        if(test.intersects(*this))
        {
            if(dist( this->midPoint(), test.intersection(*this) ) > POINT_TOLERANCE)
                ret.push_back(test.intersection(*this)) ;
            else
                ret.push_back(test.first()) ;
        }
        return ret ;
    }
    case SEGMENTED_LINE:
    {
        std::vector<Point> ret ;
        for(size_t i = 0 ; i < g->getBoundingPoints().size()-1 ; i++)
        {
            Segment test(g->getBoundingPoint(i), g->getBoundingPoint(i+1)) ;
            if(test.intersects(*this))
                ret.push_back(test.intersection(*this)) ;
	   
        }
        return ret ;
    }
    case SPHERE:
    {
        double a = vec.getX()*vec.getX() + vec.getY()*vec.getY() + vec.getZ()*vec.getZ();
        double dx = s.getX()-g->getCenter().getX() ;
        double dy = s.getY()-g->getCenter().getY() ;
        double dz = s.getZ()-g->getCenter().getZ() ;
        double b = 2.*(dx*vec.getX() + dy*vec.getY() + dz*vec.getZ());
        double c = dx*dx + dy*dy + dz*dz-g->getRadius()*g->getRadius() ;
        double delta = b*b - 4.*a*c ;

        if(std::abs(delta) < POINT_TOLERANCE)
        {
            std::vector<Point> ret ;
            Point A(s+vec*(-b/(2.*a))) ;
            if(on(A))
                ret.push_back(A) ;
            return ret ;
        }
        else if (delta >= POINT_TOLERANCE)
        {
            std::vector<Point> ret ;
            Point A(s+vec*(-b + sqrt(delta))/(2.*a)) ;
            if(on(A))
                ret.push_back(A) ;
            Point B(s+vec*(-b - sqrt(delta))/(2.*a)) ;
            if(on(B))
                ret.push_back(B) ;
            return ret ;
        }
        else
        {
            return std::vector<Point>(0) ;
        }
    }
    case HEXAHEDRON:
    {
        std::vector<Point> ret ;
        std::vector<Point> bbox = g->getBoundingBox() ;
        double maxx =  bbox[0].getX() ;
        double minx =  bbox[0].getX() ;
        double maxy =  bbox[0].getY() ;
        double miny =  bbox[0].getY() ;
        double maxz =  bbox[0].getZ() ;
        double minz =  bbox[0].getZ() ;

        for(size_t i = 1 ; i < bbox.size() ; i++)
        {
            if(bbox[i].getX() > maxx)
                maxx = bbox[i].getX() ;
            if(bbox[i].getX() < minx)
                minx = bbox[i].getX() ;
            if(bbox[i].getY() > maxy)
                maxy = bbox[i].getY() ;
            if(bbox[i].getY() < miny)
                miny = bbox[i].getY() ;
            if(bbox[i].getZ() > maxz)
                maxz = bbox[i].getZ() ;
            if(bbox[i].getZ() < minz)
                minz = bbox[i].getZ() ;
        }

        Point corner1 (minx, miny, minz) ;
        Point corner2 (minx, miny, maxz) ;
        Point corner3 (minx, maxy, minz) ;
        Point corner4 (minx, maxy, maxz) ;
        Point corner5 (maxx, miny, minz) ;
        Point corner6 (maxx, miny, maxz) ;
        Point corner7 (maxx, maxy, minz) ;
        Point corner8 (maxx, maxy, maxz) ;

        Plane plane0(corner1,corner2,corner3) ;
        if(plane0.intersects(*this))
            ret.push_back(plane0.intersection(*this));
        Plane plane1(corner2,corner3,corner4) ;
        if(plane1.intersects(*this))
            ret.push_back(plane1.intersection(*this));
        Plane plane2(corner1,corner3,corner4) ;
        if(plane2.intersects(*this))
            ret.push_back(plane2.intersection(*this));
        Plane plane3(corner1,corner2,corner4) ;
        if(plane3.intersects(*this))
            ret.push_back(plane3.intersection(*this));
        Plane plane4(corner5,corner6,corner7) ;
        if(plane4.intersects(*this))
            ret.push_back(plane4.intersection(*this));
        Plane plane5(corner6,corner7,corner8) ;
        if(plane5.intersects(*this))
            ret.push_back(plane5.intersection(*this));
        Plane plane6(corner5,corner6,corner8) ;
        if(plane6.intersects(*this))
            ret.push_back(plane6.intersection(*this));
        Plane plane7(corner5,corner7,corner8) ;
        if(plane7.intersects(*this))
            ret.push_back(plane7.intersection(*this));

        return ret ;

    }
    case TETRAHEDRON :
    {
        std::vector<Point> ret ;

        Point corner1 = g->getBoundingPoint(0) ;
        Point corner2 = g->getBoundingPoint(1) ;
        Point corner3 = g->getBoundingPoint(2) ;
        Point corner4 = g->getBoundingPoint(3) ;

        TriPoint t0(&corner1, &corner2, &corner3) ;
        if(intersects(t0))
            ret.push_back(intersection(t0)[0]);
        TriPoint t1(&corner2, &corner3, &corner4) ;
        if(intersects(t1))
            ret.push_back(intersection(t1)[0]);
        TriPoint t2(&corner1, &corner2, &corner4) ;
        if(intersects(t2))
            ret.push_back(intersection(t2)[0]);
        TriPoint t3(&corner1, &corner3, &corner4) ;
        if(intersects(t3))
            ret.push_back(intersection(t3)[0]);

        return ret ;
    }
    case POLYGON_PRISM:
    {
        std::vector<Point> ret ;
        
        if(!intersects(g))
            return ret;
        
        Point dir = first() - second();

        std::vector<bool> ins ;
        std::vector<Point> tested ;
        for(double i = 0 ; i < 1 ; i+=.1 )
        {
            Point toTest = second()+i*dir ;
            ins.push_back(g->in(toTest)) ;
            tested.push_back(toTest);
        }
        
        std::vector< std::pair<Point,Point> > bounds ;
        std::vector< std::pair<bool,bool> > boundsin ;
        Point currentFirst = tested[0] ;
        bool currentIn = ins[0] ;
        for(size_t i = 1 ; i < tested.size() ; i++)
        {
            if(ins[i] != currentIn)
            {
               bounds.push_back(std::make_pair(currentFirst, tested[i]));
               boundsin.push_back(std::make_pair(currentIn,ins[i]));
               currentFirst = tested[i] ;
               currentIn = ins[i] ;
            }
        }
        
        for(size_t i =0 ; i < bounds.size() ; i++)
        {
            for(size_t j = 0 ; j < 16 ; j++)
            {
                Point test = (bounds[i].first + bounds[i].second)*.5 ;
                bool testin = g->in(test) ;
                if(testin == boundsin[i].first)
                {
                    bounds[i].first = test ;
                }
                else
                {
                     bounds[i].second = test ;
                }
            }
            
            ret.push_back((bounds[i].first + bounds[i].second)*.5);
        }
        return ret ;
    }
    case LOFTED_POLYGON:
    {
        std::vector<Point> ret ;
        
        if(!intersects(g))
            return ret;
        
        Point dir = first() - second();

        std::vector<bool> ins ;
        std::vector<Point> tested ;
        for(double i = 0 ; i < 1 ; i+=.1 )
        {
            Point toTest = second()+i*dir ;
            ins.push_back(g->in(toTest)) ;
            tested.push_back(toTest);
        }
        
        std::vector< std::pair<Point,Point> > bounds ;
        std::vector< std::pair<bool,bool> > boundsin ;
        Point currentFirst = tested[0] ;
        bool currentIn = ins[0] ;
        for(size_t i = 1 ; i < tested.size() ; i++)
        {
            if(ins[i] != currentIn)
            {
               bounds.push_back(std::make_pair(currentFirst, tested[i]));
               boundsin.push_back(std::make_pair(currentIn,ins[i]));
               currentFirst = tested[i] ;
               currentIn = ins[i] ;
            }
        }
        
        for(size_t i = 0 ; i < bounds.size() ; i++)
        {
            for(size_t j = 0 ; j < 16 ; j++)
            {
                Point test = (bounds[i].first + bounds[i].second)*.5 ;
                bool testin = g->in(test) ;
                if(testin == boundsin[i].first)
                {
                    bounds[i].first = test ;
                }
                else
                {
                     bounds[i].second = test ;
                }
            }
            
            ret.push_back((bounds[i].first + bounds[i].second)*.5);
        }
        return ret ;
    }

    default:
        return std::vector<Point>(0) ;
    }
}


bool Segment::intersects(const Segment & l) const
{
    if (isAligned(l.first(), s, f) && isAligned(l.second(), s, f))
    {
        return l.on(f) || l.on(s) || on(l.first()) || on(l.second()) ;
    }

    if(l.on(f) || l.on(s) || on(l.first()) || on(l.second()))
        return true ;

    Matrix m(2,2) ;
    Vector v(2) ;

    if(std::abs(vec.getX()) < POINT_TOLERANCE && std::abs(l.vector().getX()) < POINT_TOLERANCE)
    {
        return false ;
    }

    if(std::abs(vec.getY()) < POINT_TOLERANCE && std::abs(l.vector().getY()) < POINT_TOLERANCE)
    {
        return false ;
    }
    if(dist(vec/vec.norm(), l.vector()/l.vector().norm()) < POINT_TOLERANCE)
    {
        return false ;
    }

    m[0][0] = -vec.getX() ;
    m[0][1] = l.vector().getX() ;
    m[1][0] = -vec.getY() ;
    m[1][1] = l.vector().getY() ;

    v[0] = -l.first().getX() + f.getX() ;
    v[1] = -l.first().getY() + f.getY() ;

    invert2x2Matrix(m) ;

    Vector fac = m * v ;

    Point intersect = f + vec*fac[0];

    return on(intersect) && l.on(intersect) ;

}

Point Segment::project(const Point & p) const
{
    Line l(f, vec) ;

    Point candidate = l.projection(p) ;
    if(on(candidate))
    {
        return candidate ;
    }


    if(squareDist3D(candidate, f) > squareDist3D(candidate, s))
        return s ;

    return f ;
}

bool Segment::intersects(const Point & a, const Point & b) const
{
    return intersects(Segment(a,b));
}

bool Segment::on(const Point &p) const
{
    //origin is s

    double ds = dist(s, p) ;
    double df = dist(f, p) ;

    if(ds < 2.*POINT_TOLERANCE || df < 2.*POINT_TOLERANCE)
        return true ;
    if(!isAligned(p, f, s))
        return false ;

    return std::abs(ds + df - dist(f, s)) < 2.*POINT_TOLERANCE ;

}

void Segment::setFirst(const Point & p)
{
    f = p ;
    mid = f*0.5 + s*0.5;
    vec = f-s ;
}

void Segment::setFirst(double x, double y)
{
    f.set(x, y) ;
    mid = f*0.5 + s*0.5;
    vec = f-s ;
}

void Segment::setSecond(const Point & p)
{
    s = p ;
    mid = f*0.5 + s*0.5;
    vec = f-s ;
}

void Segment::setSecond(double x, double y)
{
    s.set(x, y) ;
    mid = f*0.5 + s*0.5;
    vec = f-s ;
}

void Segment::set(const Point & p0, const Point & p1)
{
    f.set(p0) ;
    s.set(p1) ;
    mid = f*0.5 + s*0.5;
    vec = f-s ;
}

void Segment::set(double x0, double y0, double x1, double y1)
{
    f.set(x0, y0) ;
    s.set(x1, y1) ;
    mid = f*0.5 + s*0.5;
    vec = f-s ;
}

const Point & Segment::first() const
{
    return f ;
}

const Point & Segment::second() const
{
    return s ;
}

const Point & Segment::midPoint() const
{
    return mid ;
}

const Point & Segment::vector() const
{
    return vec ;
}


Point Segment::intersection(const Line & l) const
{
    Matrix m(2,2) ;
    Vector v(2) ;

    m[0][0] = vec.getX() ;
    m[0][1] = -l.vector().getX() ;
    m[1][0] = vec.getY() ;
    m[1][1] = -l.vector().getY() ;

    v[0] = l.origin().getX() - f.getX() ;
    v[1] = l.origin().getY() - f.getY() ;

    invert2x2Matrix(m) ;

    Vector fac = m * v ;

    return f + vec*fac[0];
}

TriPoint::TriPoint(const Point * p0, const Point * p1, const Point * p2) : point(3)
{
    point[0] = *p0 ;
    point[1] = *p1 ;
    point[2] = *p2 ;
    normal = (*p0-*p1)^(*p2-*p1) ;
    center = (*p0+*p1+*p2)/3. ;
    double n =  normal.norm() ;
    if(n > POINT_TOLERANCE)
        normal /= n ;

}

TriPoint::TriPoint(const Point & p0, const Point & p1, const Point & p2) : point(3)
{
    point[0] = p0 ;
    point[1] = p1 ;
    point[2] = p2 ;
    normal = (p0-p1)^(p2-p1) ;
    center = (p0+p1+p2)/3. ;
    double n =  normal.norm() ;
    if(n > POINT_TOLERANCE)
        normal /= n ;

}

double TriPoint::area() const
{
    return .5* ((point[0]-point[1])^(point[2]-point[1])).norm() ;
}

Vector TriPoint::normalv() const
{
    Vector ret(3) ;
    ret[0] = normal.getX() ;
    ret[1] = normal.getY() ;
    ret[2] = normal.getZ() ;
    return ret ;
}

Vector TriPoint::normalv(const Point & p) const
{
    bool sameSide = isOnTheSameSide(p, center+normal, first(), second(), third()) ;
    double sign = 1. ;
    if(!sameSide)
        sign = -1 ;
    Vector ret(3) ;
    ret[0] = normal.getX() ;
    ret[1] = normal.getY() ;
    ret[2] = normal.getZ() ;
    return ret/sign ;

}

Point TriPoint::projection(const Point & p) const
{
    Plane plane(point[0], normal) ;

    Point planProj = plane.projection(p) ;

    if(in(planProj))
        return planProj ;

    Segment s0(point[0], point[1]) ;
    Segment s1(point[0], point[2]) ;
    Segment s2(point[1], point[2]) ;

    Point sproj0 = s0.project(planProj) ;
    double d0 = squareDist3D(p, sproj0) ;
    Point sproj1 = s1.project(planProj) ;
    double d1 = squareDist3D(p, sproj1) ;
    Point sproj2 = s2.project(planProj) ;
    double d2 = squareDist3D(p, sproj2) ;

    if(d0 <= d1 && d0 <= d2)
        return sproj0 ;
    if(d1 <= d0 && d1 <= d2)
        return sproj1 ;

    return sproj2 ;

}

bool TriPoint::in(const Point & p) const
{
    Point u = point[1]-point[0] ;
    Point v = point[2]-point[0] ;

    double uu = u*u ;
    double uv = u*v ;
    double vv = v*v ;
    Point w = p - point[0] ;
    double wu = w*u ;
    double wv = w*v ;
    double d = uv*uv-uu*vv ;

    double q = (uv*wv - vv*wu) / d;

    if(q < 0. || q > 1.)
        return false ;

    double t = (uv * wu - uu * wv) / d;

    if (t < -POINT_TOLERANCE || (q + t) > 1.+POINT_TOLERANCE)
        return false;

    return true;
}

std::valarray<std::pair<Point, double> > TriPoint::getGaussPoints(bool timeDependent) const
{
    std::valarray< std::pair<Point, double> > gp(4+4*timeDependent) ;
    Point origin(point[1]) ;
    Point y(point[0]) ;
    Point x(point[2]) ;
    y -= origin ;
    x -= origin ;
    Point a = origin + x*0.2 + y*0.2 ;
    Point b = origin + x*0.6 + y*0.2 ;
    Point c = origin + x*0.2 + y*0.6 ;
    Point d = (point[0]+point[1]+point[2])/3.0 ;
//	double n = norm() ;

    double ar = area()*2. ;

    if(!timeDependent)
    {
        gp[0] = std::pair<Point, double>(a,  0.260416666666667*ar) ;
        gp[1] = std::pair<Point, double>(b,  0.260416666666667*ar) ;
        gp[2] = std::pair<Point, double>(c,  0.260416666666667*ar) ;
        gp[3] = std::pair<Point, double>(d, -0.28125*ar) ;
    }
    else
    {
        gp[0] = std::pair<Point, double>(Point(a.getX(), a.getY(), a.getZ(), -0.577350269189626), 0.260416666666667*ar) ;
        gp[1] = std::pair<Point, double>(Point(b.getX(), b.getY(), b.getZ(), -0.577350269189626), 0.260416666666667*ar) ;
        gp[2] = std::pair<Point, double>(Point(c.getX(), c.getY(), c.getZ(), -0.577350269189626), 0.260416666666667*ar) ;
        gp[3] = std::pair<Point, double>(Point(d.getX(), d.getY(), d.getZ(), -0.577350269189626), -0.28125*ar) ;
        gp[4] = std::pair<Point, double>(Point(a.getX(), a.getY(), a.getZ(), 0.577350269189626), 0.260416666666667*ar) ;
        gp[5] = std::pair<Point, double>(Point(b.getX(), b.getY(), b.getZ(), 0.577350269189626), 0.260416666666667*ar) ;
        gp[6] = std::pair<Point, double>(Point(c.getX(), c.getY(), c.getZ(), 0.577350269189626), 0.260416666666667*ar) ;
        gp[7] = std::pair<Point, double>(Point(d.getX(), d.getY(), d.getZ(), 0.577350269189626), -0.28125*ar) ;
    }


    return gp ;
}


bool Segment::intersects(const TriPoint &g) const
{

    Point v(g.point[1]-g.point[0]) ;
    Point u(g.point[2]-g.point[0]) ;

    Point dir = first() - second();
    Point w0 = second() - g.point[0];
    double a = -((u ^ v)*w0) ;
    double b = (u ^ v)*dir;
    if (std::abs(b) < POINT_TOLERANCE*POINT_TOLERANCE)
    {
        if (std::abs(a) < POINT_TOLERANCE*POINT_TOLERANCE)
            return true;

        return false ;
    }

    double r = a / b ;

    if(r < -POINT_TOLERANCE || r > 1.+POINT_TOLERANCE)
        return false ;

    Point intersect(second()+dir*r) ;
    Point w(intersect-g.point[0]) ;

    double uv = u*v;
    double wu = u*w;
    double wv = v*w;
    double uu = u*u;
    double vv = v*v;
    double d = uv*uv -uu*vv ;
    double s = (uv*wv-vv*wu)/d ;
    if (s < -POINT_TOLERANCE || s > 1.+POINT_TOLERANCE)        // I is outside T
        return false;
    double t = (uv*wu-uu*wv)/d ;
    if (t < -POINT_TOLERANCE || (s + t) > 1.+POINT_TOLERANCE)
        return false ;

    return true ;
}

Point Segment::intersection(const Segment &l) const
{
    if (isAligned(l.first(), s, f) && isAligned(l.second(), s, f))
    {
        if(on(l.first()) && on(l.second())) ;
            return l.midPoint() ;
        if(on(l.first()))
            return l.first() ;
        if(on(l.second()))
            return l.second() ;
        if(l.on(f) && l.on(s)) ;
            return midPoint() ;
        if(l.on(f))
            return f ;
        if(l.on(s))
            return s ;

    }

    Matrix m(2,2) ;
    Vector v(2) ;

    m[0][0] = vec.getX() ;
    m[0][1] = -l.vector().getX() ;
    m[1][0] = vec.getY() ;
    m[1][1] = -l.vector().getY() ;

    v[0] = l.second().getX() - s.getX() ;
    v[1] = l.second().getY() - s.getY() ;

    invert2x2Matrix(m) ;

    Vector fac = m * v ;
    return s + vec*fac[0];
}

std::vector<Point> Segment::intersection(const TriPoint &g) const
{
    std::vector<Point> ret ;
    if(isCoplanar(f, g.point[0],g.point[1],g.point[2]))
        if( g.in(f))
        {
            ret.push_back(f) ;
        }

    if(isCoplanar(s, g.point[0],g.point[1],g.point[2]))
        if( g.in(s))
        {
            ret.push_back(s) ;
        }

    Vector vec(3) ;
    vec[0] = f.getX() - g.point[0].getX() ;
    vec[1] = f.getY() - g.point[0].getY() ;
    vec[2] = f.getZ() - g.point[0].getZ() ;

    Matrix mat(3,3) ;
    mat[0][0] = f.getX()-s.getX();
    mat[0][1] = g.point[1].getX()-g.point[0].getX();
    mat[0][2] = g.point[2].getX()-g.point[0].getX();
    mat[1][0] = f.getY()-s.getY();
    mat[1][1] = g.point[1].getY()-g.point[0].getY();
    mat[1][2] = g.point[2].getY()-g.point[0].getY();
    mat[2][0] = f.getZ()-s.getZ();
    mat[2][1] = g.point[1].getZ()-g.point[0].getZ();
    mat[2][2] = g.point[2].getZ()-g.point[0].getZ();
    Vector tuv = inverse3x3Matrix(mat)*vec ;
    ret.push_back(f+ (s-f)*tuv[0]) ;
    return ret ;
}

bool isInTriangle(const Point & test, const Point&  p0, const Point & p1, const Point  &p2)
{
    return isOnTheSameSide( test, p0, p1, p2) && isOnTheSameSide(test, p1, p0, p2) && isOnTheSameSide(test, p2, p1, p2) ;
}

bool isOnTheSameSide(const Point & test, const Point & witness, const Point & f0, const Point & f1, double norm)
{
    
    Point mid = (f0+f1)*.5 ;
    return (mid-test)*(mid-witness) > 0 ;
    Point frontier(f1.getX()*norm-f0.getX()*norm,f1.getY()*norm-f0.getY()*norm,f1.getZ()*norm-f0.getZ()*norm) ;
    Point yes(witness.getX()*norm-f0.getX()*norm,witness.getY()*norm-f0.getY()*norm,witness.getZ()*norm-f0.getZ()*norm) ;
    Point perhaps(test.getX()*norm-f0.getX()*norm,test.getY()*norm-f0.getY()*norm,test.getZ()*norm-f0.getZ()*norm) ;
    return (frontier^yes).getZ()*(frontier^perhaps).getZ() > -POINT_TOLERANCE ;
}

bool isOnTheSameSide(const Point * test, const Point *witness, const Point *f0, const Point *f1, double norm)
{
    return isOnTheSameSide(*test, *witness, *f0, *f1, norm) ;

}

bool isOnTheSameSide(const Point & test, const Point & witness, const Point & f0, const Point & f1, const Point & f2, double renorm)
{
    Point mid = (f0+f1+f2)*.3333333 ;
    return (mid-test)*(mid-witness) > 0 ;
    
    Point f2test((f2.getX()-test.getX())*renorm,(f2.getY()-test.getY())*renorm,(f2.getZ()-test.getZ())*renorm) ;
    Point f1test((f1.getX()-test.getX())*renorm,(f1.getY()-test.getY())*renorm,(f1.getZ()-test.getZ())*renorm) ;
    Point f0test((f0.getX()-test.getX())*renorm,(f0.getY()-test.getY())*renorm,(f0.getZ()-test.getZ())*renorm) ;
    Point f2witness((f2.getX()-witness.getX())*renorm,(f2.getY()-witness.getY())*renorm,(f2.getZ()-witness.getZ())*renorm) ;
    Point f1witness((f1.getX()-witness.getX())*renorm,(f1.getY()-witness.getY())*renorm,(f1.getZ()-witness.getZ())*renorm) ;
    Point f0witness((f0.getX()-witness.getX())*renorm,(f0.getY()-witness.getY())*renorm,(f0.getZ()-witness.getZ())*renorm) ;
    return ((f2test.getX()*(f1test.getY()*f0test.getZ() - f0test.getY()*f1test.getZ())-f2test.getY()*(f1test.getX()*f0test.getZ() - f0test.getX()*f1test.getZ())+f2test.getZ()*(f1test.getX()*f0test.getY() - f0test.getX()*f1test.getY()))*(f2witness.getX()*(f1witness.getY()*f0witness.getZ() - f0witness.getY()*f1witness.getZ())-f2witness.getY()*(f1witness.getX()*f0witness.getZ() - f0witness.getX()*f1witness.getZ())+f2witness.getZ()*(f1witness.getX()*f0witness.getY() - f0witness.getX()*f1witness.getY())) > 0) ;
}

bool isOnTheSameSide(const Point * test, const Point * witness, const Point * f0, const Point * f1, const Point * f2, double renorm)
{
    return isOnTheSameSide(*test, *witness, *f0, *f1, *f2, renorm) ;
// 	return (((f2->getX()-test->getX())*((f1->getY()-test->getY())*(f0->getZ()-test->getZ()) - (f0->getY()-test->getY())*(f1->getZ()-test->getZ()))-(f2->getY()-test->getY())*((f1->getX()-test->getX())*(f0->getZ()-test->getZ()) - (f0->getX()-test->getX())*(f1->getZ()-test->getZ()))+(f2->getZ()-test->getZ())*((f1->getX()-test->getX())*(f0->getY()-test->getY()) - (f0->getX()-test->getX())*(f1->getY()-test->getY())))*((f2->getX()-witness->getX())*((f1->getY()-witness->getY())*(f0->getZ()-witness->getZ()) - (f0->getY()-witness->getY())*(f1->getZ()-witness->getZ()))-(f2->getY()-witness->getY())*((f1->getX()-witness->getX())*(f0->getZ()-witness->getZ()) - (f0->getX()-witness->getX())*(f1->getZ()-witness->getZ()))+(f2->getZ()-witness->getZ())*((f1->getX()-witness->getX())*(f0->getY()-witness->getY()) - (f0->getX()-witness->getX())*(f1->getY()-witness->getY()))) > 0) ;
}

double dist(const Point & v1, const Point & v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r ;
//	temp = _mm_sub_pd(v1.veczt, v2.veczt) ;
    temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    r.vec = _mm_dp_pd(temp, temp, 61) ;
//	r.vec += _mm_dp_pd(temp, temp, 62) ;
    double z = vi.veczt.getZ() - v2.veczt.getZ() ;
    return sqrt(r.val[0]+ r.val[1] + z*z);
#elif defined HAVE_SSE3
//	vecdouble rzt ;
    vecdouble rxy ;
// 	rzt.vec = _mm_sub_pd(v1.veczt, v2.veczt) ;
// 	rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    rxy.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
    double z = vi.veczt.getZ() - v2.veczt.getZ() ;
    return sqrt(z*z + rxy.val[0]+ rxy.val[1]);
#else
    double x = v1.getX()-v2.getX() ;
    double y = v1.getY()-v2.getY() ;
    double z = v1.getZ()-v2.getZ() ;
//	double t = v1.getT()-v2.getT() ;
    return sqrt(x*x+y*y+z*z) ;
#endif
}

double dist(const Point * v1, const Point * v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r ;
    temp = _mm_sub_pd(v1->veczt, v2->veczt) ;
    r.vec = _mm_dp_pd(temp, temp, 61) ;
    temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    r.vec += _mm_dp_pd(temp, temp, 62) ;
    return sqrt(r.val[0]+ r.val[1] );
#elif defined HAVE_SSE3
    vecdouble rzt ;
    vecdouble rxy ;
    rzt.vec = _mm_sub_pd(v1->veczt, v2->veczt) ;
    rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    rxy.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
    return sqrt(rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1]);
#else
    double x = v1->getX()-v2->getX() ;
    double y = v1->getY()-v2->getY() ;
    double z = v1->getZ()-v2->getZ() ;
    double t = v1->getT()-v2->getT() ;
    return sqrt(x*x+y*y+z*z+t*t) ;
#endif
}



Matrix rotateToVector (Point * toRotate, const Point & toVector)
{
    Point x= *toRotate^toVector ;
    double n = x.sqNorm() ;
    if(n < POINT_TOLERANCE)
    {
        return identity(3) ;
    }
    
    double costheta = *toRotate*toVector/sqrt(toRotate->sqNorm()*toVector.sqNorm()) ;
    double theta = acos(costheta);
    
    if(std::abs(theta-M_PI) < POINT_TOLERANCE)
    {
        x = Point(toRotate->getY(), -toRotate->getX(), 0) ;
        x /= x.norm() ;
    }
    else
        x /= sqrt(n) ;
    
    Matrix A(3,3) ;
    A[0][1] = -x.getZ(); A[0][2] = x.getY();
    A[1][0] = x.getZ() ; A[1][2] = -x.getX() ;
    A[2][0] = -x.getY(); A[2][1] = x.getX();
    
    Matrix rot = identity(3) ;
    rot += A*sin(theta) + (A*A)*(1.-costheta);
    
    *toRotate *= rot ;
    return rot ;
}

double squareDist(const  Point &v1, const Point & v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r ;
    temp = _mm_sub_pd(v1.veczt, v2.veczt) ;
    r.vec = _mm_dp_pd(temp, temp, 61) ;
    temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    r.vec += _mm_dp_pd(temp, temp, 62) ;
    return r.val[0]+ r.val[1] ;
#elif defined HAVE_SSE3
    vecdouble rzt ;
    vecdouble rxy ;
    rzt.vec = _mm_sub_pd(v1.veczt, v2.veczt) ;
    rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    rxy.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
    return rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1] ;
#else
    double x = v1.getX()-v2.getX() ;
    double y = v1.getY()-v2.getY() ;
    double z = v1.getZ()-v2.getZ() ;
    double t = v1.getT()-v2.getT() ;
    return x*x+y*y+z*z+t*t ;
#endif
}

double squareDist(const Point *v1, const Point *v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r ;
    temp = _mm_sub_pd(v1->veczt, v2->veczt) ;
    r.vec = _mm_dp_pd(temp, temp, 61) ;
    temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    r.vec += _mm_dp_pd(temp, temp, 62) ;
    return r.val[0]+ r.val[1];
#elif defined HAVE_SSE3
    vecdouble rzt ;
    vecdouble rxy ;
    rzt.vec = _mm_sub_pd(v1->veczt, v2->veczt) ;
    rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    rxy.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
    return rzt.val[0]+ rzt.val[1] + rxy.val[0]+ rxy.val[1] ;
#else
    double x = v1->getX()-v2->getX() ;
    double y = v1->getY()-v2->getY() ;
    double z = v1->getZ()-v2->getZ() ;
    double t = v1->getT()-v2->getT() ;
    return x*x+y*y+z*z+t*t ;
#endif
}

double squareDist2D(const  Point &v1, const Point & v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r ;
    temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    r.vec = _mm_dp_pd(temp, temp, 61) ;
    return r.val[0] ;
#elif defined HAVE_SSE3
    vecdouble rzt ;
    rzt.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    return rzt.val[0]+ rzt.val[1];
#else
    double x = v1.getX()-v2.getX() ;
    double y = v1.getY()-v2.getY() ;
    return x*x+y*y ;
#endif
}

double squareDist2D(const Point *v1, const Point *v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r ;
    temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    r.vec = _mm_dp_pd(temp, temp, 61) ;
    return r.val[0] ;
#elif defined HAVE_SSE3
    vecdouble rzt ;
    rzt.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    return rzt.val[0]+ rzt.val[1] ;
#else
    double x = v1->getX()-v2->getX() ;
    double y = v1->getY()-v2->getY() ;
    return x*x+y*y ;
#endif
}

double squareDist3D(const  Point &v1, const Point & v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r0 ;
    vecdouble r1 ;
    temp = _mm_sub_pd(v1.veczt, v2.veczt) ;
    r0.vec = _mm_dp_pd(temp, temp, 61) ;
    temp = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    r1.vec = _mm_dp_pd(temp, temp, 62) ;
    return r0.val[0]+ r1.val[0]+r1.val[1] ;
#elif defined HAVE_SSE3
    vecdouble rzt ;
    vecdouble rxy ;
    rzt.vec = _mm_sub_pd(v1.veczt, v2.veczt) ;
    rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    rxy.vec = _mm_sub_pd(v1.vecxy, v2.vecxy) ;
    rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
    return rzt.val[0]/*+ rzt.val[1]*/ + rxy.val[0]+ rxy.val[1] ;
#else
    double x = v1.getX()-v2.getX() ;
    double y = v1.getY()-v2.getY() ;
    double z = v1.getZ()-v2.getZ() ;
    return x*x+y*y+z*z ;
#endif
}



double squareDist3D(const Point *v1, const Point *v2)
{
#ifdef HAVE_SSE4
    __m128d temp ;
    vecdouble r0 ;
    vecdouble r1 ;
    temp = _mm_sub_pd(v1->veczt, v2->veczt) ;
    r0.vec = _mm_dp_pd(temp, temp, 61) ;
    temp = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    r1.vec = _mm_dp_pd(temp, temp, 62) ;
    return r0.val[0]+ r1.val[0]+r1.val[1] ;
#elif defined HAVE_SSE3
    vecdouble rzt ;
    vecdouble rxy ;
    rzt.vec = _mm_sub_pd(v1->veczt, v2->veczt) ;
    rzt.vec = _mm_mul_pd(rzt.vec, rzt.vec) ;
    rxy.vec = _mm_sub_pd(v1->vecxy, v2->vecxy) ;
    rxy.vec = _mm_mul_pd(rxy.vec, rxy.vec) ;
    return rzt.val[0]/*+ rzt.val[1]*/ + rxy.val[0]+ rxy.val[1] ;
#else
    double x = v1->getX()-v2->getX() ;
    double y = v1->getY()-v2->getY() ;
    double z = v1->getZ()-v2->getZ() ;
    return x*x+y*y+z*z ;
#endif
}


ConvexPolygon* convexHull(const std::vector<Point *> * points)
{
    Point *pivot = (*points)[0]  ;
    for(size_t i = 1 ; i < points->size() ; i++)
    {
        if (pivot->getY() < (*points)[i]->getY())
            pivot = (*points)[i] ;
    }
    std::cout << "pivot = " << pivot->getX() << ", " << pivot->getY() << std::endl ;

    //!then we build a map ordered by the angle to the pivot.

    std::map< double, Point* >  pointSet ;

    for(size_t i = 0 ; i < points->size() ; i++)
    {
        if ((*(*points)[i]) != (*pivot))
        {
            double angle = atan2(pivot->getY()-(*points)[i]->getY(), pivot->getX()-(*points)[i]->getX()) ;

            if(pointSet.find(angle) != pointSet.end())
            {
                std::cout << "already there..." << std::endl ;
                if(squareDist(pointSet[angle], pivot) < squareDist((*points)[i], pivot))
                {
                    pointSet[angle] = (*points)[i] ;
                }
            }
            else
                pointSet[angle] = (*points)[i] ;
        }
    }

    pointSet[-10] = pivot ;

    //!we build the convex hull : we are on the convex hull if and only if we only turned left.

    std::vector<Point*> temphull ;
    std::vector<Point*> orderedset ;

    for(auto i = pointSet.begin() ; i!= pointSet.end() ; ++i)
    {
        std::cout << "point " <<  i->second->getX() << ", " <<  i->second->getY() << std::endl ;
        orderedset.push_back(i->second) ;
    }

    temphull.push_back(orderedset[0] ) ;
    temphull.push_back(orderedset[1] ) ;

    for(auto i = orderedset.begin()+2 ; i != orderedset.end() ; ++i)
    {
        for(size_t j = 0 ; j < temphull.size() ; j++)
            std::cout << "(" << temphull[j]->getX() << ", " << temphull[j]->getY() << ")" << std::endl ;

        //! this is a usual cross product of two vectors...
        if(  ((*i)->getY() - (*temphull.rbegin())->getY())*((*temphull.rbegin())->getX() - temphull[temphull.size()-2]->getX() ) -
                ((*i)->getX() - (*temphull.rbegin())->getX())*((*temphull.rbegin())->getY() - temphull[temphull.size()-2]->getY() ) >
                1e-6 )
        {
            temphull.push_back(*i) ;
            std::cout << "new point in hull = " <<  (*i)->getX() << ", " <<  (*i)->getY() << std::endl ;
        }
        else
        {
            while( !(((*i)->getY() - (*temphull.rbegin())->getY())*((*temphull.rbegin())->getX() - temphull[temphull.size()-2]->getX() ) -
                     ((*i)->getX() - (*temphull.rbegin())->getX())*((*temphull.rbegin())->getY() -temphull[temphull.size()-2]->getY() ) >
                     1e-6))
            {
                std::cout << "out of the hull = " <<  (*temphull.rbegin())->getX() << ", " <<  (*temphull.rbegin())->getY() << std::endl ;
                temphull.pop_back();
            }
            temphull.push_back(*i) ;
            std::cout << "new point in hull = " <<  (*i)->getX() << ", " <<  (*i)->getY() << std::endl ;
        }


    }

    ConvexPolygon *hull = new ConvexPolygon(temphull.size()) ;

    std::copy(temphull.begin(), temphull.end(),hull->begin()) ;

    for(size_t i = 0 ; i < hull->size() ; i++)
        std::cout << "(" << hull->getPoint(i)->getX() << ", " << hull->getPoint(i)->getY() << ")" << std::endl ;

    return hull ;
}

OrientableCircle::OrientableCircle(double r,double x, double y, double z, Point n )
{
    gType = ORIENTABLE_CIRCLE ;
    this->center = Point(x,y,z) ;
    if(n.norm() > POINT_TOLERANCE)
    {
        this->normal = n/n.norm() ;
        this->radius = r ;
    }
    else
    {
        this->radius = 0 ;
        normal = Point(1, 0, 0) ;
    }
}

OrientableCircle::OrientableCircle(double r,const Point * p0, Point n )
{
    gType = ORIENTABLE_CIRCLE ;
    this->center = *p0 ;
    if(n.norm() > POINT_TOLERANCE)
    {
        this->normal = n/n.norm() ;
        this->radius = r ;
    }
    else
    {
        this->radius = 0 ;
        normal = Point(1, 0, 0) ;
    }
}

OrientableCircle::OrientableCircle(double r,const Point p0, Point n )
{
    gType = ORIENTABLE_CIRCLE ;
    this->center = p0 ;
    if(n.norm() > POINT_TOLERANCE)
    {
        this->normal = n/n.norm() ;
        this->radius = r ;
    }
    else
    {
        this->radius = 0 ;
        normal = Point(1, 0, 0) ;
    }
}

OrientableCircle::OrientableCircle()
{
    gType = ORIENTABLE_CIRCLE ;
    this->center = Point() ;
    this->normal = Point(0,0,1) ;
    this->radius = 1 ;
}


std::vector<Point> OrientableCircle::getSamplingBoundingPoints(size_t num_points) const
{
    Vector start(3) ;
    start[0] = -normal.getY()*radius ;
    start[1] = normal.getX()*radius ;
    start[2] = 0 ;
    std::vector<Point> ret ;
    if(num_points == 0)
        return ret ;

    if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]))
    {
        start[0] = 0 ;
        start[1] = normal.getZ()*radius ;
        start[2] = -normal.getY()*radius ;
    }

    if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
    {
        start[0] = normal.getZ()*radius ;
        start[1] = 0 ;
        start[2] = -normal.getX()*radius ;
    }


    double t = (1./(double)num_points)*M_PI ;

    double q0 = cos(t) ;
    double st = sin(t) ;
    double q1 = st * normal.getX() ;
    double q2 = st * normal.getY() ;
    double q3 = st * normal.getZ() ;

    Matrix R(3,3) ;

    R[0][0] =  q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R[0][1] = 2.*(q1*q2 - q0*q3) ;
    R[0][2] = 2.*(q1*q3 + q0*q2) ;
    R[1][0] = 2.*(q2*q1 + q0*q3) ;
    R[1][1] = (q0*q0 - q1*q1 + q2*q2 - q3*q3) ;
    R[1][2] = 2.*(q2*q3 - q0*q1) ;
    R[2][0] = 2.*(q3*q1 - q0*q2) ;
    R[2][1] = 2.*(q3*q2 + q0*q1) ;
    R[2][2] = (q0*q0 - q1*q1 - q2*q2 + q3*q3) ;


    for(size_t i = 0 ; i < num_points ; i++)
    {
        start=R*start ;
        ret.push_back(Point( start[0] + center.getX(),
                             start[1] + center.getY(),
                             start[2] + center.getZ())) ;
    }

    return ret ;
}

void OrientableCircle::sampleBoundingSurface(size_t num_points)
{
    if(num_points == 0)
        return ;

    Vector start(3) ;
    start[0] = -normal.getY()*radius ;
    start[1] = normal.getX()*radius ;
    start[2] = 0 ;

    if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]))
    {
        start[0] = 0 ;
        start[1] = normal.getZ()*radius ;
        start[2] = -normal.getY()*radius ;
    }

    if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
    {
        start[0] = normal.getZ()*radius ;
        start[1] = 0 ;
        start[2] = -normal.getX()*radius ;
    }


    double t = (1./(double)num_points)*M_PI ;

    double q0 = cos(t) ;
    double st = sin(t) ;
    double q1 = st * normal.getX() ;
    double q2 = st * normal.getY() ;
    double q3 = st * normal.getZ() ;

    Matrix R(3,3) ;

    R[0][0] =  q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R[0][1] = 2.*(q1*q2 - q0*q3) ;
    R[0][2] = 2.*(q1*q3 + q0*q2) ;
    R[1][0] = 2.*(q2*q1 + q0*q3) ;
    R[1][1] = (q0*q0 - q1*q1 + q2*q2 - q3*q3) ;
    R[1][2] = 2.*(q2*q3 - q0*q1) ;
    R[2][0] = 2.*(q3*q1 - q0*q2) ;
    R[2][1] = 2.*(q3*q2 + q0*q1) ;
    R[2][2] = (q0*q0 - q1*q1 - q2*q2 + q3*q3) ;

    this->boundingPoints.resize(num_points) ;
    for(size_t i = 0 ; i < num_points ; i++)
    {
        start=R*start ;
        boundingPoints[i] = new Point( start[0] + center.getX(),
                                       start[1] + center.getY(),
                                       start[2] + center.getZ()) ;
    }

}

void OrientableCircle::sampleSurface(size_t num_points)
{
    if(num_points == 0)
        return ;

    sampleBoundingSurface(num_points) ;

    std::vector<Point> toAdd ;

    double dr = 4.*M_PI*radius/num_points ;
    double rad = radius-dr ;

    while(rad > POINT_TOLERANCE)
    {

        size_t nPoints = (size_t)round((double)num_points*(rad/radius)) ;
        if(nPoints < 3)
            break ;

        Vector start(3) ;
        start[0] = -normal.getY()*rad ;
        start[1] = normal.getX()*rad ;
        start[2] = 0 ;

        if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
        {
            start[0] = 0 ;
            start[1] = normal.getZ()*rad ;
            start[2] = -normal.getY()*rad ;
        }

        if(std::abs(start[0]) < POINT_TOLERANCE && std::abs(start[1]) < POINT_TOLERANCE && std::abs(start[2]) < POINT_TOLERANCE)
        {
            start[0] = normal.getZ()*rad ;
            start[1] = 0 ;
            start[2] = -normal.getX()*rad ;
        }

        double t = (1./(double)nPoints)*M_PI ;
        double st = sin(t) ;
        double q0 = cos(t) ;
        double q1 = st * normal.getX() ;
        double q2 = st * normal.getY() ;
        double q3 = st * normal.getZ() ;

        Matrix R(3,3) ;

        R[0][0] =  q0*q0 + q1*q1 - q2*q2 - q3*q3;
        R[0][1] = 2.*(q1*q2 - q0*q3) ;
        R[0][2] = 2.*(q1*q3 + q0*q2) ;
        R[1][0] = 2.*(q2*q1 + q0*q3) ;
        R[1][1] = (q0*q0 - q1*q1 + q2*q2 - q3*q3) ;
        R[1][2] = 2.*(q2*q3 - q0*q1) ;
        R[2][0] = 2.*(q3*q1 - q0*q2) ;
        R[2][1] = 2.*(q3*q2 + q0*q1) ;
        R[2][2] = (q0*q0 - q1*q1 - q2*q2 + q3*q3) ;

        for(size_t i = 0 ; i < nPoints ; i++)
        {
            start=R*start ;
            toAdd.push_back( Point( start[0] + center.getX(),
                                    start[1] + center.getY(),
                                    start[2] + center.getZ())) ;
        }
        rad -= dr ;
    }
    toAdd.push_back(center) ;
    this->inPoints.resize(toAdd.size()) ;
    for(size_t i = 0 ; i < toAdd.size() ; i++)
    {
        inPoints[i] = new Point(toAdd[i]) ;
    }
}

bool OrientableCircle::in(const Point & v) const
{
    return false ;
}

double OrientableCircle::area() const
{
    return radius*radius*M_PI ;
}

double OrientableCircle::volume() const
{
    return 0 ;
}

void OrientableCircle::project(Point * p) const
{
    Plane plane(center, normal) ;
    Point pproj = plane.projection(*p) ;
    if(squareDist3D(center, pproj) > POINT_TOLERANCE*POINT_TOLERANCE)
    {
        Point vec = pproj-center ;
        vec/=vec.norm() ;
        p->set( center + vec*radius) ;
        return ;
    }

    Point p_((double)rand()/RAND_MAX,(double)rand()/RAND_MAX,(double)rand()/RAND_MAX ) ;
    project(&p_);
    *p = p_ ;
}

void OrientableCircle::computeCenter()
{
    return ;
}

double OrientableCircle::getRadius() const
{
    return radius ;
}

bool isCoplanar(const Point *test, const Point *f0, const Point *f1,const Point *f2, double renorm)
{
    return isCoplanar(*test, *f0, *f1, *f2, renorm) ;
} 

double signedAlignement(const Point &test, const Point &f0, const Point &f1)
{
    Point a(f1) ;
    a -= test ;
    Point b(f0) ;
    b -= test ;

    Point n(a^b) ;
    double d = n.norm() ;
    if(d  < POINT_TOLERANCE)
        return 0 ;
    n /= d ;

    return (a^b)*n ;
}

bool isAligned(const Point &test, const Point &f0, const Point &f1)
{

    if(test == f1 || test == f0)
        return true ;
//	if(std::abs((f0 - test) * (f1 - test)) < POINT_TOLERANCE)
//		return false ;

//	Point centre ;
    Point centre(test.getX()+f0.getX()+f1.getX(), test.getY()+f0.getY()+f1.getY(), test.getZ()+f0.getZ()+f1.getZ()) ;
    centre *=.3333333333333333333333333 ;
    Point f0_(f0.getX()-centre.getX(), f0.getY()-centre.getY(), f0.getZ()-centre.getZ()) ;
    Point f1_(f1.getX()-centre.getX(), f1.getY()-centre.getY(), f1.getZ()-centre.getZ()) ;
    Point test_(test.getX()-centre.getX(), test.getY()-centre.getY(), test.getZ()-centre.getZ()) ;
    double scale = sqrt(4.*POINT_TOLERANCE/(std::max(std::max(f0_.sqNorm(), f1_.sqNorm()), test_.sqNorm()))) ;
    f0_   *= scale ;
    f1_   *= scale ;
    test_ *= scale ;
//	if (std::abs(signedAlignement(test_, f0_, f1_)) > 2.*POINT_TOLERANCE)
//		return false ;

    double na = sqrt((f0_.getX()-f1_.getX())*(f0_.getX()-f1_.getX())+(f0_.getY()-f1_.getY())*(f0_.getY()-f1_.getY())+(f0_.getZ()-f1_.getZ())*(f0_.getZ()-f1_.getZ())) ;
    double nb = sqrt((f0_.getX()-test_.getX())*(f0_.getX()-test_.getX())+(f0_.getY()-test_.getY())*(f0_.getY()-test_.getY())+(f0_.getZ()-test_.getZ())*(f0_.getZ()-test_.getZ())) ;
    double nc = sqrt((f1_.getX()-test_.getX())*(f1_.getX()-test_.getX())+(f1_.getY()-test_.getY())*(f1_.getY()-test_.getY())+(f1_.getZ()-test_.getZ())*(f1_.getZ()-test_.getZ())) ;

    if(na >= nb && na >= nc)
    {
        Line l(f0_,Point((f1_.getX()-f0_.getX())/na,(f1_.getY()-f0_.getY())/na,(f1_.getZ()-f0_.getZ())/na)) ;
        Sphere s(POINT_TOLERANCE, test_) ;
        return l.intersects(&s) ;
    }
    if(nb >= na && nb >= nc)
    {
        Line l(f0_,Point((test_.getX()-f0_.getX())/nb,(test_.getY()-f0_.getY())/nb,(test_.getZ()-f0_.getZ())/nb)) ;
        Sphere s(POINT_TOLERANCE, f1_) ;
        return l.intersects(&s) ;
    }
    Line l(f1_,Point((test_.getX()-f1_.getX())/nc,(test_.getY()-f1_.getY())/nc,(test_.getZ()-f1_.getZ())/nc)) ;
    Sphere s(POINT_TOLERANCE, f0_) ;
    return l.intersects(&s) ;
}

bool isAligned(const Point *test, const Point *f0, const Point *f1)
{
    return isAligned(*test, *f0, *f1) ;
} 


int coplanarCount( Point *const* pts, int numpoints, const Point &f0, const Point &f1, const Point &f2, double renorm)
{
    int count = 0 ;
    Point centre = (**(pts)+f0+f1+f2)*.25 ;
    Point normal = (Point(f1.getX()-f0.getX(),f1.getY()-f0.getY(),f1.getZ()-f0.getZ())^Point(f2.getX()-f0.getX(),f2.getY()-f0.getY(),f2.getZ()-f0.getZ()))*renorm ;

    Point f0_(f0.getX()*renorm-centre.getX()*renorm,f0.getY()*renorm-centre.getY()*renorm,f0.getZ()*renorm-centre.getZ()*renorm) ;
    Point f1_(f1.getX()*renorm-centre.getX()*renorm,f1.getY()*renorm-centre.getY()*renorm,f1.getZ()*renorm-centre.getZ()*renorm) ;
    Point f2_(f2.getX()*renorm-centre.getX()*renorm,f2.getY()*renorm-centre.getY()*renorm,f2.getZ()*renorm-centre.getZ()*renorm) ;

    Point AB = Point(f0_.getX()-f1_.getX(),f0_.getY()-f1_.getY(),f0_.getZ()-f1_.getZ())^Point(f2_.getX()-f1_.getX(),f2_.getY()-f1_.getY(),f2_.getZ()-f1_.getZ()) ;
    for(int i = 0 ; i < numpoints ; i++)
    {
        Point test_(**(pts+i)*renorm-centre*renorm) ;

        if(test_ == f1_ || test_ == f0_ || test_ == f2_)
        {
            count++ ;
            continue ;
        }

        double c0 = AB*(f2_-test_) ;
        if(c0*c0 > POINT_TOLERANCE)
            continue ;

        double c1 = AB.getX()*(f2_.getX()-test_.getX()-normal.getX())+AB.getY()*(f2_.getY()-test_.getY()-normal.getY())+AB.getZ()*(f2_.getZ()-test_.getZ()-normal.getZ()) ;
        double c2 = AB.getX()*(f2_.getX()-test_.getX()+normal.getX())+AB.getY()*(f2_.getY()-test_.getY()+normal.getY())+AB.getZ()*(f2_.getZ()-test_.getZ()+normal.getZ()) ;

        bool positive = c0 > 0 || c1 > 0 || c2 > 0 ;
        bool negative = c0 < 0 || c1 < 0 || c2 < 0 ;
        if( positive && negative )
            count++ ;
    }

    return count ;
}

bool isCoplanar(const Point &test, const Point &f0, const Point &f1, const Point &f2, double renorm)
{

    Point centre = (test+f0+f1+f2)*.25 ;
    Point f0_(f0.getX()*renorm-centre.getX()*renorm,f0.getY()*renorm-centre.getY()*renorm,f0.getZ()*renorm-centre.getZ()*renorm) ;
    Point f1_(f1.getX()*renorm-centre.getX()*renorm,f1.getY()*renorm-centre.getY()*renorm,f1.getZ()*renorm-centre.getZ()*renorm) ;
    Point f2_(f2.getX()*renorm-centre.getX()*renorm,f2.getY()*renorm-centre.getY()*renorm,f2.getZ()*renorm-centre.getZ()*renorm) ;
    Point test_(test.getX()*renorm-centre.getX()*renorm,test.getY()*renorm-centre.getY()*renorm,test.getZ()*renorm-centre.getZ()*renorm) ;

    if(test_ == f1_)
        return true ;
    if(test_ == f0_)
        return true ;
    if(test_ == f2_)
        return true ;

    double c0 = signedCoplanarity(test_, f0_, f1_, f2_) ;
    double c02 = c0*c0 ;
    if(c02 > POINT_TOLERANCE)
        return false ;

    Point normal = Point(f1_.getX()-f0_.getX(),f1_.getY()-f0_.getY(),f1_.getZ()-f0_.getZ())^Point(f2_.getX()-f0_.getX(),f2_.getY()-f0_.getY(),f2_.getZ()-f0_.getZ()) ;
//	normal /= scale ;

    Point a(test_.getX()+normal.getX(),test_.getY()+normal.getY(),test_.getZ()+normal.getZ()) ;
    Point b(test_.getX()-normal.getX(),test_.getY()-normal.getY(),test_.getZ()-normal.getZ()) ;

    double c1 = signedCoplanarity(a, f0_, f1_, f2_) ;
    double c2 = signedCoplanarity(b, f0_, f1_, f2_) ;

    bool positive = c0 > 0 || c1 > 0 || c2 > 0 ;
    bool negative = c0 < 0 || c1 < 0 || c2 < 0 ;
    return  positive && negative ;
} 

double coplanarity(const Point *test, const Point *f0, const Point *f1,const Point *f2)
{
    return coplanarity(*test, *f0, *f1, *f2) ;
} 

double coplanarity(const Point &test, const Point &f0, const Point &f1, const Point &f2)
{

    Point A(f0-f1) ;
    Point B(f2-f1) ;
    Point C(f2-test) ;

    return  std::abs(triProduct(A, B, C))  ;
} 

double signedCoplanarity(const Point &test, const Point &f0, const Point &f1, const Point &f2)
{

    Point A(f0-f1) ;
    Point B(f2-f1) ;
    Point C(f2-test) ;

    return  triProduct(A, B, C)  ;
} 

double signedCoplanarity(const Point *test, const Point *f0, const Point *f1,const Point *f2)
{
    return signedCoplanarity(*test, *f0, *f1, *f2) ;
} 

double triProduct(const Point &A, const Point &B, const Point &C)
{
    return (A^B)*C ;
}

}
