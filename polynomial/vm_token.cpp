//
// C++ Implementation: vm_token
//
// Description:
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "vm_token.h"
#include "../elements/elements.h"

using namespace Amie ;


void InHomogeneousProjectionOperation::eval ( double * a, double * b, double * c ) const
{
    double res = 0;
    Point test ( *a, *b ) ;
    Point gProj ( test ) ;
    inGeo->project ( &gProj );
    if ( inGeo->in ( test ) )
    {
        double totaldist ;
        std::vector<double> weight ;
        std::vector<Point> projs ;
        for ( size_t i = 0 ; i < inProjector.size() ; i++ )
        {
            Line l ( inProjector[i].first(), inProjector[i].vector() ) ;

            Point proj = l.projection ( test ) ;

            projs.push_back ( proj );
            double d = dist ( proj, test ) ;
            totaldist += d ;
            weight.push_back ( d );
        }
        projs.push_back ( gProj );
        double d = dist ( gProj, test ) ;
        totaldist += d ;
        weight.push_back ( d );

        double maxn = 0 ;
        double renorm = 0 ;
        for ( size_t i = 0 ; i < inProjector.size() ; i++ )
        {
            double n = inProjector[i].norm() ;
            if ( weight[i] < POINT_TOLERANCE )
            {
                *c = dist ( projs[i], inProjector[i].second() ) /n ;
                return ;
            }
            maxn = std::max ( maxn, n ) ;
            res += ( 1./weight[i] ) *dist ( projs[i], inProjector[i].second() ) /n ;
            renorm += 1./weight[i] ;
        }
        if ( weight.back() < POINT_TOLERANCE )
        {
            *c = 1 ;
            return ;
        }
        res += ( 1./weight.back() ) * ( 1.-dist ( projs.back(), test ) /inGeo->getRadius() ) ;
        renorm += 1./weight.back() ;

        if ( renorm > POINT_TOLERANCE )
        {
            res /= renorm ;
        }
        else
        {
            res = 0 ;
        }
    }
    else
    {
        double totaldist ;
        std::vector<double> weight ;
        std::vector<Point> projs ;
        for ( size_t i = 0 ; i < outProjector.size() ; i++ )
        {
            Line l ( outProjector[i].first(), outProjector[i].vector() ) ;
            Point proj = l.projection ( test ) ;
            projs.push_back ( proj );
            double d = dist ( proj, test ) ;
            totaldist += d ;
            weight.push_back ( d );
        }
        projs.push_back ( gProj );
        double d = dist ( gProj, test ) ;
        totaldist += d ;
        weight.push_back ( d );

        double maxn = 0 ;
        double renorm = 0 ;
        for ( size_t i = 0 ; i < outProjector.size() ; i++ )
        {
            double n = outProjector[i].norm()  ;
            if ( weight[i] < POINT_TOLERANCE )
            {
                *c = dist ( projs[i], outProjector[i].second() ) /n ;
                return ;
            }
            maxn = std::max ( maxn, n ) ;
            res += ( 1./weight[i] ) *dist ( projs[i], outProjector[i].second() ) /n ;
            renorm += 1./weight[i] ;
        }
        if ( weight.back() < POINT_TOLERANCE )
        {
            *c = 1 ;
            return ;
        }

        res += ( 1./weight.back() ) * ( 1.-dist ( projs.back(), test ) /inGeo->getRadius() ) ;
        renorm += 1./weight.back() ;

        if ( renorm > POINT_TOLERANCE )
        {
            res /= renorm ;
        }
        else
        {
            res = 0 ;
        }
    }

    *c =  res;
}


double sign ( const double t )
{
    if ( t < 0 )
    {
        return -1 ;
    }
    if ( t > 0 )
    {
        return 1 ;
    }
    return 0 ;
}

double positivity ( const double t )
{
    return t >= 0 ;
}

double negativity ( const double t )
{
    return t <= 0 ;
}

double interpolate ( const double a, const double b )
{
    return a/ ( b + a ) ;
}



Memory::Memory()
{
} 

Memory::~Memory()
{
} 


Context::Context ( const double & x_, const double & y_ , const double & z_ , const double & t_ , const double & u_ , const double & v_ , const double & w_ ) : x ( x_ ), y ( y_ ), z ( z_ ), t ( t_ ), u ( u_ ), v ( v_ ), w ( w_ )
{

} 

Context::Context ( const double & x_, const double & y_, const double & z_, const double & t_, const double & u_, const double & v_ ) : x ( x_ ), y ( y_ ), z ( z_ ), t ( t_ ), u ( u_ ), v ( v_ ), w ( 0 )
{

} 

Context::Context ( const double & x_, const double & y_, const double & z_, const double & t_, const double & u_ ) : x ( x_ ), y ( y_ ), z ( z_ ), t ( t_ ), u ( u_ ), v ( 0 ), w ( 0 )
{

} 

Context::Context ( const double & x_, const double & y_, const double & z_, const double & t_ ) : x ( x_ ), y ( y_ ), z ( z_ ), t ( t_ ), u ( 0 ), v ( 0 ), w ( 0 )
{

} 


Context::Context ( const double & x_, const double & y_, const double & z_ ) : x ( x_ ), y ( y_ ), z ( z_ ), t ( 0 ), u ( 0 ), v ( 0 ), w ( 0 )
{

} 

Context::Context ( const double & x_, const double & y_ ) : x ( x_ ), y ( y_ ), z ( 0 ), t ( 0 ), u ( 0 ), v ( 0 ), w ( 0 )
{

} 


Context::Context ( const double & x_ ) : x ( x_ ), y ( 0 ), z ( 0 ), t ( 0 ), u ( 0 ), v ( 0 ), w ( 0 )
{

} 

Context::Context() : x ( 0 ), y ( 0 ), z ( 0 ), t ( 0 ), u ( 0 ), v ( 0 ), w ( 0 )
{

} 

void Context::set ( const double & x_, const double & y_ , const double & z_ , const double & t_ , const double & u_ , const double & v_ , const double & w_ )
{
    x = x_;
    y = y_;
    z = z_;
    t = t_;
    u = u_;
    v = v_;
    w = w_;
}

void Context::set ( const double & x_, const double & y_, const double & z_ )
{
    x = x_;
    y = y_;
    z = z_;
}

void Context::set ( const double & x_, const double & y_ )
{
    x = x_;
    y = y_;
}

void Context::set ( const double & x_ )
{
    x = x_;
}

GeometryOperation::GeometryOperation() { } 


PositionOperation::PositionOperation ( const Segment & s_ )
{
    Point vector = s_.normal ( Point ( -1000, -1000 ) );
    w= ( s_.midPoint() +vector*10. ) ;
    s.push_back ( s_ )  ;
}

PositionOperation::PositionOperation ( const std::vector<Segment> & s_ )
{
    for ( size_t i = 0 ; i < s_.size() ; i++ )
    {
        s.push_back ( s_[i] )  ;
    }
    Point vector= s_[0].normal ( Point ( -1000, -1000 ) );
    w= ( s_[0].midPoint() +vector*10. ) ;
}
void PositionOperation::eval ( double * a, double * b, double * c ) const
{
    Point test ( *a,*b ) ;
    int intersections = 0 ;
    for ( size_t i = 0 ; i < s.size() ; i++ )
    {
        if ( s[i].intersects ( test, w ) )
        {
            intersections++ ;
        }
    }

    *c = ( intersections & 1 ) * 2 - 1;
// 		if(intersections%2 == 1)
// 		    *context.memory.getT()op_pos = -1 ;
// 		else
// 		    *context.memory.getT()op_pos = 1 ;
}
GeometryOperation * PositionOperation::getCopy() const
{
    return new PositionOperation ( s ) ;
}

int PositionOperation::adressOffset() const
{
    return -1 ;
}

LineDistanceOperation::LineDistanceOperation ( const Line & l_ ) : l ( l_ )
{
}

void LineDistanceOperation::eval ( double * a, double * b, double * c ) const
{

    Point test ( *a, *b ) ;
    *c=  sqrt ( squareDist2D ( test,l.projection ( test ) ) );
}

GeometryOperation * LineDistanceOperation::getCopy() const
{
    return new LineDistanceOperation ( l ) ;
}

int LineDistanceOperation::adressOffset() const
{
    return -1 ;
}


InHomogeneousProjectionOperation::InHomogeneousProjectionOperation ( Geometry * inGeo, const std::vector<Segment> & inProjector, const std::vector<Segment> &outProjector ) : inGeo ( inGeo ), inProjector ( inProjector ), outProjector ( outProjector )
{
}

void InHomogeneousProjectionOperation::eval ( double * a, double * b, double * c ) const ;

GeometryOperation * InHomogeneousProjectionOperation::getCopy() const
{
    return new InHomogeneousProjectionOperation ( inGeo, inProjector, outProjector ) ;
}

int InHomogeneousProjectionOperation::adressOffset() const
{
    return -1 ;
}


DomainOperation::DomainOperation ( const Geometry * g )
{
    geo = g ;
}

void DomainOperation::eval ( double * a, double * b, double * c ) const
{
    Point p ( *a,*b,*c ) ;
    if ( geo->in ( p ) )
    {
        *c = 1 ;
    }
    else
    {
        *c = -1 ;
    }
}

GeometryOperation * DomainOperation::getCopy() const
{
    return new DomainOperation ( geo ) ;
}

int DomainOperation::adressOffset() const
{
    return -2 ;
}



DomainBinaryOperation::DomainBinaryOperation ( const Geometry * g ) : geo ( g )
{
}

void DomainBinaryOperation::eval ( double * a, double * b, double * c ) const
{

    Point p ( *a, *b ) ;
    if ( geo->in ( p ) )
    {
        *c = 1 ;
    }
    else
    {
        *c = -1 ;
    }
}

GeometryOperation * DomainBinaryOperation::getCopy() const
{
    return new DomainBinaryOperation ( geo ) ;
}

int DomainBinaryOperation::adressOffset() const
{
    return -1 ;
}

PointDistanceBinaryOperation::PointDistanceBinaryOperation ( const Point & p ) :  x0 ( p.getX() ), y0 ( p.getY() )
{
}

void PointDistanceBinaryOperation::eval ( double * a, double * b, double * c ) const
{
    double x = *a-x0 ;
    double y = *b-y0 ;
    *c = sqrt ( x*x+y*y ) ;

}

GeometryOperation * PointDistanceBinaryOperation::getCopy() const
{
    return new PointDistanceBinaryOperation ( Point ( x0, y0 ) ) ;
}

int PointDistanceBinaryOperation::adressOffset() const
{
    return -1 ;
}


PointDistanceTrinaryOperation::PointDistanceTrinaryOperation ( const Point & p ) : x0 ( p.getX() ), y0 ( p.getY() ), z0 ( p.getZ() )
{
}

void PointDistanceTrinaryOperation::eval ( double * a, double * b, double * c ) const
{
    double x = *a-x0 ;
    double y = *b-y0 ;
    double z = *c-z0 ;
    *c = sqrt ( x*x+y*y+z*z ) ;

}

GeometryOperation * PointDistanceTrinaryOperation::getCopy() const
{
    return new PointDistanceTrinaryOperation ( Point ( x0, y0, z0 ) ) ;
}

int PointDistanceTrinaryOperation::adressOffset() const
{
    return -2 ;
}

RotationBinaryOperation::RotationBinaryOperation ( double a ) : cangle ( cos ( a ) ), sangle ( sin ( a ) )
{
}

void RotationBinaryOperation::eval ( double * a, double * b, double * c ) const
{

    double x = *a ;
    double y =  *b ;
    *a = x*cangle + y*sangle ;
    *b = -x*sangle + y*cangle ;

}
GeometryOperation * RotationBinaryOperation::getCopy() const
{
    return new RotationBinaryOperation ( acos ( cangle ) ) ;
}

int RotationBinaryOperation::adressOffset() const
{
    return 0 ;
}

AngleBinaryOperation::AngleBinaryOperation ( double a, const Point & p ) :cangle ( cos ( a ) ), sangle ( sin ( a ) ), pivot ( p.getX() *cos ( a ) +p.getY() *sin ( a ), -p.getX() *sin ( a ) +p.getY() *cos ( a ) )
{
}

void AngleBinaryOperation::eval ( double * a, double * b, double * c ) const
{

    double x = *a ;
    double y =  *b ;
    double x_t = x*cangle + y*sangle ;
    double y_t = -x*sangle + y*cangle ;
    *c = atan2 ( y_t-pivot.getY(), x_t-pivot.getX() ) ;

}

GeometryOperation * AngleBinaryOperation::getCopy() const
{
    return new AngleBinaryOperation ( acos ( cangle ), pivot ) ;
}

int AngleBinaryOperation::adressOffset() const
{
    return -1 ;
}

PointSquareDistanceBinaryOperation::PointSquareDistanceBinaryOperation ( const Point & p )
{
    base = p ;
}

void PointSquareDistanceBinaryOperation::eval ( double * a, double * b, double * c ) const
{

    Point p ( *a, *b ) ;

    *c = squareDist2D ( p, base ) ;

}

GeometryOperation * PointSquareDistanceBinaryOperation::getCopy() const
{
    return new PointSquareDistanceBinaryOperation ( base ) ;
}

int PointSquareDistanceBinaryOperation::adressOffset() const
{
    return -1 ;
}



LineOfSightOperation::LineOfSightOperation ( const Point & p,  const Geometry * o ) :  base ( p ), obstruction ( o )
{
}

void LineOfSightOperation::eval ( double * a, double * b, double * c ) const
{

    Point p ( *a, *b ) ;


    *c = Segment ( base, p ).intersects ( obstruction ) ;

}

GeometryOperation * LineOfSightOperation::getCopy() const
{
    return new LineOfSightOperation ( base, obstruction ) ;
}

int LineOfSightOperation::adressOffset() const
{
    return -1 ;
}


ProjectionOperation3D::ProjectionOperation3D ( Segment s_ ) : s ( s_ )
{
}

void ProjectionOperation3D::eval ( double * a, double * b, double * c ) const
{

    *c = dist ( Point ( *a, *b, *c ), s.project ( Point ( *a, *b, *c ) ) ) ;
}

GeometryOperation * ProjectionOperation3D::getCopy() const
{
    return new ProjectionOperation3D ( s ) ;
}

int ProjectionOperation3D::adressOffset() const
{
    return -2 ;
}


ProjectionOperation2D::ProjectionOperation2D ( Segment s_ ) : s ( s_ )
{
}

void ProjectionOperation2D::eval ( double * a, double * b, double * c ) const
{

    *c = dist ( Point ( *a, *b ), s.project ( Point ( *a, *b ) ) ) ;
}

GeometryOperation * ProjectionOperation2D::getCopy() const
{
    return new ProjectionOperation2D ( s ) ;
}

int ProjectionOperation2D::adressOffset() const
{
    return -1 ;
}


ProjectionBinaryOperation::ProjectionBinaryOperation ( const Geometry * s_ ) :  g ( s_ )
{
}

void ProjectionBinaryOperation::eval ( double * a, double * b, double * c ) const
{
    Point p ( *a, *b ) ;
    Point p_ ( p ) ;
    g->project ( &p_ ) ;

    *c = sqrt ( squareDist2D ( p, p_ ) ) ;

}

GeometryOperation * ProjectionBinaryOperation::getCopy() const
{
    return new ProjectionBinaryOperation ( g ) ;
}

int ProjectionBinaryOperation::adressOffset() const
{
    return -1 ;
}



HatEnrichmentAlt::HatEnrichmentAlt(const Geometry * g , const Point & head, const Point & p0, const Point & p1) : g(g), head(head),  p0(p0), p1(p1) { }

void HatEnrichmentAlt::eval(double * a, double * b, double * c) const
{
    Point position ( *a, *b ) ;
    
    if(g->in(position) == g->in(p0))
    {
        Line l(p0, position-p0) ;
        std::vector<Point> intersgeo = l.intersection(g) ;
        if(intersgeo.empty())
        {
            *c = 0 ;
            return ;
        }
        
        Triangle t(p0, p1, head) ;
        Point pmin = intersgeo[0] ;
        if(!t.inCircumCircle(pmin))
            pmin = intersgeo[1] ;
        
        double distTot = dist(p0, pmin);
        double distPos = dist(position, p0) ;
        *c = distPos/distTot ;
        
        return ;
    }
    
    Line l(p1, position-p1) ;
    std::vector<Point> intersgeo = l.intersection(g) ;
    if(intersgeo.empty())
    {
        *c = 0 ;
        return ;
    }
    
    Triangle t(p0, p1, head) ;
    Point pmin = intersgeo[0] ;
    if(!t.inCircumCircle(pmin))
        pmin = intersgeo[1] ;
    
    double distTot = dist(p1, pmin);
    double distPos = dist(position, p1) ;
    *c = distPos/distTot ;
    
    return ;
}

GeometryOperation * HatEnrichmentAlt::getCopy() const
{
    return new HatEnrichmentAlt(g, head, p0, p1) ;
}

int HatEnrichmentAlt::adressOffset() const
{
    return -2 ;
}


HatEnrichment::HatEnrichment(const Geometry * g , const Point & p, const Segment & s) :g(g), p(p), s(s) 
{
}

void HatEnrichment::eval(double * a, double * b, double * c) const
{
    Point position ( *a, *b ) ;
    
    if(g->in(position))
    {
        *c = 0 ;
        return ;
    }
    
    if(p == position)
    {
       *c = 0 ;
       return ;
    }
    if(s.on(position))
    {
       *c = 0 ;
       return ;
    }
    
    Line l(p, position-p) ;
    Line ls(s.first(), s.second()-s.first()) ;
    Point interseg = l.intersection(ls) ;
    
    std::vector<Point> intersgeo = l.intersection(g) ;
    if(intersgeo.size() < 2)
    {
        *c = 0 ;
        return ;
    }

    Triangle t(p, s.first(), s.second()) ;
    Point pmin = intersgeo[0] ;
    
    if(squareDist2D(t.getCircumCenter(), pmin) > squareDist2D(t.getCircumCenter(), intersgeo[1]))
        pmin = intersgeo[1] ;
        
    
   
    double basis = dist(interseg, pmin) ;
    if(g->in(p) == g->in(position))
    { 
        double distTotPoint = dist(p, pmin);
        if(distTotPoint < default_derivation_delta)
        {
            *c = basis*basis ;
            return ;
        }
        double distPos = dist(position, p) ;
        *c = basis*basis*distPos/distTotPoint ;
        return ;
    }
    
    double distTotSeg = dist(interseg, pmin);   
    if(distTotSeg < default_derivation_delta)
    {
        *c = basis*basis ;
        return ;
    }
    
    double distPos = dist(position, interseg) ;
    
    *c =  basis*basis*distPos/distTotSeg ;

}

GeometryOperation * HatEnrichment::getCopy() const 
{
    return new HatEnrichment(g, p, s) ;
}

int HatEnrichment::adressOffset() const 
{ 
    return -2 ;
}


HatEnrichmentDerivative::HatEnrichmentDerivative(const Geometry * g , const Point & p, const Segment & s, Variable v): g(g), p(p), s(s), v(v) { }

void HatEnrichmentDerivative::eval(double * a, double * b, double * c) const
{
    
    Point position ( *a, *b ) ;
    if(g->in(position))
    {
        *c = 0 ;
        return ;
    }
    
    Point positionm(position) ;
    Point positionp(position) ;
    
    double d = 0.0001*std::min(g->getRadius(), 0.5*dist(s.midPoint(), p)) ;
    if(v == XI)
    {
        positionm.getX() -= d ;
        positionp.getX() += d ;
    }
    else if(v == ETA)
    {
        positionm.getY() -= d ;
        positionp.getY() += d ;
    }

    Line lm(p, positionm-p) ;
    Line lp(p, positionp-p) ;    
    Line ls(s.first(), s.second()-s.first()) ;
    
    std::vector<Point> intersgeom = lm.intersection(g) ;
    std::vector<Point> intersgeop = lp.intersection(g) ;
    if(intersgeom.size() < 2 || intersgeop.size() < 2 )
    {
        *c = 0 ;
        return ;
    }
    
    Point intersegm = lm.intersection(ls) ;
    Point intersegp = lp.intersection(ls) ;

    Point pminm = intersgeom[0] ;
    Point pminp = intersgeop[0] ;
    
    Triangle t(p, s.first(), s.second()) ;
    
    if(squareDist2D(t.getCircumCenter(), pminm) > squareDist2D(t.getCircumCenter(), intersgeom[1]))
        pminm = intersgeom[1] ;
    if(squareDist2D(t.getCircumCenter(), pminp) > squareDist2D(t.getCircumCenter(), intersgeop[1]))
        pminp = intersgeop[1] ;
    
    double basis = dist((intersegm+intersegp)*.5, p) ;

    if(g->in(p) == g->in(position))
    {
            
        double distTotmPoint = dist(p, pminm);
        if(distTotmPoint < default_derivation_delta)
        {
            *c = 0 ;
            return ;
        }
        
        double distTotpPoint = dist(p, pminp);
        double distPosm = dist(positionm, p) ;
        double distPosp = dist(positionp, p) ;  
        *c = basis*basis*0.5*(distPosp/distTotpPoint-distPosm/distTotmPoint)/d ;

        return ;
    }
    
    double distTotmSeg = dist(intersegm, pminm);
    if(distTotmSeg < default_derivation_delta)
    {
        *c = 0 ;
        return ;
    }
    double distTotpSeg = dist(intersegp, pminp);                   
    double distPosm = dist(positionm, intersegm) ;
    double distPosp = dist(positionp, intersegp) ;
    *c = basis*basis*0.5*(distPosp/distTotpSeg-distPosm/distTotmSeg)/d ;
    
}

GeometryOperation * HatEnrichmentDerivative::getCopy() const 
{
    return new HatEnrichmentDerivative(g, p, s, v) ;
}

int HatEnrichmentDerivative::adressOffset() const
{
    return -2 ;
}



HatEnrichment3D::HatEnrichment3D(const Geometry * g , const Point & p, const TriPoint & s) :g(g), p(p), s(s) {}

void HatEnrichment3D::eval(double * a, double * b, double * c) const
{
//     std::cout << *a << "  " << *b << "  " << *c << std::endl ;
    Point position ( *a, *b, *c ) ;
    Line l(p, p-position) ;
    Tetrahedron t (p, s.first(), s.second(), s.third()) ;
    
    std::vector<Point> interseg = l.intersection(s) ;
    if(interseg.empty())
    {
//        s.first().print();
//        s.second().print();
//        s.third().print();
//        p.print();
//        position.print();
       *a = 0 ;
       return ;
    }
    std::vector<Point> intersgeo = l.intersection(g) ;
    for(size_t i = 0 ;  i < intersgeo.size() ; i++)
    {
        if(t.in(intersgeo[i]))
        {
            if(g->in(p) == g->in(position))
            {
                double distTot =dist(p, intersgeo[i]) ;
                double distPos = dist(position, p) ;
                *a = distPos/distTot ;
                return ;
            }
            else
            {                
                double distTot = dist(interseg.front(), intersgeo[i]) ;
                double distPos = dist(position, interseg.front()) ;
                *a = distPos/distTot ;
                return ;
            }
        }
    }
    
    *a = 0 ;
}

GeometryOperation * HatEnrichment3D::getCopy() const 
{
    return new HatEnrichment3D(g, p, s) ;
}

int HatEnrichment3D::adressOffset() const 
{ 
    return 0 ;
}


