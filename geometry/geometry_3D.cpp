// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution

#include "geometry_3D.h"

using namespace Mu ;

Tetrahedron::Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3): ConvexGeometry(4)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(4) ;
	boundingPoints[0] = p0 ;
	boundingPoints[1] = p1 ;
	boundingPoints[2] = p2 ;
	boundingPoints[3] = p3 ;
	
	if(this->volume() < 0 )
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
	
	assert(this->volume() >0 ) ;
	computeCenter() ;	
}

Tetrahedron::Tetrahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7): ConvexGeometry(8)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(8) ;
	boundingPoints[0] = p0 ;
	boundingPoints[1] = p1 ;
	boundingPoints[2] = p2 ;
	boundingPoints[3] = p3 ;
	boundingPoints[4] = p4 ;
	boundingPoints[5] = p5 ;
	boundingPoints[6] = p6 ;
	boundingPoints[7] = p7 ;
	
	if(this->volume() < 0 )
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
	assert(this->volume() >0 ) ;
	computeCenter() ;	
}


Tetrahedron::Tetrahedron(Point p0, Point p1, Point p2, Point p3): ConvexGeometry(4)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(4) ;
	boundingPoints[0] = &p0 ;
	boundingPoints[1] = &p1 ;
	boundingPoints[2] = &p2 ;
	boundingPoints[3] = &p3 ;
	if(this->volume() < 0 )
	{
		boundingPoints[0] = &p1; 
		boundingPoints[1] = &p0;
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
}

Tetrahedron::Tetrahedron(): ConvexGeometry(4)
{
	gType = TETRAHEDRON ;
	
	this->boundingPoints.resize(4) ;
	boundingPoints[0] = new Point(1,0,0) ;
	boundingPoints[1] = new Point(0,1,0) ;
	boundingPoints[2] = new Point(0,0,1) ;
	boundingPoints[3] = new Point(0,0,0) ;
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
}

void Tetrahedron::sampleBoundingSurface(size_t num_points)
{
	//! \todo make it do something
	return ;
}

void Tetrahedron::computeCenter()
{
	this->center = Point(0,0,0) ;
	for(size_t i = 0  ; i < boundingPoints.size() ; i++)
	{
		this->center += (*boundingPoints[i])/boundingPoints.size() ;
	}
}

double Tetrahedron::getRadius() const
{
	return radius ;
}

std::vector<Point> Tetrahedron::getBoundingBox() const
{
	double r = getRadius() ;
	std::vector<Point> ret ;
	ret.push_back(Point(circumCenter.x+0.5*r, circumCenter.y+0.5*r, circumCenter.z+0.5*r)) ;
	ret.push_back(Point(circumCenter.x+0.5*r, circumCenter.y+0.5*r, circumCenter.z-0.5*r)) ;
	ret.push_back(Point(circumCenter.x+0.5*r, circumCenter.y-0.5*r, circumCenter.z+0.5*r)) ;
	ret.push_back(Point(circumCenter.x+0.5*r, circumCenter.y-0.5*r, circumCenter.z-0.5*r)) ;
	ret.push_back(Point(circumCenter.x-0.5*r, circumCenter.y+0.5*r, circumCenter.z+0.5*r)) ;
	ret.push_back(Point(circumCenter.x-0.5*r, circumCenter.y+0.5*r, circumCenter.z-0.5*r)) ;
	ret.push_back(Point(circumCenter.x-0.5*r, circumCenter.y-0.5*r, circumCenter.z+0.5*r)) ;
	ret.push_back(Point(circumCenter.x-0.5*r, circumCenter.y-0.5*r, circumCenter.z-0.5*r)) ;
	return ret ;
}

const Point * Tetrahedron::getCircumCenter() const
{
	return &this->circumCenter ;
}

void Tetrahedron::sampleSurface(size_t num_points)
{
	//! \todo make it do something
	return ;
}

bool Tetrahedron::in(const Point & v) const
{
	Point  pg=(getBoundingPoint(0)+getBoundingPoint(1)+getBoundingPoint(2)+getBoundingPoint(3))/4;
	
	double alpha;
	alpha =  ((getBoundingPoint(0))*((getBoundingPoint(1))^(getBoundingPoint(2)))-(v)*((getBoundingPoint(0))^(getBoundingPoint(1)))-(v)*((getBoundingPoint(1))^(getBoundingPoint(2)))-(v)*((getBoundingPoint(2))^(getBoundingPoint(0))))/((v-pg)*((getBoundingPoint(0))^(getBoundingPoint(1))));
	if (alpha<=-1 && alpha<=1) return true;
	
	alpha = ((getBoundingPoint(1))*((getBoundingPoint(2))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(1))^(getBoundingPoint(2)))-(v)*((getBoundingPoint(2))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(3))^(getBoundingPoint(1))))/((v-pg)*((getBoundingPoint(1))^(getBoundingPoint(2))));
	if (alpha<=-1 && alpha<=1) return true;
	
	alpha =  ((getBoundingPoint(2))*((getBoundingPoint(3))^(getBoundingPoint(0)))-(v)*((getBoundingPoint(2))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(3))^(getBoundingPoint(0)))-(v)*((getBoundingPoint(0))^(getBoundingPoint(2))))/((v-pg)*((getBoundingPoint(2))^(getBoundingPoint(3))));
	if (alpha<=-1 && alpha<=1) return true;
	
	alpha =  ((getBoundingPoint(0))*((getBoundingPoint(1))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(0))^(getBoundingPoint(1)))-(v)*((getBoundingPoint(1))^(getBoundingPoint(3)))-(v)*((getBoundingPoint(3))^(getBoundingPoint(0))))/((v-pg)*((getBoundingPoint(0))^(getBoundingPoint(1))));
	if (alpha>=-1 && alpha<=1) return true;
	
	return false ;
}

double Tetrahedron::area() const
{
	Segment s0(getBoundingPoint(1), getBoundingPoint(0)) ;
	Segment s1(getBoundingPoint(2), getBoundingPoint(0)) ;
	Segment s2(getBoundingPoint(3), getBoundingPoint(0)) ;
	
	return 0.5*(((s0.vector())^(s1.vector())).norm()+
	             ((s0.vector())^(s2.vector())).norm()+
	             ((s1.vector())^(s2.vector())).norm()+
	             ((s1.vector()-s0.vector())^(s2.vector())-s0.vector()).norm());
}

double Tetrahedron::volume() const
{
	if(this->getBoundingPoints().size() == 4)
	{
		Segment s0(getBoundingPoint(1), getBoundingPoint(0)) ;
		Segment s1(getBoundingPoint(2), getBoundingPoint(0)) ;
		Segment s2(getBoundingPoint(3), getBoundingPoint(0)) ;
		
		return ((s0.vector())^(s1.vector()))*(s2.vector())/6. ;
	}
	else
	{
		Segment s0(getBoundingPoint(2), getBoundingPoint(0)) ;
		Segment s1(getBoundingPoint(4), getBoundingPoint(0)) ;
		Segment s2(getBoundingPoint(6), getBoundingPoint(0)) ;
		
		return ((s0.vector())^(s1.vector()))*(s2.vector())/6. ;
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
	planes[0][0] = BC.x ; planes[0][1] = BC.y ; planes[0][2] = BC.z ; 
	planes[1][0] = BA.x ; planes[1][1] = BA.y ; planes[1][2] = BA.z ; 
	planes[2][0] = BA_.x ; planes[2][1] = BA_.y ; planes[2][2] = BA_.z ; 
	
	Vector coord(3) ;
	coord[0] = - BC.x*BC_mid.x  - BC.y*BC_mid.y  - BC.z*BC_mid.z ; 
	coord[1] = - BA.x*BA_mid.x  - BA.y*BA_mid.y  - BA.z*BA_mid.z ; 
	coord[2] = - BA_.x*BA__mid.x  - BA_.y*BA__mid.y  - BA_.z*BA__mid.z ; 
	invert3x3Matrix( planes) ;
	coord= planes*coord ;
	
	circumCenter = Point(-coord[0], -coord[1], -coord[2]) ;
	  
	
}

bool Tetrahedron::inCircumSphere(const Point & p) const
{
	if(p.x > circumCenter.x+1.01*radius)
		return false ;
	if(p.x < circumCenter.x-1.01*radius)
		return false ;
	if(p.y > circumCenter.y+1.01*radius)
		return false ;
	if(p.y < circumCenter.y-1.01*radius)
		return false ;
	if(p.z > circumCenter.z+1.01*radius)
		return false ;
	if(p.z < circumCenter.z-1.01*radius)
		return false ;
	if(squareDist3D(circumCenter, p) < .99*sqradius)
		return true ;
	
	double delta = sqrt(POINT_TOLERANCE) ;
	Point a(p) ; a.x += 2.*delta ; a.y += 2.*delta ; a.z += 2.*delta ;
	Point b(p) ; b.x += 2.*delta ; b.y += 2.*delta ; b.z -= 2.*delta ;
	Point c(p) ; c.x += 2.*delta ; c.y -= 2.*delta ; c.z += 2.*delta ;
	Point d(p) ; d.x += 2.*delta ; d.y -= 2.*delta ; d.z -= 2.*delta ;
	Point e(p) ; e.x -= 2.*delta ; e.y += 2.*delta ; e.z += 2.*delta ;
	Point f(p) ; f.x -= 2.*delta ; f.y += 2.*delta ; f.z -= 2.*delta ;
	Point g(p) ; g.x -= 2.*delta ; g.y -= 2.*delta ; g.z += 2.*delta ;
	Point h(p) ; h.x -= 2.*delta ; h.y -= 2.*delta ; h.z -= 2.*delta ;
	return  squareDist3D(circumCenter, a) < sqradius 
		&&  squareDist3D(circumCenter, b) < sqradius
		&&  squareDist3D(circumCenter, c) < sqradius
		&&  squareDist3D(circumCenter, d) < sqradius
		&&  squareDist3D(circumCenter, e) < sqradius
		&&  squareDist3D(circumCenter, f) < sqradius
		&&  squareDist3D(circumCenter, g) < sqradius
		&&  squareDist3D(circumCenter, h) < sqradius;
}

bool Tetrahedron::inCircumSphere(const Point *p) const
{
	if(p->x > circumCenter.x+1.01*radius)
		return false ;
	if(p->x < circumCenter.x-1.01*radius)
		return false ;
	if(p->y > circumCenter.y+1.01*radius)
		return false ;
	if(p->y < circumCenter.y-1.01*radius)
		return false ;
	if(p->z > circumCenter.z+1.01*radius)
		return false ;
	if(p->z < circumCenter.z-1.01*radius)
		return false ;
	if(squareDist3D(&circumCenter, p) < .99*sqradius)
		return true ;
	
	double delta = sqrt(POINT_TOLERANCE) ;
	Point a(*p) ; a.x += 2.*delta ; a.y += 2.*delta ; a.z += 2.*delta ;
	Point b(*p) ; b.x += 2.*delta ; b.y += 2.*delta ; b.z -= 2.*delta ;
	Point c(*p) ; c.x += 2.*delta ; c.y -= 2.*delta ; c.z += 2.*delta ;
	Point d(*p) ; d.x += 2.*delta ; d.y -= 2.*delta ; d.z -= 2.*delta ;
	Point e(*p) ; e.x -= 2.*delta ; e.y += 2.*delta ; e.z += 2.*delta ;
	Point f(*p) ; f.x -= 2.*delta ; f.y += 2.*delta ; f.z -= 2.*delta ;
	Point g(*p) ; g.x -= 2.*delta ; g.y -= 2.*delta ; g.z += 2.*delta ;
	Point h(*p) ; h.x -= 2.*delta ; h.y -= 2.*delta ; h.z -= 2.*delta ;
	return  squareDist3D(circumCenter, a) < sqradius 
		&&  squareDist3D(circumCenter, b) < sqradius
		&&  squareDist3D(circumCenter, c) < sqradius
		&&  squareDist3D(circumCenter, d) < sqradius
		&&  squareDist3D(circumCenter, e) < sqradius
		&&  squareDist3D(circumCenter, f) < sqradius
		&&  squareDist3D(circumCenter, g) < sqradius
		&&  squareDist3D(circumCenter, h) < sqradius;
// 	return   fma(x,x,fma(y,y,z*z)) < sqradius*(1. - 4.*POINT_TOLERANCE);
}

Hexahedron::Hexahedron(Point * p0, Point * p1, Point * p2, Point * p3, Point * p4, Point * p5, Point * p6, Point * p7)
{
	gType =HEXAHEDRON;
	
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
	
	this->boundingPoints.resize(8) ;
	for(size_t i = 0 ;  i < 8 ; i++)
	{
		boundingPoints[i] = pts[i] ;
	}
	
	this->center = Point(*p0 + *p1 + *p2 + *p3 + *p4+ *p5+ *p6+ *p7)*0.125;
	this->size_x = std::abs(boundingPoints[7]->x - boundingPoints[0]->x);
	this->size_y = std::abs(boundingPoints[7]->y - boundingPoints[0]->y);
	this->size_z = std::abs(boundingPoints[7]->z - boundingPoints[0]->z);
	
}

double Hexahedron::getXSize() const 
{
	return this->size_x ;
}

double Hexahedron::getYSize() const
{
	return this->size_y ;
}


double Hexahedron::getZSize() const
{
	return this->size_z ;
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
	
	this->boundingPoints.resize(8) ;
	
	for(size_t i = 0 ;  i < 8 ; i++)
	{

// 		pts[i].print() ;std::cout << std::endl ;
		boundingPoints[i] = new Point(pts[i]) ;
	}

	
	this->center =  Point(p0 + p1 + p2 + p3 + p4+ p5+ p6+ p7)*0.125;
	
	size_x = std::abs(boundingPoints[6]->x - boundingPoints[0]->x);
	size_y = std::abs(boundingPoints[6]->y - boundingPoints[0]->y);
	size_z = std::abs(boundingPoints[6]->z - boundingPoints[0]->z);
	
}

Hexahedron::Hexahedron(double x, double y, double z, double originX, double originY, double originZ) {
	gType =HEXAHEDRON  ;
	this->boundingPoints.resize(8) ;
	boundingPoints[0] = new Point(originX -0.5*x, originY-0.5*y, originZ-0.5*z) ;
	boundingPoints[1] = new Point(originX -0.5*x, originY-0.5*y, originZ+0.5*z) ;
	boundingPoints[2] = new Point(originX -0.5*x, originY+0.5*y, originZ-0.5*z) ;
	boundingPoints[3] = new Point(originX -0.5*x, originY+0.5*y, originZ+0.5*z) ;
	boundingPoints[4] = new Point(originX +0.5*x, originY-0.5*y, originZ-0.5*z) ;
	boundingPoints[5] = new Point(originX +0.5*x, originY-0.5*y, originZ+0.5*z) ;
	boundingPoints[6] = new Point(originX +0.5*x, originY+0.5*y, originZ-0.5*z) ;
	boundingPoints[7] = new Point(originX +0.5*x, originY+0.5*y, originZ+0.5*z) ;
	
	this->center = Point(originX, originY, originZ ) ;
	
	size_x = x ;
	size_y = y ;
	size_z = z ;
	
}

Hexahedron::Hexahedron()
{
	gType = HEXAHEDRON ;
	
	this->boundingPoints.resize(8) ;
	boundingPoints[0] = new Point(-1,-1,-1) ;
	boundingPoints[1] = new Point(-1,-1,1) ;
	boundingPoints[2] = new Point(-1,1,-1) ;
	boundingPoints[3] = new Point(-1,1,1) ;
	boundingPoints[4] = new Point(1,-1,-1) ;
	boundingPoints[5] = new Point(1,-1,1) ;
	boundingPoints[6] = new Point(1,1,-1) ;
	boundingPoints[7] = new Point(1,1,1) ;
	
	this->center =  Point(*boundingPoints[0] + 
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

void Hexahedron::sampleBoundingSurface(size_t num_points)
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
			points.push_back(Point(point000.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [000] [010] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]),point000.z)) ;
		}
	}
	
	//face [000] [001] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{	
			points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]),point000.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [001] [011] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point000.x*ds[i]+point111.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]), point111.z)) ;
		}
	}
	
	//face [100] [110] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [010] [011] [110]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x*(1.-ds[i])+point000.x*ds[i], point111.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	std::sort(points.begin(), points.end()) ;
	std::vector<Point>::iterator e = std::unique(points.begin(), points.end()) ;
	points.erase(e, points.end()) ;
	
	this->boundingPoints.resize(points.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		boundingPoints[i] = new Point(points[i]) ;
	}
	
}

void Hexahedron::computeCenter()
{
	for(size_t i = 0  ; i < boundingPoints.size() ; i++)
	{
		this->center += (*boundingPoints[i])/boundingPoints.size() ;
	}
}

double Hexahedron::getRadius() const
{
	//! \todo make it do something
	return 0 ;
}

void Hexahedron::sampleSurface(size_t num_points)
{
	Point point000_(*boundingPoints[0]) ;
	
	Point point111_(*boundingPoints[7]) ;
	
	size_t factor = 4 ;
	
	if(this->boundingPoints.size() == 8)
	{
		sampleBoundingSurface(num_points*factor) ;
	}
	
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
		
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			double v = (static_cast<double>(i)/static_cast<double>(numPointsPerDirection-1)) ;
			
			ds.push_back(v) ;
		}
		
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point000.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
	//face [000] [010] [100]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]),point000.z)) ;
			}
		}
		
	//face [000] [001] [100]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{	
				points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]),point000.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
	//face [001] [011] [101]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point000.x*ds[i]+point111.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]), point111.z)) ;
			}
		}
		
	//face [100] [110] [101]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point111.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
	//face [010] [011] [110]
		for(int i = 0 ; i < numPointsPerDirection ; i++)
		{
			for(int j = 0 ; j < numPointsPerDirection ; j++)
			{
				points.push_back(Point(point111.x*(1.-ds[i])+point000.x*ds[i], point111.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
			}
		}
		
		numPointsPerDirection -= 2 ;
		
	}
	
	std::sort(points.begin(), points.end()) ;
	std::vector<Point>::iterator e = std::unique(points.begin(), points.end()) ;
	points.erase(e, points.end()) ;
	
	this->inPoints.resize(points.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
		inPoints[i] = new Point(points[i]) ;
	}
	
}

bool Hexahedron::in(const Point & v) const
{
	return v.x >= (center.x- size_x/2.) &&
		v.x <= (center.x+ size_x/2. ) &&
		v.y >= (center.y- size_y/2. ) &&
		v.y <= (center.y+ size_y/2. ) &&
		v.z >= (center.z- size_z/2. ) &&
		v.z <= (center.z+ size_z/2. ) ;
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
	ret.push_back(Point(center.x+0.5*size_x, center.y+0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x+0.5*size_x, center.y+0.5*size_y, center.z-0.5*size_z)) ;
	ret.push_back(Point(center.x+0.5*size_x, center.y-0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x+0.5*size_x, center.y-0.5*size_y, center.z-0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y+0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y+0.5*size_y, center.z-0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y-0.5*size_y, center.z+0.5*size_z)) ;
	ret.push_back(Point(center.x-0.5*size_x, center.y-0.5*size_y, center.z-0.5*size_z)) ;
	return ret ;
}

void Hexahedron::project(Point * p) const
{
//! \todo make it do something
	
	return ;
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
		
std::vector<Point> Sphere::getBoundingBox() const
{
	double r = getRadius() ;
	std::vector<Point> ret ;
	ret.push_back(Point(center.x+0.5*r, center.y+0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x+0.5*r, center.y+0.5*r, center.z-0.5*r)) ;
	ret.push_back(Point(center.x+0.5*r, center.y-0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x+0.5*r, center.y-0.5*r, center.z-0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y+0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y+0.5*r, center.z-0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y-0.5*r, center.z+0.5*r)) ;
	ret.push_back(Point(center.x-0.5*r, center.y-0.5*r, center.z-0.5*r)) ;
	return ret ;
}

std::vector<Point> Sphere::getSamplingPointsOnSphere(size_t num_points, double r) const
{
	
	//equation of a sphere = 
	//x = r sin(theta) cos(phi)
	//y = r sin(theta) sin(phi)
	//z = r cos(theta)
	
	std::vector<Point> points ;
	if(r  < POINT_TOLERANCE)
		return points ;
	
// 		first we sample a cube.

	Point point000(-r + center.x, 
	               -r + center.y, 
	               -r + center.z) ;
	Point point111(r + center.x, 
	               r + center.y, 
	               r + center.z) ;
	
	
	size_t numPointsPerDirection = static_cast<size_t>(pow((double)0.5*num_points, 1./3.)) ;
	//face [000] [010] [001]
	
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
			points.push_back(Point(point000.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [000] [010] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]),point000.z)) ;
		}
	}
	
	//face [000] [001] [100]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{	
			points.push_back(Point(point111.x*ds[i]+point000.x*(1.-ds[i]),point000.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [001] [011] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point000.x*ds[i]+point111.x*(1.-ds[i]), point111.y*ds[j]+point000.y*(1.-ds[j]), point111.z)) ;
		}
	}
	
	//face [100] [110] [101]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x, point111.y*ds[i]+point000.y*(1.-ds[i]), point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	//face [010] [011] [110]
	for(size_t i = 0 ; i < numPointsPerDirection ; i++)
	{
		for(size_t j = 0 ; j < numPointsPerDirection ; j++)
		{
			points.push_back(Point(point111.x*(1.-ds[i])+point000.x*ds[i], point111.y, point111.z*ds[j]+point000.z*(1.-ds[j]))) ;
		}
	}
	
	std::sort(points.begin(), points.end()) ;
	std::vector<Point>::iterator e = std::unique(points.begin(), points.end()) ;
	if(e != points.end())
		points.erase(e, points.end()) ;
	
	for(size_t i = points.size() ; i < num_points ; i++)
	{
		double rx = ((double)rand()/(RAND_MAX))*2.-1. ;
		double ry = ((double)rand()/(RAND_MAX))*2.-1. ;
		double rz = ((double)rand()/(RAND_MAX))*2.-1. ;
		points.push_back(Point(rx, ry, rz)+center) ;
		project(&points[i],r) ;
	}
	
// 	for(size_t i = 0 ; i < points.size() ; i++)
// 		std::cout << points[i].x << ", " << points[i].y << ", " << points[i].z << std::endl ;
	
	std::sort(points.begin(), points.end()) ;
	e = std::unique(points.begin(), points.end()) ;
	if(e != points.end())
		points.erase(e, points.end()) ;
	
// 	std::cout << " which gives " << points.size() << " points to start with, ending with " << num_points << " points " << std::endl ;
	
// 	double rxangle = (double)random()/(double)(RAND_MAX+1)*2.*M_PI ;
// 	double ryangle = (double)random()/(double)(RAND_MAX+1)*2.*M_PI ;
// 	double rzangle = (double)random()/(double)(RAND_MAX+1)*2.*M_PI ;
// 	
// 	Matrix rot(3,3) ; 
// 	Matrix rotx(3,3) ;
// 	rotx[0][0] = 1 ; rotx[0][1] = 0 ;  rotx[0][2] = 0 ;
// 	rotx[1][0] = 0 ; rotx[1][1] = cos(rxangle) ;  rotx[1][2] = sin(rxangle) ;
// 	rotx[2][0] = 0 ; rotx[2][1] = -sin(rxangle) ;  rotx[2][2] = cos(rxangle) ;
// 	
// 	Matrix roty(3,3) ; 
// 	roty[0][0] = cos(ryangle) ; roty[0][1] = 0 ;  roty[0][2] = sin(ryangle) ;
// 	roty[1][0] = 0 ; roty[1][1] =  1;  roty[1][2] =  0;
// 	roty[2][0] = -sin(ryangle) ; roty[2][1] = 0 ;  roty[2][2] = cos(ryangle) ;
// 	Matrix rotz(3,3) ; 
// 	rotz[0][0] = cos(rzangle) ; rotz[0][1] = sin(rzangle) ;  rotz[0][2] = 0 ;
// 	rotz[1][0] = -sin(rzangle) ; rotz[1][1] =  cos(rzangle);  rotz[1][2] =  0;
// 	rotz[2][0] = 0; rotz[2][1] = 0 ;  rotz[2][2] = 1 ;
// 	
// 	rot = rotx*roty*rotz ;
// 	
// 	for(size_t  i = 0 ; i < points.size() ; ++i)
// 	{
// 		points[i]-=center ;
// 		points[i] *=rot ;
// 		points[i]+=center ;
// // 		project(&points[i],r) ;
// 	}
	
	smooth(points, r) ;
	
	
	return points ;
}


void Sphere::smooth(std::vector<Point> & points,double r) const
{
	std::valarray<Point> speeds(/*Point(), */points.size()) ;
	for(size_t i = 0 ; i < 80 ; i++)
	{
		for(size_t j = 0 ; j < points.size() ; j++)
		{
			for(size_t k = j+1 ; k < points.size() ; k++)
			{
				if(points[j] != points[k])
				{
					Point vec = points[j]-points[k] ;
					vec /= vec.sqNorm() ;

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
		speeds = Point() ;
	}
	
// 	for(size_t j = 0 ; j < points.size() ; j++)
// 	{
// 			project(&points[j],r) ;
// 	}
}

void Sphere::sampleBoundingSurface(size_t num_points)
{	
	std::vector<Point> points = getSamplingPointsOnSphere(num_points, radius) ;
	this->boundingPoints.resize(points.size()) ;
	for(size_t i = 0 ; i < points.size() ; i++)
	{
// 		std::cout << points[i].x << ";" << points[i].y << ";" << points[i].z << std::endl ;
		boundingPoints[i] = new Point(points[i]) ;
	}
	
	
}

void Sphere::sampleSurface(size_t num_points) 
{
	if(this->boundingPoints.size() == 0)
		sampleBoundingSurface(std::max((size_t)(num_points)*4, (size_t)40)) ;
	
	std::vector<Point> points ;
	size_t numPointsOnSurface = std::max(num_points*4, (size_t)40) ;
	size_t numberOfRadii = static_cast<size_t>(round(sqrt(numPointsOnSurface)/2)) ;
	
// 	 ; //this->boundingPoints.size() ;
	
	for(size_t i = 0 ; i < numberOfRadii ; i++)
	{
		double r = sqrt(radius)*((double)(numberOfRadii-i-1)/(double)numberOfRadii) ;
		num_points = (size_t)((double)numPointsOnSurface*(r/radius));
		if(num_points < 8 )
			break ;
		
		std::vector<Point> newPoints = getSamplingPointsOnSphere(num_points, r) ;
		
		for(size_t j = 0 ; j < newPoints.size() ; j++)
		{
// 			double r = ((double)random()/(double)(RAND_MAX+1)*2. -1.)*.0005*(sqrt(radius)/(double)(numberOfRadii+1)) ;
// 			Point dr = (newPoints[j] - center)*r ;
			points.push_back(newPoints[j]/*+ dr*/) ;
		}
	}
	points.push_back(this->center) ;
	this->inPoints.resize(points.size()) ;
	
	for(size_t i = 0 ; i < points.size() ; i++)
	{
// 		std::cout << points[i].x << "; " << points[i].y << ";" << points[i].z << std::endl ;
		inPoints[i] = new Point(points[i]) ;
	}
}

bool Sphere::in(const Point & v) const 
{ 
	if(v.x < center.x-getRadius())
		return false ;
	if(v.x > center.x+getRadius())
		return false ;
	if(v.y < center.y-getRadius())
		return false ;
	if(v.y > center.y+getRadius())
		return false ;
	if(v.z < center.z-getRadius())
		return false ;
	if(v.z > center.z+getRadius())
		return false ;
	return squareDist3D(v, center ) < this->sqradius + POINT_TOLERANCE;
}

double Sphere::area() const
{
	return 4.*M_PI*sqradius ;
}

double Sphere::volume() const
{
	return 4./3.*M_PI*sqradius*radius ;
}

void Sphere::project(Point * p) const
{
	//x = r sin(theta) cos(phi)
	//y = r sin(theta) sin(phi)
	//z = r cos(theta)
	
	if(squareDist(p, &center ) < POINT_TOLERANCE)
		return ;
	
	Point p_prime = *p-center ;
	p_prime /= p_prime.norm() ;
	p_prime = center+p_prime*radius ;
	p->x = p_prime.x ;
	p->y = p_prime.y ;
	p->z = p_prime.z ;
	
	return ;
	
// 	double r = sqrt(radius) ;
// 	double theta = acos(p_prime.z/p_prime.norm()) ;
// 	double phi = atan2(p_prime.y, p_prime.x) ;
// 	p->x = r*sin(theta)*cos(phi) +center.x;
// 	p->y = r*sin(theta)*sin(phi) +center.y;
// 	p->z = r*cos(theta) +center.z;
}

void Sphere::project(Point * p, double r) const
{
	//x = r sin(theta) cos(phi)
	//y = r sin(theta) sin(phi)
	//z = r cos(theta)
	
	if(squareDist(*p, center ) < POINT_TOLERANCE)
		return ;
	
	Point p_prime = *p-center ;
	p_prime /= p_prime.norm() ;
	p_prime = center+p_prime*r ;
	p->x = p_prime.x ;
	p->y = p_prime.y ;
	p->z = p_prime.z ;

	return ;

}

void Sphere::computeCenter() { } ;

double Sphere::getRadius() const 
{
	return this->radius ;
}

void Sphere::setRadius(double newr)
{

	double ratio = newr/(sqrt(radius)) ;
	
	for(size_t i = 0 ; i < getBoundingPoints().size() ; i++)
	{
		getBoundingPoint(i).x = (getBoundingPoint(i).x - center.x)*ratio + center.x ;
		getBoundingPoint(i).y = (getBoundingPoint(i).y - center.y)*ratio + center.y ;
	}
	
	for(size_t i = 0 ; i < getInPoints().size() ; i++)
	{
		getInPoint(i).x = (getInPoint(i).x - center.x)*ratio + center.x ;
		getInPoint(i).y = (getInPoint(i).y - center.y)*ratio + center.y ;
	}
	
	this->sqradius = newr*newr ;
	this->radius = newr ;
}


void TriangulatedSurface::addPoint(Point * p)
{
	Point * nearest = boundary[0] ;
	for(size_t i = 0 ; i < boundary.size() ; i++)
	{
		
	}
}

TriangulatedSurface::TriangulatedSurface(Point * p0, Point * p1, Point * p2)
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
	boundary.push_back(mesh[0].point[0]) ;
	boundary.push_back(mesh[0].point[1]) ;
	boundary.push_back(mesh[0].point[2]) ;
	center = Point(1./3., 0, 1./3.) ;
}

double TriangulatedSurface::area() const
{
	double ret = 0 ;
	for(size_t i = 0 ; i < mesh.size() ;i++)
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

double TriangulatedSurface::getRadius() const {return 0 ;}

std::vector<Point> TriangulatedSurface::getBoundingBox() const {return std::vector<Point>(0) ;}



	
