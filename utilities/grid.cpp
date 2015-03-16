// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011

#include "grid.h"
#include <set>

using namespace Amie ;

Voxel::Voxel()
{
	tlf.print() ;
	filled = false ;
}

Voxel::Voxel(double x, double y, double z ,double s) : tlf(x-s*.5, y+s*.5, z+s*.5), 
                                                       trf(x+s*.5, y+s*.5, z+s*.5), 
                                                       blf(x-s*.5, y-s*.5, z+s*.5), 
                                                       brf(x+s*.5, y-s*.5, z+s*.5),
                                                       tlb(x-s*.5, y+s*.5, z-s*.5), 
                                                       trb(x+s*.5, y+s*.5, z-s*.5), 
                                                       blb(x-s*.5, y-s*.5, z-s*.5), 
                                                       brb(x+s*.5, y-s*.5, z-s*.5), filled(false)
{
}

const std::vector<const Geometry *> & Voxel::getFeatures() const
{
	return this->features ;
}

std::vector<const Geometry *> & Voxel::getFeatures()
{
	return this->features ;
}

bool Voxel::in(const Point & p) const
{
	return (p.getX() >= tlf.getX())  && (p.getX() <= brf.getX()) && (p.getY() >= brf.getY()) && (p.getY() <= tlf.getY()) && (p.getZ() >= brb.getZ()) && (p.getZ() <= tlf.getZ());
}

void Voxel::coOccuringFeatures(std::vector<const Geometry *> &f , const Geometry * inc) const
{
	if(pixels.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	for(size_t i = 0 ; i < 8 ; i++)
	{
		if(pixels[i]->coOccur(inc))
		{
			pixels[i]->coOccuringFeatures(f, inc) ;
		}
	}
}

void Voxel::coOccuringFeatures(std::vector<const Geometry *> &f , const Point & p) const
{
	if(pixels.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<Geometry *> ret ;
	
	for(size_t i = 0 ; i < 8 ; i++)
	{
		if(pixels[i]->coOccur(p))
		{
			pixels[i]->coOccuringFeatures(f, p) ;
		}
	}
}

bool Voxel::coOccur(const Geometry * inc) const
{
    std::vector<Point> bbox = inc->getBoundingBox() ;
    Hexahedron test((tlf.getX()+trf.getX())*.5, (tlf.getY()+blb.getY())*.5, (tlf.getZ()+blb.getZ())*.5, trf.getX()-tlf.getX(), tlb.getY()-blf.getY(), blf.getZ()-tlb.getZ()) ;
    bool ret = inc->in(tlf) 
        || inc->in(trf) 
        || inc->in(brf) 
        || inc->in(blf) 
        || inc->in(tlb) 
        || inc->in(trb) 
        || inc->in(brb) 
        || inc->in(blb)
        || inc->in((trf+blb)*.5)
        || in(inc->getCenter()+Point(inc->getRadius(), 0, 0))
        || in(inc->getCenter())
        || in(inc->getCenter()+Point(-inc->getRadius(), 0, 0)) 
        || in(inc->getCenter()+Point(0,inc->getRadius(), 0)) 
        || in(inc->getCenter()+Point(0,-inc->getRadius(), 0))
        || in(inc->getCenter()+Point(0, 0,inc->getRadius())) 
        || in(inc->getCenter()+Point(0, 0,-inc->getRadius()))
        || in(bbox[0])
        || in(bbox[1])
        || in(bbox[2]) 
        || in(bbox[3])  
        || in(bbox[4])
        || in(bbox[5])
        || in(bbox[6]) 
        || in(bbox[7])  
        || test.intersects(inc);
    return ret ;
    
    
    
	if(inc->getGeometryType() == TETRAHEDRON)
	{
		const Tetrahedron * t = dynamic_cast<const Tetrahedron *>(inc) ;
		size_t n = t->getBoundingPoints().size() ;
/*		t->getBoundingPoint(0).print() ;
		t->getBoundingPoint(1).print() ;
		t->getBoundingPoint(2).print() ;
		t->getBoundingPoint(3).print() ;
		tlf.print() ;
		trf.print() ;
		blf.print() ;
		brf.print() ;
		tlb.print() ;
		trb.print() ;
		blb.print() ;
		brb.print() ;
		std::cout << std::endl ;*/
		return /*inc->in(tlf) 
		|| inc->in(trf) 
		|| inc->in(brf) 
		|| inc->in(blf) 
		|| inc->in(tlb) 
		|| inc->in(trb) 
		|| inc->in(brb) 
		|| inc->in(blb) 
		||*/ in(inc->getCenter())
		|| in(t->getBoundingPoint(0)) 
		|| in(t->getBoundingPoint(n/4)) 
		|| in(t->getBoundingPoint(n*2/4)) 
		|| in(t->getBoundingPoint(n*3/4)) ;
	}

#warning this requires Hexahedron-tetraheron implementation
	return inc->in(tlf) 
		|| inc->in(trf) 
		|| inc->in(brf) 
		|| inc->in(blf) 
		|| inc->in(tlb) 
		|| inc->in(trb) 
		|| inc->in(brb) 
		|| inc->in(blb) 
		|| /*Hexahedron(tlf.getX()-brb.getX(), tlf.getY()-brb.getY(), tlf.getZ()-brb.getZ(), (tlf.getX()+brb.getX())*.5, (tlf.getY()+brb.getY())*.5, (tlf.getZ()+brb.getZ())*.5).intersects(inc) ||*/ in(inc->getCenter()) ;
}

void Voxel::remove(const Geometry * inc)
{

	auto e = std::find(features.begin(), features.end(), inc) ;
	if(e != features.end())
		features.erase(e) ;
	
	if(!pixels.empty())
		for(size_t i = 0 ; i < 8 ; i++)
			pixels[i]->remove(inc) ;
	
	if(features.empty() && !pixels.empty())
	{
		pixels.clear() ;
	}
	
	filled = false ;
}

bool Voxel::isFilled() const
{
	return filled ;
}

Point Voxel::center() const
{
	return (tlf + brb)*.5 ;
}

void Voxel::refine()
{
	if(pixels.empty())
	{
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.25, blf.getY()+(tlf.getY()-blf.getY())*.75, tlb.getZ()+(trf.getZ()-tlb.getZ())*.25, (trf.getX()-tlf.getX())*.5)) ;
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.75, blf.getY()+(tlf.getY()-blf.getY())*.75, tlb.getZ()+(trf.getZ()-tlb.getZ())*.25, (trf.getX()-tlf.getX())*.5)) ;
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.25, blf.getY()+(tlf.getY()-blf.getY())*.25, tlb.getZ()+(trf.getZ()-tlb.getZ())*.25, (trf.getX()-tlf.getX())*.5)) ;
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.75, blf.getY()+(tlf.getY()-blf.getY())*.25, tlb.getZ()+(trf.getZ()-tlb.getZ())*.25, (trf.getX()-tlf.getX())*.5)) ;
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.25, blf.getY()+(tlf.getY()-blf.getY())*.75, tlb.getZ()+(trf.getZ()-tlb.getZ())*.75, (trf.getX()-tlf.getX())*.5)) ;
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.75, blf.getY()+(tlf.getY()-blf.getY())*.75, tlb.getZ()+(trf.getZ()-tlb.getZ())*.75, (trf.getX()-tlf.getX())*.5)) ;
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.25, blf.getY()+(tlf.getY()-blf.getY())*.25, tlb.getZ()+(trf.getZ()-tlb.getZ())*.75, (trf.getX()-tlf.getX())*.5)) ;
		pixels.push_back(new Voxel(tlf.getX()+(trf.getX()-tlf.getX())*.75, blf.getY()+(tlf.getY()-blf.getY())*.25, tlb.getZ()+(trf.getZ()-tlb.getZ())*.75, (trf.getX()-tlf.getX())*.5)) ;
	}
	else
		return ;
	
	for(size_t i = 0 ; i < 8 ; i++)
	{
		for(size_t j = 0 ; j < features.size() ; j++)
			if(pixels[i]->coOccur(features[j]))
				pixels[i]->forceAdd(features[j]) ;
	}

	features.clear() ;
}

bool Voxel::add(const Geometry * inc)
{
	
	if(filled)
		return false;
	
	if(!pixels.empty())
	{
		bool ret = true ;
		#pragma omp parallel for schedule(runtime)
		for(int i = 0 ; i < 8 ; i++)
			if(pixels[i]->coOccur(inc))
				ret = ret && pixels[i]->add(inc) ;
		
		if(!ret)
			for(size_t i = 0 ; i < 8 ; i++)
				pixels[i]->remove(inc) ;
		
		if(ret)
			forceAdd(inc) ;
		return ret ;
	}
	
	if(!features.empty())
	{
		for(size_t i = 0 ; i < this->features.size() ; i++)
		{
			if(this->features[i]->intersects(inc))
				return false;
		}
		
		this->features.push_back(inc) ;
		
// 		if(features.size() > 64 )
// 			refine() ;
// 		
		return true;
	}
	else
	{
		if(inc->in(tlf) 
			&& inc->in(trf) 
			&& inc->in(brf) 
			&& inc->in(blf)
			&& inc->in(tlb) 
			&& inc->in(trb) 
			&& inc->in(brb) 
			&& inc->in(blb)
			)
			filled = true ;
		
		this->features.push_back(inc) ;
		
		return true ;
	}
}

void Voxel::forceAdd(const Geometry * inc)
{
	if(pixels.empty())
		this->features.push_back(inc) ;
	
	if(pixels.empty())
		if(inc->in(tlf) 
		&& inc->in(trf) 
		&& inc->in(brf) 
		&& inc->in(blf)
		&& inc->in(tlb) 
		&& inc->in(trb) 
		&& inc->in(brb) 
		&& inc->in(blb)
		)
			filled = true ;
	
	if(!pixels.empty())
	{
//		#pragma omp parallel for
		for(int i = 0 ; i < 8 ; i++)
		{
			if(pixels[i]->coOccur(inc))
				pixels[i]->forceAdd(inc) ;
		}
	}
	
// 	if(features.size() > 64 )
// 	{
// 		refine() ;
// 	}
}

void Voxel::print() const
{
// 	for(size_t i = 0 ; i < features.size() ; i++)
// 		features[i]->print() ;
}

Voxel::~Voxel()
{
	if(!pixels.empty())
		for(size_t i = 0 ; i < 8 ; i++)
			delete pixels[i] ;
}

bool Voxel::coOccur(const Point & p) const
{
	return p.getX() >= tlf.getX()-POINT_TOLERANCE 
	   &&  p.getX() <= trf.getX()+POINT_TOLERANCE 
	   &&  p.getY() >= blf.getY()-POINT_TOLERANCE 
	   &&  p.getY() <= tlf.getY()+POINT_TOLERANCE 
	   &&  p.getZ() >= blb.getZ()-POINT_TOLERANCE 
	   &&  p.getZ() <= tlf.getZ()+POINT_TOLERANCE;
}

int Voxel::computeFillFactor() const
{
	if(!pixels.empty())
		return 0 ;
	
	int count  = 0 ;
	for(size_t i = 0 ; i < 10 ; i++)
	{
		for(size_t j = 0 ; j < 10 ; j++)
		{
			for(size_t k = 0 ; k < 10 ; k++)
			{
				Point test(blf.getX()+(double)(i+1)*.1*(brf.getX()-blf.getX()), blf.getY()+(double)(j+1)*.1*(tlf.getY()-blf.getY()), blf.getY()+(double)(k+1)*.1*(tlf.getY()-blb.getY())) ;
				
				for(size_t l = 0 ; l < features.size() ; l++)
				{
					if(features[l]->in(test))
					{
						count++ ;
						break ;
					}
				}
			}
		}
	}
	
	return count ;
}

Pixel::Pixel() : pixels(0)
{
	filled = false ;
}

Pixel::Pixel(double x, double y, double s) : tl(x-s*.5, y+s*.5), tr(x+s*.5, y+s*.5), bl(x-s*.5, y-s*.5), br(x+s*.5, y-s*.5), filled(false) {} ;

const std::vector<const Geometry *> & Pixel::getFeatures() const
{
	return this->features ;
}

std::vector<const Geometry *> & Pixel::getFeatures()
{
	return this->features ;
}

bool Pixel::in(const Point & p) const
{
	return (p.getX() >= tl.getX())  && (p.getX() <= br.getX()) && (p.getY() >= br.getY()) && (p.getY() <= tl.getY());
}

void Pixel::refine()
{
	
	if(pixels.size() == 0)
	{
		pixels.resize(4) ;
		pixels[0] = Pixel(tl.getX()+(tr.getX()-tl.getX())*.25, bl.getY()+(tl.getY()-bl.getY())*.75, (tr.getX()-tl.getX())*.5) ;
		pixels[1] = Pixel(tl.getX()+(tr.getX()-tl.getX())*.75, bl.getY()+(tl.getY()-bl.getY())*.75, (tr.getX()-tl.getX())*.5) ;
		pixels[2] = Pixel(tl.getX()+(tr.getX()-tl.getX())*.25, bl.getY()+(tl.getY()-bl.getY())*.25, (tr.getX()-tl.getX())*.5) ;
		pixels[3] = Pixel(tl.getX()+(tr.getX()-tl.getX())*.75, bl.getY()+(tl.getY()-bl.getY())*.25, (tr.getX()-tl.getX())*.5) ;
		
		for(size_t i = 0 ; i < 4 ; i++)
		{
			for(size_t j = 0 ; j < features.size() ; j++)
				if(pixels[i].coOccur(features[j]))
					pixels[i].forceAdd(features[j]) ;
		}
	}
	else
		return ;
}

Pixel::~Pixel()
{
}

bool Pixel::coOccur(const Geometry * const inc) const
{
	std::vector<Point> bbox = inc->getBoundingBox() ;
	Rectangle test((tl.getX()+br.getX())*.5, (tl.getY()+br.getY())*.5, tr.getX()-bl.getX(), tr.getY()-bl.getY()) ;
	bool ret = inc->in(tl) 
		|| inc->in(tr) 
		|| inc->in(br) 
		|| inc->in(bl) 
		|| inc->in(tr*.5+bl*.5)
		|| in(inc->getCenter()+Point(inc->getRadius(), 0))
		|| in(inc->getCenter())
		|| in(inc->getCenter()+Point(-inc->getRadius(), 0)) 
		|| in(inc->getCenter()+Point(0,inc->getRadius())) 
		|| in(inc->getCenter()+Point(0,-inc->getRadius()))
		|| in(bbox[0])
		|| in(bbox[1])
		|| in(bbox[2]) 
		|| in(bbox[3])  
		|| test.intersects(inc);
	return ret ;
}

bool Pixel::coOccur(const Point & p) const
{
	return p.getX() >= tl.getX() - POINT_TOLERANCE 
	   &&  p.getX() <= tr.getX() + POINT_TOLERANCE 
	   &&  p.getY() >= bl.getY() - POINT_TOLERANCE 
	   &&  p.getY() <= tl.getY() + POINT_TOLERANCE ;
}

void Pixel::remove(const Geometry * inc)
{
	auto e = std::find(features.begin(), features.end(), inc) ;
	if(e != features.end())
		features.erase(e) ;
	
	if(pixels.size())
		for(size_t i = 0 ; i < 4 ; i++)
			pixels[i].remove(inc) ;
	
	if(features.empty() && pixels.size())
	{
		pixels.resize(0) ;
	}
	
	filled = false ;
}

int Pixel::computeFillFactor() const
{
	if(!pixels.size())
		return 0 ;
	
	int count  = 0 ;
	for(size_t i = 0 ; i < 10 ; i++)
	{
		for(size_t j = 0 ; j < 10 ; j++)
		{
			Point test(bl.getX()+(double)(i+1)*.1*(br.getX()-bl.getX()), bl.getY()+(double)(j+1)*.1*(tl.getY()-bl.getY())) ;
			
			for(size_t k = 0 ; k < features.size() ; k++)
			{
				if(features[k]->in(test))
				{
					count++ ;
					break ;
				}
			}
		}
	}
	
	return count ;
}

void Pixel::coOccuringFeatures(std::vector<const Geometry *> &f , const Geometry * inc) const
{
	if(!pixels.size() && !features.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<const Geometry *> ret ;

	
	if(pixels.size())
	{
		for(size_t i = 0 ; i < 4 ; i++)
		{
			if(pixels[i].coOccur(inc))
			{
				pixels[i].coOccuringFeatures(f, inc) ;
			}
		}
	}
}

void Pixel::coOccuringFeatures(std::vector<const Geometry *> &f , const Point & p) const
{
	if(!pixels.size() && !features.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<const Geometry *> ret ;
	
	if(pixels.size())
	{
		for(size_t i = 0 ; i < 4 ; i++)
		{
			if(pixels[i].coOccur(p))
			{
				pixels[i].coOccuringFeatures(f, p) ;
			}
		}
	}
}

bool Pixel::add(const Geometry * inc)
{
	if(filled)
		return false;
	
	if(pixels.size())
	{
		bool ret = true ;
		for(size_t i = 0 ; i < 4 ; i++)
		{
			if(pixels[i].coOccur(inc))
				ret = ret && pixels[i].add(inc) ;
		}
		
		if(!ret)
			for(size_t i = 0 ; i < 4 ; i++)
				pixels[i].remove(inc) ;
		
		if(ret)
			forceAdd(inc) ;
		
		return ret ;
	}

	if(!features.empty())
	{
		for(size_t i = 0 ; i < this->features.size() ; i++)
		{
			if(this->features[i]->intersects(inc))
			{
				return false;
			}
		}
		this->features.push_back(inc) ;
		
// 		if(features.size() > 32 )
// 			refine() ;
		
		return true;
	}
	else
	{
		if(inc->in(tl) && inc->in(tr) && inc->in(br) && inc->in(bl))
			filled = true ;
		
		this->features.push_back(inc) ;
		
		return true ;
	}
}

void Pixel::forceAdd(const Geometry * inc)
{	
	this->features.push_back(inc) ;
	
	if(inc->in(tl) && inc->in(tr) && inc->in(br) && inc->in(bl))
		filled = true ;
	
	if(pixels.size())
	{
		for(size_t i = 0 ; i < 4 ; i++)
			if(pixels[i].coOccur(inc))
				pixels[i].forceAdd(inc) ;
	}
	

	
// 	if(!filled && features.size() > 32 )
// 	{
// 		refine() ;
// 	}
}

void Pixel::print() const
{
	tl.print() ;
	tr.print() ;
	bl.print() ;
	br.print() ;
// 	for(size_t i = 0 ; i < features.size() ; i++)
// 		features[i]->print() ;
}


Grid3D::Grid3D(double sizeX, double sizeY, double sizeZ, int div, const Point & center ): x(std::abs(sizeX)), y(std::abs(sizeY)) , z(std::abs(sizeZ)), c(center)
{
	dirtyCounter = 0;
	if(x < y && x < z)
	{
		lengthX = div ;
		lengthY = div*(y/x) ;
		lengthZ = div*(z/x) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)nullptr,lengthZ),lengthY)) ;
	}
	else if(y < x && y < z)
	{
		lengthY = div ;
		lengthX = div*(x/y) ;
		lengthZ = div*(z/y) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)nullptr,lengthZ),lengthY)) ;
	}
	else if(z < x && z < y)
	{
		lengthZ = div ;
		lengthX = div*(x/z) ;
		lengthY = div*(y/z) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)nullptr,lengthZ),lengthY)) ;
	}
	else
	{
		lengthX = div ;
		lengthY = div ;
		lengthZ = div ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)nullptr,lengthZ),lengthY)) ;
	}
	
	psize = x/lengthX;
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			for(size_t k = 0 ; k < lengthZ ; k++)
			{
				
// 				pixels[i][j][k] = new Voxel(x*(double)(i)/(double)lengthX+psize*.5+center.getX(),
// 				                            y*(double)(j)/(double)lengthY+psize*.5+center.getY(),
// 				                            z*(double)(k)/(double)lengthZ+psize*.5+center.getZ(),  psize) ;
// 				freepixel.push_back(pixels[i][j][k]) ;
				pixels[i][j][k] = new Voxel(x*(double)(i)/(double)lengthX+psize*.5+center.getX()-x/2.,
				                            y*(double)(j)/(double)lengthY+psize*.5+center.getY()-y/2.,
				                            z*(double)(k)/(double)lengthZ+psize*.5+center.getZ()-z/2.,  psize) ;
				unfilledpixel.push_back(pixels[i][j][k]) ;
			}
		}
	}
	
// 	std::sort(freepixel.begin(), freepixel.end()) ;
	std::sort(unfilledpixel.begin(), unfilledpixel.end()) ;
}

Point Grid3D::randomFreeCenter() const 
{
	return Point(x*((double)rand()/(RAND_MAX)-.5)+c.getX(), 
	             y*((double)rand()/(RAND_MAX)-.5)+c.getY(), 
	             z*((double)rand()/(RAND_MAX)-.5)+c.getZ()) ;
}

Grid3D Grid3D::getGrid(int div) const
{
	Hexahedron all(x, y, z, c) ;
	std::vector<const Geometry *> features = coOccur(&all) ;
	Grid3D ret(x,y, z,div,c) ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		ret.forceAdd(features[i]) ;
	}
	
	return ret ;
}

// std::vector<HexahedralElement *> elementsFromGrid(int div)
// {
// 	Grid3D newGrid(div) ;
// }

Grid3D::~Grid3D()
{
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			for(size_t k = 0 ; k < lengthZ ; k++)
			{
				delete pixels[i][j][k] ;
			}
		}
	}
}

double Grid3D::fraction() const
{
	return (double)dirtyCounter/(lengthX*lengthY*lengthZ) ;
}

bool Grid3D::add(const Geometry * inc)
{
	
	std::vector<const Geometry *> toTest = coOccur(inc);
	for(size_t i = 0 ; i < toTest.size() ; i++)
		if(
            inc->intersects(toTest[i])
            || inc->in(toTest[i]->getCenter())
            || toTest[i]->in(inc->getCenter())
        )
			return false ;
	
	forceAdd(inc) ;
	
	return true ;
	
// 	bool ret = true ;
// 	std::vector<Voxel *> cleanup ;
// 	
// 	double startX = .5*x + inc->getCenter().getX()-inc->getRadius() ;
// 	int startI = std::max(0., startX/psize - 2) ;
// 	
// 	double endX =  startX+2.*inc->getRadius();
// 	int endI = std::min(endX/psize + 2, (double)lengthX);
// 	
// 	double startY = .5*y + inc->getCenter().getY()-inc->getRadius() ;
// 	int startJ = std::max(0., startY/psize - 2) ;
// 	
// 	double endY =  startY+2.*inc->getRadius();
// 	int endJ = std::min(endY/psize + 2, (double)lengthY);
// 	
// 	double startZ = .5*y + inc->getCenter().getZ()-inc->getRadius() ;
// 	int startK = std::max(0., startZ/psize - 2) ;
// 	
// 	double endZ =  startZ+2.*inc->getRadius();
// 	int endK = std::min(endZ/psize + 2, (double)lengthZ);
// 	
// 	for(int i = startI ; i < endI ; i++)
// 	{
// 		for(int j = startJ ; j < endJ ; j++)
// 		{
// 			for(int k = startK ; k < endK ; k++)
// 			{
// 
// 				if(pixels[i][j][k]->coOccur(inc))
// 				{
// 					if(pixels[i][j][k]->add(inc))
// 					{
// 						cleanup.push_back(pixels[i][j][k]) ;
// 					}
// 					else
// 					{
// 						for(size_t l = 0 ; l < cleanup.size() ; l++)
// 						{
// 							cleanup[l]->remove(inc) ;
// 						}
// 						return false ;
// 					}
// 					
// 				}
// 			}
// 		}
// 		
// 	}
// 	
// 	return ret ;
}

void Grid3D::forceAdd(const Geometry * inc)
{
	double startX = .5*x-c.getX() + inc->getCenter().getX()-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)pixels.size());
	
	double startY = .5*y-c.getY() + inc->getCenter().getY()-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)pixels[0].size());
	
	double startZ = .5*z-c.getZ() + inc->getCenter().getZ()-inc->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*inc->getRadius();
	int endK = std::min(endZ/psize + 2, (double)pixels[0][0].size());
	std::vector<Voxel *> cleanup ;
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
			{
				if(pixels[i][j][k]->coOccur(inc))
				{
					pixels[i][j][k]->forceAdd(inc) ;
				}
			}
		}
	}

}

bool Grid3D::remove(const Geometry* inc)
{

	for(size_t i = 0 ; i < pixels.size() ; i++)
	{
		for(size_t j = 0 ; j < pixels[i].size() ; j++)
		{
			for(size_t k = 0 ; k < pixels[i][j].size() ; k++)
			{
					pixels[i][j][k]->remove(inc) ;
			}
		}
	}

	return true ;
}

std::vector<const Geometry *> Grid3D::coOccur(const Geometry * geo) const
{
	std::vector<const Geometry *> ret ;
	double startX = .5*x-c.getX() + geo->getCenter().getX()-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)pixels.size());
	
	double startY = .5*y-c.getY() + geo->getCenter().getY()-geo->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/psize + 2, (double)pixels[0].size());
	
	double startZ = .5*z-c.getZ() + geo->getCenter().getZ()-geo->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*geo->getRadius();
	int endK = std::min(endZ/psize + 2, (double)pixels[0][0].size());

/*	if(geo->getGeometryType() == TETRAHEDRON)
	{
		const Tetrahedron * t = dynamic_cast<const Tetrahedron *>(geo) ;
		startX = .5*x + t->getCircumCenter().getX()-t->getRadius()*1.1 ;
		startI = std::max(0., startX/psize - 2) ;
		
		endX =  startX+2.2*geo->getRadius();
		endI = std::min(endX/psize + 2, (double)pixels.size());
		
		startY = .5*y + t->getCircumCenter().getY()-t->getRadius()*1.1 ;
		startJ = std::max(0., startY/psize - 2) ;
		
		endY =  startY+2.2*t->getRadius();
		endJ = std::min(endY/psize + 2, (double)pixels[0].size());
		
		startZ = .5*z + t->getCircumCenter().getZ()-t->getRadius()*1.1 ;
		startK = std::max(0., startZ/psize - 2) ;
		
		endZ =  startZ+2.2*t->getRadius();
		endK = std::min(endZ/psize + 2, (double)pixels[0][0].size());
	}*/
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
			{
// 	for(int i = 0 ; i < lengthX ; i++)
// 	{
// 		for(int j = 0 ; j < lengthY ; j++)
// 		{
// 			for(int k = 0 ; k < lengthZ ; k++)
// 			{
				if(pixels[i][j][k]->coOccur(geo))
				{
					pixels[i][j][k]->coOccuringFeatures(ret,geo) ;
				}
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	auto e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

std::vector<const Geometry *> Grid3D::coOccur(const Point & p) const 
{
	std::vector<const Geometry *> ret ;
	double startX = x*.5-c.getX() + p.getX() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+.05*x;
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.getY() + p.getY() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+.05*y;
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = z*.5-c.getZ() + p.getZ() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+.05*z;
	int endK = std::min(endZ/psize + 2, (double)lengthZ);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			for(int k = startK ; k < endK ; k++)
			{
				if(pixels[i][j][k]->coOccur(p))
				{
					pixels[i][j][k]->coOccuringFeatures(ret,p) ;
				}
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	auto e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;

}

Grid::~Grid()
{
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			delete pixels[i][j] ;
		}
	}
}

Grid::Grid(double sizeX, double sizeY, int div, const Point & center ) : x(sizeX), y(sizeY), c(center)
{
	if(x>y)
	{
		
		lengthX = div ;
		lengthY = std::max((int)round(div*y/x), 1) ;
		pixels.resize(lengthX,std::valarray<Pixel *>((Pixel *)nullptr,lengthY)) ;
	}
	else
	{
		lengthX = std::max((int)round(div*x/y), 1) ;
		lengthY = div ;
		pixels.resize(lengthX,std::valarray<Pixel *>((Pixel *)nullptr,lengthY)) ;
	}
	
	psize = std::max(std::abs(x/lengthX), std::abs(y/lengthY));

	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			this->pixels[i][j] = new Pixel(x*(double)(i)/(double)lengthX+psize*.5+center.getX()-.5*x,
			                         y*(double)(j)/(double)lengthY+psize*.5+center.getY()-.5*y, psize) ;
		}
	}
}

std::vector<const Geometry *> Grid::coOccur(const Geometry * geo) const
{
	std::vector<const Geometry *> ret ;
	double startX = x*.5-c.getX() + geo->getCenter().getX()-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.getY() + geo->getCenter().getY()-geo->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
		
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(geo))
			{
				pixels[i][j]->coOccuringFeatures(ret,geo) ;
			}
		}
	}

	std::stable_sort(ret.begin(), ret.end());
	auto e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

 std::vector<const Geometry *> Grid::coOccur(const Point & p) const 
{
	std::vector<const Geometry *> ret ;
	double startX = x*.5-c.getX() + p.getX() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+.05*x;
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.getY() + p.getY() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+.05*y;
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(p))
			{
				pixels[i][j]->coOccuringFeatures(ret,p) ;
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	auto e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}


void Grid::forceAdd(const Geometry * inc)
{
	
	double startX = x*.5-c.getX() + inc->getCenter().getX()-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.getY() + inc->getCenter().getY()-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);

//	std::cout << "grid force add: " << psize << "\" << startI << "," << endI << "\t" << startJ << "," << endJ << std::endl ;
	
	bool done = false ;
	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(inc))
			{
				pixels[i][j]->forceAdd(inc) ;
				done = true ;
			}
		}
	}
	
	
	if(!done)
		pixels[0][0]->forceAdd(inc) ;
}

bool Grid::remove(const Geometry * inc)
{
	for(size_t i = 0 ; i < pixels.size() ; i++)
	{
		for(size_t j = 0 ; j < pixels[i].size() ; j++)
		{
				pixels[i][j]->remove(inc) ;
		}
	}
	
	return true ;
}

Grid Grid::getGrid(int div) const
{
	Rectangle all(x, y, c) ;
	std::vector<const Geometry *> features = coOccur(&all) ;
	Grid ret(x,y,div,c) ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		ret.forceAdd(features[i]) ;
	}
	
	return ret ;
}



bool Grid::add(const Geometry * inc)
{
	std::vector<const Geometry *> toTest = coOccur(inc);

	for(size_t i = 0 ; i < toTest.size() ; i++)
	{
		if(inc->intersects(toTest[i]))
			return false ;
		if(toTest[i]->in(inc->getCenter()))
			return false ;
	}
	forceAdd(inc) ;
	
	return true ;
	
}
