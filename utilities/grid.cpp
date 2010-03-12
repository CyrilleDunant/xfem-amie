#include "grid.h"

using namespace Mu ;

Voxel::Voxel()
{
	tlf.print() ;
	filled = false ;
}

Voxel::Voxel(double x, double y, double z ,double s) : tlf(x-s*.5, y+s*.5, z+s*.5), trf(x+s*.5, y+s*.5, z+s*.5), blf(x-s*.5, y-s*.5, z+s*.5), brf(x+s*.5, y-s*.5, z+s*.5),tlb(x-s*.5, y+s*.5, z-s*.5), trb(x+s*.5, y+s*.5, z-s*.5), blb(x-s*.5, y-s*.5, z-s*.5), brb(x+s*.5, y-s*.5, z-s*.5), filled(false)
{
}

const std::vector<Geometry *> & Voxel::getFeatures() const
{
	return this->features ;
}

std::vector<Geometry *> & Voxel::getFeatures()
{
	return this->features ;
}

bool Voxel::in(const Point & p) const
{
	return (p.x >= tlf.x)  && (p.x <= brf.x) && (p.y >= brf.y) && (p.y <= tlf.y) && (p.z >= brb.z) && (p.z <= tlf.z);
}

void Voxel::coOccuringFeatures(std::vector<Geometry *> &f , const Geometry * inc) const
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

void Voxel::coOccuringFeatures(std::vector<Geometry *> &f , const Point & p) const
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
	return inc->in(tlf) 
		|| inc->in(trf) 
		|| inc->in(brf) 
		|| inc->in(blf) 
		|| inc->in(tlb) 
		|| inc->in(trb) 
		|| inc->in(brb) 
		|| inc->in(blb) 
		|| Hexahedron(tlf.x-brb.x, tlf.y-brb.y, tlf.z-brb.z, (tlf.x+brb.x)*.5, (tlf.y+brb.y)*.5, (tlf.z+brb.z)*.5).intersects(inc) || in(inc->getCenter());
}

void Voxel::remove(Geometry * inc)
{

	std::vector<Geometry *>::iterator e = std::find(features.begin(), features.end(), inc) ;
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
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.25, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.75, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.25, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
		pixels.push_back(new Voxel(tlf.x+(trf.x-tlf.x)*.75, blf.y+(tlf.y-blf.y)*.25, tlb.z+(trf.z-tlb.z)*.75, (trf.x-tlf.x)*.5)) ;
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

bool Voxel::add(Geometry * inc)
{
	
	if(filled)
		return false;
	
	if(!pixels.empty())
	{
		bool ret = true ;
		#pragma omp parallel for
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
		
		if(features.size() > 64 )
			refine() ;
		
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

void Voxel::forceAdd(Geometry * inc)
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
		#pragma omp parallel for
		for(int i = 0 ; i < 8 ; i++)
		{
			if(pixels[i]->coOccur(inc))
				pixels[i]->forceAdd(inc) ;
		}
	}
	
	if(features.size() > 64 )
	{
		refine() ;
	}
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
	return p.x >= tlf.x &&  p.x <= trf.x && p.y >= blf.y &&  p.y <= tlf.y && p.z >= blb.z &&  p.z <= tlf.z;
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
				Point test(blf.x+(double)(i+1)*.1*(brf.x-blf.x), blf.y+(double)(j+1)*.1*(tlf.y-blf.y), blf.y+(double)(k+1)*.1*(tlf.y-blb.y)) ;
				
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

const std::vector<Geometry *> & Pixel::getFeatures() const
{
	return this->features ;
}

std::vector<Geometry *> & Pixel::getFeatures()
{
	return this->features ;
}

bool Pixel::in(const Point & p) const
{
	return (p.x >= tl.x)  && (p.x <= br.x) && (p.y >= br.y) && (p.y <= tl.y);
}

void Pixel::refine()
{
	
	if(pixels.size() == 0)
	{
		pixels.resize(4) ;
		pixels[0] = Pixel(tl.x+(tr.x-tl.x)*.25, bl.y+(tl.y-bl.y)*.75, (tr.x-tl.x)*.5) ;
		pixels[1] = Pixel(tl.x+(tr.x-tl.x)*.75, bl.y+(tl.y-bl.y)*.75, (tr.x-tl.x)*.5) ;
		pixels[2] = Pixel(tl.x+(tr.x-tl.x)*.25, bl.y+(tl.y-bl.y)*.25, (tr.x-tl.x)*.5) ;
		pixels[3] = Pixel(tl.x+(tr.x-tl.x)*.75, bl.y+(tl.y-bl.y)*.25, (tr.x-tl.x)*.5) ;
		
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

bool Pixel::coOccur(const Geometry * inc) const
{
	std::vector<Point> bbox = inc->getBoundingBox() ;
	return inc->in(tl) 
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
		|| in(bbox[3])  ;
}

bool Pixel::coOccur(const Point & p) const
{
	return p.x > tl.x &&  p.x < tr.x && p.y > bl.y &&  p.y < tl.y ;
}

void Pixel::remove(Geometry * inc)
{
	std::vector<Geometry *>::iterator e = std::find(features.begin(), features.end(), inc) ;
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
			Point test(bl.x+(double)(i+1)*.1*(br.x-bl.x), bl.y+(double)(j+1)*.1*(tl.y-bl.y)) ;
			
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

void Pixel::coOccuringFeatures(std::vector<Geometry *> &f , const Geometry * inc) const
{
	if(!pixels.size() && !features.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<Geometry *> ret ;

	
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

void Pixel::coOccuringFeatures(std::vector<Geometry *> &f , const Point & p) const
{
	if(!pixels.size() && !features.empty())
	{
		f.insert(f.end(), features.begin(), features.end()) ;
		return ;
	}
	
	std::vector<Geometry *> ret ;
	
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

bool Pixel::add(Geometry * inc)
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
				return false;
		}
		this->features.push_back(inc) ;
		
		if(features.size() > 32 )
			refine() ;
		
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

void Pixel::forceAdd(Geometry * inc)
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
	

	
	if(!filled && features.size() > 32 )
	{
		refine() ;
	}
}

void Pixel::print() const
{
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
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	else if(y < x && y < z)
	{
		lengthY = div ;
		lengthX = div*(x/y) ;
		lengthZ = div*(z/y) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	else if(z < x && z < y)
	{
		lengthZ = div ;
		lengthX = div*(x/z) ;
		lengthY = div*(y/z) ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	else
	{
		lengthX = div ;
		lengthY = div ;
		lengthZ = div ;
		pixels.resize(lengthX,std::valarray<std::valarray<Voxel *> >(std::valarray<Voxel *>((Voxel *)NULL,lengthZ),lengthY)) ;
	}
	
	psize = x/lengthX;
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			for(size_t k = 0 ; k < lengthZ ; k++)
			{
				
// 				pixels[i][j][k] = new Voxel(x*(double)(i)/(double)lengthX+psize*.5+center.x,
// 				                            y*(double)(j)/(double)lengthY+psize*.5+center.y,
// 				                            z*(double)(k)/(double)lengthZ+psize*.5+center.z,  psize) ;
// 				freepixel.push_back(pixels[i][j][k]) ;
				pixels[i][j][k] = new Voxel(x*(double)(i)/(double)lengthX+psize*.5+center.x-x/2.,
				                            y*(double)(j)/(double)lengthY+psize*.5+center.y-y/2.,
				                            z*(double)(k)/(double)lengthZ+psize*.5+center.z-z/2.,  psize) ;
				unfilledpixel.push_back(pixels[i][j][k]) ;
			}
		}
	}
	
// 	std::sort(freepixel.begin(), freepixel.end()) ;
	std::sort(unfilledpixel.begin(), unfilledpixel.end()) ;
}

Point Grid3D::randomFreeCenter() const 
{
	return Point(x*((double)rand()/(RAND_MAX)-.5)+c.x, 
	             y*((double)rand()/(RAND_MAX)-.5)+c.y, 
	             z*((double)rand()/(RAND_MAX)-.5)+c.z) ;
}

Grid3D Grid3D::getGrid(int div) const
{
	Hexahedron all(x, y, z, c) ;
	std::vector<Geometry *> features = coOccur(&all) ;
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

bool Grid3D::add(Geometry * inc)
{
	
	std::vector<Geometry *> toTest = coOccur(inc);
	for(size_t i = 0 ; i < toTest.size() ; i++)
		if(inc->intersects(toTest[i]))
			return false ;
	
	forceAdd(inc) ;
	
	return true ;
	
// 	bool ret = true ;
// 	std::vector<Voxel *> cleanup ;
// 	
// 	double startX = .5*x + inc->getCenter().x-inc->getRadius() ;
// 	int startI = std::max(0., startX/psize - 2) ;
// 	
// 	double endX =  startX+2.*inc->getRadius();
// 	int endI = std::min(endX/psize + 2, (double)lengthX);
// 	
// 	double startY = .5*y + inc->getCenter().y-inc->getRadius() ;
// 	int startJ = std::max(0., startY/psize - 2) ;
// 	
// 	double endY =  startY+2.*inc->getRadius();
// 	int endJ = std::min(endY/psize + 2, (double)lengthY);
// 	
// 	double startZ = .5*y + inc->getCenter().z-inc->getRadius() ;
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

void Grid3D::forceAdd(Geometry * inc)
{
	double startX = .5*x-c.x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y-c.y + inc->getCenter().y-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = .5*z-c.z + inc->getCenter().z-inc->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*inc->getRadius();
	int endK = std::min(endZ/psize + 2, (double)lengthZ);
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

std::vector<Geometry *> Grid3D::coOccur(const Geometry * geo) const
{
	std::vector<Geometry *> ret ;
	double startX = .5*x-c.x + geo->getCenter().x-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = .5*y-c.y + geo->getCenter().y-geo->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*geo->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = .5*z-c.z + geo->getCenter().z-geo->getRadius() ;
	int startK = std::max(0., startZ/psize - 2) ;
	
	double endZ =  startZ+2.*geo->getRadius();
	int endK = std::min(endZ/psize + 2, (double)lengthZ);

	if(geo->getGeometryType() == TETRAHEDRON)
	{
		const Tetrahedron * t = dynamic_cast<const Tetrahedron *>(geo) ;
		startX = .5*x + t->getCircumCenter()->x-t->getRadius()*1.1 ;
		startI = std::max(0., startX/psize - 2) ;
		
		endX =  startX+2.2*geo->getRadius();
		endI = std::min(endX/psize + 2, (double)lengthX);
		
		startY = .5*y + t->getCircumCenter()->y-t->getRadius()*1.1 ;
		startJ = std::max(0., startY/psize - 2) ;
		
		endY =  startY+2.2*t->getRadius();
		endJ = std::min(endY/psize + 2, (double)lengthY);
		
		startZ = .5*z + t->getCircumCenter()->z-t->getRadius()*1.1 ;
		startK = std::max(0., startZ/psize - 2) ;
		
		endZ =  startZ+2.2*t->getRadius();
		endK = std::min(endZ/psize + 2, (double)lengthZ);
	}
	
	bool foundPixel = false ;
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
					foundPixel = true ;
					pixels[i][j][k]->coOccuringFeatures(ret,geo) ;
					std::stable_sort(ret.begin(), ret.end());
					std::vector<Geometry *>::iterator e = std::unique(ret.begin(), ret.end()) ;
				}
			}
		}
	}
	
	std::stable_sort(ret.begin(), ret.end());
	std::vector<Geometry *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

std::vector<Geometry *> Grid3D::coOccur(const Point & p) const 
{
	std::vector<Geometry *> ret ;
	double startX = x*.5-c.x + p.x ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+.05*x;
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.y + p.y ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+.05*y;
	int endJ = std::min(endY/psize + 2, (double)lengthY);
	
	double startZ = z*.5-c.z + p.z ;
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
	std::vector<Geometry *>::iterator e = std::unique(ret.begin(), ret.end()) ;
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
		
		lengthX = div*(x/y) ;
		lengthY = div ;
		pixels.resize(lengthX,std::valarray<Pixel *>((Pixel *)NULL,lengthY)) ;
	}
	else
	{
		lengthX = div ;
		lengthY = div*(y/x) ;
		pixels.resize(lengthX,std::valarray<Pixel *>((Pixel *)NULL,lengthY)) ;
	}
	
	psize = std::max(std::abs(x/lengthX), std::abs(y/lengthY));
	
	for(size_t i = 0 ; i < lengthX ; i++)
	{
		for(size_t j = 0 ; j < lengthY ; j++)
		{
			pixels[i][j] = new Pixel(x*(double)(i)/(double)lengthX+psize*.5+center.x-.5*x,
			                         y*(double)(j)/(double)lengthY+psize*.5+center.y-.5*y, psize) ;
		}
	}
}

std::vector<Geometry *> Grid::coOccur(const Geometry * geo) const
{
	std::vector<Geometry *> ret ;
	double startX = x*.5-c.x + geo->getCenter().x-geo->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*geo->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.y + geo->getCenter().y-geo->getRadius() ;
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
				std::stable_sort(ret.begin(), ret.end());
				std::vector<Geometry *>::iterator e = std::unique(ret.begin(), ret.end()) ;
				ret.erase(e, ret.end()) ;
			}
		}
	}

	std::stable_sort(ret.begin(), ret.end());
	std::vector<Geometry *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}

 std::vector<Geometry *> Grid::coOccur(const Point & p) const 
{

	std::vector<Geometry *> ret ;
	double startX = x*.5-c.x + p.x ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+.05*x;
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.y + p.y ;
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
	std::vector<Geometry *>::iterator e = std::unique(ret.begin(), ret.end()) ;
	ret.erase(e, ret.end()) ;
	return ret ;
}


void Grid::forceAdd(Geometry * inc)
{
	
	double startX = x*.5-c.x + inc->getCenter().x-inc->getRadius() ;
	int startI = std::max(0., startX/psize - 2) ;
	
	double endX =  startX+2.*inc->getRadius();
	int endI = std::min(endX/psize + 2, (double)lengthX);
	
	double startY = y*.5-c.y + inc->getCenter().y-inc->getRadius() ;
	int startJ = std::max(0., startY/psize - 2) ;
	
	double endY =  startY+2.*inc->getRadius();
	int endJ = std::min(endY/psize + 2, (double)lengthY);

	for(int i = startI ; i < endI ; i++)
	{
		for(int j = startJ ; j < endJ ; j++)
		{
			if(pixels[i][j]->coOccur(inc))
			{
				pixels[i][j]->forceAdd(inc) ;
			}
		}
		
	}


}

Grid Grid::getGrid(int div) const
{
	Rectangle all(x, y, c) ;
	std::vector<Geometry *> features = coOccur(&all) ;
	Grid ret(x,y,div,c) ;
	for(size_t i = 0 ; i < features.size() ; i++)
	{
		ret.forceAdd(features[i]) ;
	}
	
	return ret ;
}

bool Grid::add(Geometry * inc)
{
	std::vector<Geometry *> toTest = coOccur(inc);

	for(size_t i = 0 ; i < toTest.size() ; i++)
	{

		if(inc->intersects(toTest[i]))
			return false ;
	}
	forceAdd(inc) ;
	
	return true ;
	
}
