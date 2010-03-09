#ifndef GRID_H
#define GRID_H
#include<valarray>

#include "../geometry/geometry_base.h"
#include "../geometry/geometry_2D.h"
#include "../geometry/geometry_3D.h"

namespace Mu
{
	
	/** \brief Member of an acess grid. Contains references to finer grid level and/or Feature s*/
class Pixel
{
protected:
	std::vector<Geometry *> features ;
	Point tl ;
	Point tr ;
	Point bl ;
	Point br ;
	bool filled ;
	std::valarray<Pixel> pixels ;
	void refine() ;
	int computeFillFactor() const;
// 	const short level ;
// 	const short levels ;
public:
	/** \brief default constructor*/
	Pixel();
	~Pixel();
	/** \brief Construct a pixel centred at x,y of size s
	 *
	 * @param x center x
	 * @param y center y
	 * @param s size
	*/
	Pixel(double x, double y, double s) ;

	/** \brief return features stored in this pixel*/
	const std::vector<Geometry *> & getFeatures() const;
	
	/** \brief return features stored in this pixel*/
	std::vector<Geometry *> & getFeatures();
	
	/** \brief return true if the argument lies in the pixel*/
	bool in(const Point & p) const;

	/** \brief return true if the Geometry overlaps the pixel*/
	bool coOccur(const Geometry * inc) const;
	
	/** \brief return true if the argument lies in the pixel*/
	bool coOccur(const Point & p) const;

	/** \brief Given a Geometry, return a list of stored Feature pointers overlapping the Geometry
*
* The search goes recursively through all the sub-pixels
* @param  ret Vector in which the result should be stored
* @param inc Geometry to check for overlaps
*/
	void coOccuringFeatures(std::vector<Geometry *>& ret, const Geometry * inc) const ;

	/** \brief Given a Point, return a list of stored Feature pointers in which the point lies
*
* The search goes recursively through all the sub-pixels
* @param  ret Vector in which the result should be stored
* @param inc Point to check for overlaps
*/
	void coOccuringFeatures(std::vector<Geometry *>& ret , const Point & p) const ;

/** \brief remove the argument from the Feature list*/
	void remove(Geometry * inc);
	
/** \brief add the argument to the Feature list if it does not overlap with another already present Feature*/
	bool add(Geometry * inc);

/** \brief add the argument to the Feature list unconditionnally*/
	void forceAdd(Geometry * inc) ;

	void print() const ;

} ;


/** \brief Member of an acess grid. Contains references to finer grid level and/or Feature s*/
class Voxel
{
protected:
	std::vector<Geometry *> features ;
	Point tlf ;
	Point trf ;
	Point blf ;
	Point brf ;
	Point tlb ;
	Point trb;
	Point blb ;
	Point brb ;
	bool filled ;

	std::vector<Voxel *> pixels ;
	void refine() ;
	int computeFillFactor() const;
// 	const short level ;
// 	const short levels ;
public:
	/** \brief default constructor*/
	Voxel();

		/** \brief Construct a voxel centred at x,y,z of size s
	 *
	 * @param x center x
	 * @param y center y
	 * @param z center z
	 * @param s size
	*/
	Voxel(double x, double y, double z ,double s) ;

	~Voxel();

	/** \brief return features stored in this voxel*/
	const std::vector<Geometry *> & getFeatures() const;
	
	/** \brief return features stored in this voxel*/
	std::vector<Geometry *> & getFeatures();
	
	/** \brief return true if the argument lies in the voxel*/
	bool in(const Point & p) const;

	/** \brief return true if the Geometry overlaps the voxel*/
	bool coOccur(const Geometry * inc) const;

	/** \brief return true if the argument lies in the voxel*/
	bool coOccur(const Point & p) const;

	/** \brief Given a Geometry, return a list of stored Feature pointers overlapping the Geometry
*
* The search goes recursively through all the sub-voxels
* @param  ret Vector in which the result should be stored
* @param inc Geometry to check for overlaps
*/
	void coOccuringFeatures(std::vector<Geometry *>&, const Geometry * inc) const ;

	/** \brief Given a Point, return a list of stored Feature pointers in which the point lies
*
* The search goes recursively through all the sub-pixels
* @param  ret Vector in which the result should be stored
* @param inc Point to check for overlaps
*/
	void coOccuringFeatures(std::vector<Geometry *>&, const Point & p) const ;

/** \brief remove the argument from the Feature list*/
	void remove(Geometry * inc);
	
/** \brief add the argument to the Feature list if it does not overlap with another already present Feature*/
	bool add(Geometry * inc);

/** \brief add the argument to the Feature list unconditionnally*/
	void forceAdd(Geometry * inc) ;

	void print() const ;
	
	bool isFilled() const ;
	
	Point center() const ;

} ;

	
/** \brief access grid for Features*/
class Grid
{
protected:
	
	double x ;
	double y ;
	Point c ;
	size_t lengthX ;
	size_t lengthY ;
	
	double psize ;
public:
	std::valarray< std::valarray<Pixel *> > pixels;
/** \brief Copnstruct a grid from a size, an initial number of divisions and a center
*
* @param sizeX Length of the access grid
* @param sizeY width of the access grid
* @param div maximum number of spacial divisions to use
* @param center center of the grid
 */
	Grid(double sizeX, double sizeY, int div, const Point & center );
	
	~Grid() ;

/** \brief Add a Feature if it does not overlap with another allready present Feature
*
* @param inc Feature to add
* @return true if insertion was successful
*/
	bool add(Geometry * inc);

/** \brief Add a Feature unconditionnally*/
	void forceAdd(Geometry * inc) ;

/** \brief Return the lis of Features overlapping the argument*/
	std::vector<Geometry *> coOccur(const Geometry * geo) const ;

/** \brief Return the list of Features containing the argument*/
	std::vector<Geometry *> coOccur(const Point & p) const ;

/** \brief Get a new grid, given a number of divisions*/
	Grid getGrid(int div) const;
} ;

/** \brief access grid for Features*/
class Grid3D
{
protected:
	
	std::vector<Voxel *> unfilledpixel ;
	double x ;
	double y ;
	double z ;
	Point c ;
	size_t lengthX ;
	size_t lengthY ;
	size_t lengthZ ;
	
	double psize ;
	int dirtyCounter ;
public:
	std::valarray<std::valarray< std::valarray<Voxel *> > > pixels;
/** \brief Copnstruct a grid from a size, an initial number of divisions and a center
*
* @param sizeX Length of the access grid
* @param sizeY width of the access grid
* @param sizeZ breadth of the access grid
* @param div maximum number of spacial divisions to use
* @param center center of the grid
 */
	Grid3D(double sizeX, double sizeY, double sizeZ, int div, const Point & center );

/** \brief Return a Point contained in no allready placed Feature in the grid*/
	Point randomFreeCenter() const ;
	~Grid3D() ;

/** \brief Add a Feature if it does not overlap with another allready present Feature
*
* @param inc Feature to add
* @return true if insertion was successful
*/
	bool add(Geometry * inc);

/** \brief Add a Feature unconditionnally*/
	void forceAdd(Geometry * inc) ;

/** \brief Return the lis of Features overlapping the argument*/
	std::vector<Geometry *> coOccur(const Geometry * geo) const ;

/** \brief Return the list of Features containing the argument*/
	std::vector<Geometry *> coOccur(const Point & p) const;
	double fraction() const ;

/** \brief Get a new grid, given a number of divisions*/
	Grid3D getGrid(int div) const ;
} ;
} ;















#endif //GRID_H
