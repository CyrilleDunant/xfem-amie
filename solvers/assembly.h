
//
// C++ Interface: assembly
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __ASSEMBLY_NEW_H_
#define __ASSEMBLY_NEW_H_

#include <map>
#include <set>
#include <complex>
#include <iomanip>

#include "../delaunay.h"
#include "../delaunay_3d.h"
#include "../elements/elements.h"
#include "../physics/physics.h"
#include "../sparse/sparse_matrix.h"

namespace Mu
{


// class Assembly
// {
// protected:
// 	std::map< std::pair<size_t, size_t>, double > morseMatrix ;
// 	std::vector<std::vector<double > > mother ;
// 	bool cached ;
// 	Form *l ;
// 	
// public:
// 	Assembly( Form *w) ;
// 	
// 	virtual ~Assembly() ;
// 	
// 	void operator +=(const ElementarySurface * e) ;
// 	
// 	Matrix finalise() const ;
// 	
// 	void print() const ;
// 	
// } ;

typedef enum 
{
	GENERAL = 0,
	FIX_ALONG_XI,
	SET_ALONG_XI,
	FIX_ALONG_ETA,
	SET_ALONG_ETA,
	FIX_ALONG_ZETA,
	SET_ALONG_ZETA,
	FIX_ALONG_XI_ETA,
	SET_ALONG_XI_ETA,
	FIX_ALONG_XI_ZETA,
	SET_ALONG_XI_ZETA,
	FIX_ALONG_ETA_ZETA,
	SET_ALONG_ETA_ZETA,
	SET_FORCE_XI = 16,
	SET_FORCE_ETA = 17,
	SET_FORCE_ZETA = 18 ,
} LagrangeMultiplierType ;

class LagrangeMultiplier
{
public:
	std::valarray<unsigned int> ids ;
	Vector coefs ;
	int id ;
	double value ;
	std::vector<std::pair<unsigned int, double> > hints ;
public:
	LagrangeMultiplier(std::valarray<unsigned int> i, Vector c, double v ,int my_id  = -1) ;
	LagrangeMultiplier(int _id){ id  = id ;} 
	
	CoordinateIndexedIncompleteSparseMatrix getMatrix() const ;
	void setId(const int) ;
	int getId() const ;
	
	const std::valarray<unsigned int> getDofIds() const ;
	double getValue() const  ;
	std::vector<std::pair<unsigned int, double> > getHints() const ;
	void setHint(std::pair<unsigned int, double>) ;
	
	bool operator < (const LagrangeMultiplier & m) const
	{
		return id < m.getId() ;
	}
	
	bool operator == (const LagrangeMultiplier & m) const
	{
		return id == m.getId() ;
	}
	
	void print() const
	{
		std::cout << "Multiplier affectig ID " << id << ", val = " << value << std::endl ;
	}
	
	LagrangeMultiplierType type ;
} ;

struct MultiplierHasIdSupEq
{
	int m_id ;
	MultiplierHasIdSupEq(int id) {m_id = id ;}
	bool operator()(const LagrangeMultiplier & m)
	{
		return m_id >= m.getId() ;
	}
} ;

struct MultiplierHasIdInfEq
{
	int m_id ;
	MultiplierHasIdInfEq(int id) {m_id = id ;}
	bool operator()(const LagrangeMultiplier & m)
	{
		return m_id <= m.getId() ;
	}
} ;

struct MultiplierHasId
{
	int m_id ;
	MultiplierHasId(int id) {m_id = id ;}
	bool operator()(const LagrangeMultiplier & m)
	{
		return m_id == m.getId() ;
	}
} ;

class Assembly
{
protected:
		
	std::vector<ElementarySurface *> element2d ;
	std::vector<ElementaryVolume *> element3d ;
	std::vector<LagrangeMultiplier> multipliers ;
	bool has3Dims ; 

	Vector externalForces ;
	Vector displacements ;
	Vector nonLinearExternalForces ;
	CoordinateIndexedSparseMatrix * coordinateIndexedMatrix ;
	CoordinateIndexedIncompleteSparseMatrix * nonLinearPartialMatrix ;
	std::map<std::pair<size_t, size_t>, double > * boundaryMatrix ;
	
	void make_final() ;
	
	
	size_t multiplier_offset ;

public:
	Assembly() ;
	
	virtual ~Assembly() ;
	
	void operator +=(ElementarySurface * e) ;
	void operator +=(ElementaryVolume * e) ;
	void add(ElementarySurface * e) ;
	void add(ElementaryVolume * e) ;
	void update(ElementarySurface * e, double dt, Vector * params) ;
	void update(ElementaryVolume * e, double dt, Vector * params) ;
	bool is2D() const;
	bool is3D() const;	
	void set3D() ;
	void set2D() ;
	void print() ;
	void printDiag() const ;
	
	void setBoundaryConditions() ;
	
	bool nonLinearStep() ;
// 	Vector solve(size_t maxit = 1000) ;
	Vector & solve(Vector x, size_t maxit = 1000, const bool verbose = false) ;
	Vector & cgsolve(Vector x0 = Vector(0), size_t maxit = 10000) ;
	Vector cgnpsolve(Vector b, size_t maxit) ;
	
	CoordinateIndexedSparseMatrix & getMatrix() ;
	const CoordinateIndexedSparseMatrix & getMatrix() const;
	CoordinateIndexedIncompleteSparseMatrix & getNonLinearMatrix() ;
	Vector & getForces() ;
	Vector & getNonLinearForces() ;
	Vector & getDisplacements() {return displacements ;}
// 	const Vector & getForces() const;
	
	std::vector<ElementarySurface *> getElements2d() const ;
	std::vector<ElementarySurface *> & getElements2d()  ;
	ElementarySurface * getElement2d(const size_t ) const ;
	ElementarySurface * getElement2d(const size_t )  ;
	
	std::vector<ElementaryVolume *> getElements3d() const ;
	std::vector<ElementaryVolume *> & getElements3d()  ;
	ElementaryVolume * getElement3d(const size_t ) const ;
	ElementaryVolume * getElement3d(const size_t )  ;
	
	void fixPoint(size_t id) ;
	void setPoint(double ex, size_t id) ;
	void setPoint(double ex, double ey, size_t id) ;
	void setPoint(double ex, double ey, double ez, size_t id) ;
/** Set unknown of a given node in a given direction to a given value
	 * @param var Variable (XI, ETA, or ZETA) along which displacement is fixed
	 * @param id ID of the node 
	 * @param val value at which unknown is set
 */
	void setPointAlong(Variable, double val, size_t id) ;
	
	/** Apply a force on given node
	 * @param var Variable (XI, ETA, or ZETA) along which force is applied 
	 * @param id ID of the node 
	 * @param val intensity of the force
	 */
	void setForceOn(Variable var, double val, size_t id) ;
	void setForceOn(double val, size_t id) ;
	void fixPoint(size_t id, Mu::Variable v) ;
	void setPeriodicPoint(size_t id0, size_t id1) ;
		
	double froebeniusNorm() ;
	
	void fix() ;
	void clear() ;
	
	
} ;


} ;

std::map<std::pair<size_t, size_t>, Mu::Matrix> incompleteCholeskyDecomposition(std::map<std::pair<size_t, size_t>, Mu::Matrix>  & morseMatrix) ;
Vector operator *(const std::map< std::pair<size_t, size_t>, Mu::Matrix > A , const Vector x) ;
Vector operator *(const std::map< std::pair<size_t, size_t>, double > A , const Vector x) ;
std::map< std::pair<size_t, size_t>, Mu::Matrix > operator *(const std::map< std::pair<size_t, size_t>, Mu::Matrix > A , const std::map< std::pair<size_t, size_t>, Mu::Matrix > B) ;

#endif

 
