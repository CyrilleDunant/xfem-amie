#ifndef PHASE_H
#define PHASE_H

#include "../../features/feature_base.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../../elements/integrable_entity.h"

namespace Amie
{

typedef enum{

INCLUSION_IS_SPHERE,
INCLUSION_IS_ELLIPSOID,
INCLUSION_IS_CYLINDER
} InclusionGeometryType ;
    
/** Simple class containing basic information about behaviours*/
struct Phase
{

private:
	Form * behaviour ;

public:
	Matrix C ;
	Vector beta ;
	std::vector<Vector> lambda ;
	double volume ;

	Matrix A ;
    
    InclusionGeometryType t;
    double a;
    double b;
    double c;

public:
	Phase(InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	Phase(DelaunayTriangle * tri, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	Phase(DelaunayTetrahedron * tet, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	Phase(Feature * f, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	Phase(Feature * f, DelaunayTriangle * tri, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	Phase(Form * behaviour, double f, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	Phase(const Phase & p, InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;

	virtual void apply( InclusionGeometryType t = INCLUSION_IS_SPHERE, double a=1, double b=1, double c=1) ;
	Form * getBehaviour() ;
	
	void print() ;

	Phase & operator = (const Phase & p) ;

private:
	void stiffnessFromBehaviour(SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;
	void expansionFromBehaviour(SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;
	void ruptureFromBehaviour(SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;

};

}

#endif // PHASE_H
