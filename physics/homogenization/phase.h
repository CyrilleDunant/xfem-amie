#ifndef PHASE_H
#define PHASE_H

#include "../../features/feature_base.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../../elements/integrable_entity.h"

namespace Amie
{

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

public:
	Phase() ;
	Phase(DelaunayTriangle * tri) ;
	Phase(DelaunayTetrahedron * tet) ;
	Phase(Feature * f) ;
	Phase(Feature * f, DelaunayTriangle * tri) ;
	Phase(Form * behaviour, double f, SpaceDimensionality dim = SPACE_THREE_DIMENSIONAL) ;
	Phase(const Phase & p) ;

	virtual void apply() ;
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
