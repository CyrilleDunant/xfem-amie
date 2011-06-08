#ifndef PHASE_H
#define PHASE_H

#include "../../features/feature_base.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
#include "../../elements/integrable_entity.h"

namespace Mu
{

/** Simple class containing basic information about behaviours*/
struct Phase
{
private:
    Form * behaviour ;

public:
    Matrix C ;
    Vector beta ;
    double volume ;

    Matrix A ;

public:
    Phase() ;
    Phase(DelaunayTriangle * tri) ;
    Phase(DelaunayTetrahedron * tet) ;
    Phase(Feature * f) ;
    Phase(const Phase & p) ;

    Form * getBehaviour() ;

private:
    void stiffnessFromBehaviour() ;
    void expansionFromBehaviour() ;

};

}

#endif // PHASE_H
