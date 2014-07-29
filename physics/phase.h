#ifndef PHASE_H
#define PHASE_H

namespace Amie
{

struct Phase
{
public:
    Form * behaviour ;
    double volume ;

public:
    Phase(DelaunayTriangle * tri) ;
    Phase(DelaunayTetrahedron * tet) ;
    Phase(Feature * f) ;
};

}


#endif // PHASE_H
