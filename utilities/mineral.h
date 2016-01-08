#ifndef MINERAL_H
#define MINERAL_H

#include "tensor.h"

namespace Amie
{

struct Mineral
{
protected:
    bool valid ;

public:
    std::string name ;
    SymmetryType symmetry ;

    std::map<std::string, std::string> information ;
    std::map<std::string, std::vector<std::string> > lists ;
    std::map<std::string, std::vector<double> > components ; // stiffness components according to different sources

    Matrix stiffness ;

public:
    Mineral(SymmetryType sym = SYMMETRY_CUBIC, std::map<std::string, double> cij = std::map<std::string, double>() ) ;
    Mineral(std::string file, std::string sep = ".-", int index = -1, double factor = 1., bool force = false) ;

    void set(std::string key, double v ) ;
    void reset(std::string key, double v ) ;
    void multiply(std::string key, double v, int index = -1 ) ;

    bool isValid() const { return valid ; }
    bool check() ;

    double get( std::string key, int index = -1 ) ;
    Vector getElementaryComponents( int index = -1 ) ;    
    Matrix getStiffnessMatrix( int index = -1, bool force = false ) ;
} ;



} 
#endif // MINERAL_H
