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
    std::map<std::string, Vector> components ; // stiffness components according to different sources

    Matrix stiffness ;

public:
    Mineral(SymmetryType sym, std::map<std::string, double> cij = std::map<std::string, double>() ) ;
    Mineral(std::string file, std::string sep = ".-", int index = -1, double factor = 1., bool force = false) ;

    void set(std::string key, Vector val ) { components[key] = val ; }
    void set(std::string key, std::vector<double> v ) ;
    void set(std::string key, double v ) ;

    bool isValid() const { return valid ; }
    bool check() ;

    double getElementaryComponent( std::string key, int index = -1 ) ;
    Vector getElementaryComponents( int index = -1 ) ;    
    Matrix getStiffnessMatrix( int index = -1, bool force = false ) ;
} ;



} 
#endif // MINERAL_H
