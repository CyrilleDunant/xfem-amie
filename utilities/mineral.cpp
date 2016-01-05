#include "mineral.h"
#include "enumeration_translator.h"
#include <fstream>

using namespace Amie ;

// cij must contain the elastic tensor components as "c11" "c12" etc
Mineral::Mineral(SymmetryType sym, std::map<std::string, double> cij ) 
{
   symmetry = sym ;
   valid = false ;
   for(auto c = cij.begin() ; c != cij.end() ; c++)
       add( c->first, c->second ) ;
   stiffness = getStiffnessMatrix() ;
   if(stiffness.array().max() > 0)
       valid = true ;
}


// reads a mineral from a file, assuming that the file contains information such as
// mineral.symmetry = "triclinic" // or whatever
// mineral.c11 = [ 15 17 16 ] // values from different references are stored in a vector
Mineral::Mineral( std::string path, std::string sep, int index, double factor, bool force )
{
    symmetry = SYMMETRY_CUBIC ;
    valid = false ;

    std::string filename = path.substr( path.rfind("/")+1 ) ;
    name = filename.substr(0,filename.find_first_of(sep)) ;
    
    std::fstream reader ;
    reader.open( path.c_str(), std::ios::in ) ;
    if(!reader.fail())
    {
        for( std::string line ; getline( reader, line ) ; )
        {
            if(line.find( name ) == 0)
            {
                std::string key = line.substr( line.find(".")+1, line.find(" ")-line.find(".")-1) ;
                if(line.find("\"") != std::string::npos)
                {
                    size_t start = line.find("\"") ;
                    size_t end = line.rfind("\"") ;
                    std::string value = line.substr( start+1, end-start-1 ) ;
                    information[key] = value ;
                }
                else if(line.find("[") != std::string::npos)
                {
                    size_t start = line.find("[") ;
                    size_t end = line.find("]", start+1) ;
                    std::string all = line.substr( start+1, end-start-1 ) ;
                    while(all.size() > 0 && all[0] == ' ') { all = all.substr(1) ; }
                    std::vector<double> data ;
                    while(all.size() > 0)
                    {
                        size_t sep = all.find(" ") ;
                        std::string token = all ;
                        if(sep != std::string::npos)
                        {
                            token = all.substr( 0, sep ) ;
                            all = all.substr(sep+1) ;
                            while(all.size() > 0 && all[0] == ' ') { all = all.substr(1) ; }
                        }
                        if(sep == std::string::npos)
                        {
                            all = std::string() ;
                        }
                        data.push_back( atof(token.c_str() ) ) ;
                    }
                    set( key, data ) ;
                }
            }
        }
        if( information.find("symmetry") != information.end() )
        {
            bool ok = false ;
            std::string sym = "SYMMETRY_"+information["symmetry"] ;
            for(size_t i = 8 ; i < sym.size() ; i++)
                sym[i] = toupper( sym[i] ) ;
            symmetry = Enum::getSymmetryType( sym, &ok ) ;
            if(ok)
            {
                stiffness = getStiffnessMatrix( index, force ) * factor ;
                if(stiffness.array().max() > 0)
                    valid = true ;
            }
        }
    }
}

void Mineral::set(std::string key, std::vector<double> v ) 
{
    Vector val(v.size()) ;
    for(size_t i = 0 ; i < v.size() ; i++)
       val[i] = v[i] ;
    set( key, val ) ;
}

void Mineral::add(std::string key, double v ) 
{
    Vector val ;
    if(components.find(key) == components.end())
    {
        val.resize( 1 ) ;
    }
    else
    {
        val.resize( components[key].size()+1 ) ;
        for(size_t i = 0 ; i < components[key].size() ; i++)
            val[i] = components[key][i] ;
    }
    val[val.size()-1] = v ;
    set( key, val ) ;
}



double Mineral::getElementaryComponent( std::string key, int index ) 
{
    if( components.find( key) == components.end() ) { return 0 ; }
    if(index < 0 || index >= (int) components[key].size() ) { return components[key].sum()/components[key].size() ; }
    return components[key][index] ;
}

Vector Mineral::getElementaryComponents( int index )
{
   Vector c ;
   switch(symmetry)
   {
        case SYMMETRY_CUBIC:
        {
             c.resize(3) ;
             c[0] = getElementaryComponent("c11", index) ;
             c[1] = getElementaryComponent("c44", index) ;
             c[2] = getElementaryComponent("c12", index) ;
             break ;
        }
        case SYMMETRY_HEXAGONAL:
        {
             c.resize(5) ;
             c[0] = getElementaryComponent("c11", index) ;
             c[1] = getElementaryComponent("c33", index) ;
             c[2] = getElementaryComponent("c44", index) ;
             c[3] = getElementaryComponent("c12", index) ;
             c[4] = getElementaryComponent("c13", index) ;
             break ;
        }
        case SYMMETRY_MONOCLINIC:
        {
             c.resize(13) ;
             c[0] = getElementaryComponent("c11", index) ;
             c[1] = getElementaryComponent("c22", index) ;
             c[2] = getElementaryComponent("c33", index) ;
             c[3] = getElementaryComponent("c44", index) ;
             c[4] = getElementaryComponent("c55", index) ;
             c[5] = getElementaryComponent("c66", index) ;
             c[6] = getElementaryComponent("c12", index) ;
             c[7] = getElementaryComponent("c13", index) ;
             c[8] = getElementaryComponent("c16", index) ;
             c[9] = getElementaryComponent("c23", index) ;
             c[10] = getElementaryComponent("c26", index) ;
             c[11] = getElementaryComponent("c36", index) ;
             c[12] = getElementaryComponent("c45", index) ;
             break ;
        }
        case SYMMETRY_ORTHORHOMBIC:
        {
             c.resize(9) ;
             c[0] = getElementaryComponent("c11", index) ;
             c[1] = getElementaryComponent("c22", index) ;
             c[2] = getElementaryComponent("c33", index) ;
             c[3] = getElementaryComponent("c44", index) ;
             c[4] = getElementaryComponent("c55", index) ;
             c[5] = getElementaryComponent("c66", index) ;
             c[6] = getElementaryComponent("c12", index) ;
             c[7] = getElementaryComponent("c13", index) ;
             c[8] = getElementaryComponent("c23", index) ;
             break ;
        }
        case SYMMETRY_TETRAGONAL:
        {
             c.resize(7) ;
             c[0] = getElementaryComponent("c11", index) ;
             c[1] = getElementaryComponent("c33", index) ;
             c[2] = getElementaryComponent("c44", index) ;
             c[3] = getElementaryComponent("c66", index) ;
             c[4] = getElementaryComponent("c12", index) ;
             c[5] = getElementaryComponent("c13", index) ;
             c[6] = getElementaryComponent("c16", index) ;
             break ;
        }
        case SYMMETRY_TRIGONAL:
        {
             c.resize(7) ;
             c[0] = getElementaryComponent("c11", index) ;
             c[1] = getElementaryComponent("c33", index) ;
             c[2] = getElementaryComponent("c44", index) ;
             c[3] = getElementaryComponent("c66", index) ;
             c[4] = getElementaryComponent("c12", index) ;
             c[5] = getElementaryComponent("c13", index) ;
             c[6] = getElementaryComponent("c14", index) ;
             break ;
        }
        case SYMMETRY_TRICLINIC:
        {
             c.resize(21) ;
             c[0] = getElementaryComponent("c11", index) ;
             c[1] = getElementaryComponent("c22", index) ;
             c[2] = getElementaryComponent("c33", index) ;
             c[3] = getElementaryComponent("c44", index) ;
             c[4] = getElementaryComponent("c55", index) ;
             c[5] = getElementaryComponent("c66", index) ;
             c[6] = getElementaryComponent("c12", index) ;
             c[7] = getElementaryComponent("c13", index) ;
             c[8] = getElementaryComponent("c14", index) ;
             c[9] = getElementaryComponent("c15", index) ;
             c[10] = getElementaryComponent("c16", index) ;
             c[11] = getElementaryComponent("c23", index) ;
             c[12] = getElementaryComponent("c24", index) ;
             c[13] = getElementaryComponent("c25", index) ;
             c[14] = getElementaryComponent("c26", index) ;
             c[15] = getElementaryComponent("c34", index) ;
             c[16] = getElementaryComponent("c35", index) ;
             c[17] = getElementaryComponent("c36", index) ;
             c[18] = getElementaryComponent("c45", index) ;
             c[19] = getElementaryComponent("c46", index) ;
             c[20] = getElementaryComponent("c56", index) ;
             break ;
        }
   }
   return c ;
}

Matrix Mineral::getStiffnessMatrix( int index, bool force ) 
{
   return Tensor::orthotropicCauchyGreen( getElementaryComponents(index), symmetry, force ) ;
}

