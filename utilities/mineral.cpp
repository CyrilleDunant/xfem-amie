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
       set( c->first, c->second ) ;
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
                std::string key = line.substr( line.find(".")+1, line.find_first_of(" =")-line.find(".")-1) ;
                if(line.find("\"") != std::string::npos)
                {
                    size_t start = line.find("\"") ;
                    size_t end = line.rfind("\"") ;
                    std::string value = line.substr( start+1, end-start-1 ) ;
                    information[key] = value ;
                }
                else if(line.find("[") != std::string::npos)
                {
                    if( key[0] == 'c' && key.size() == 3 ) { key[0] = 'C' ; }
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
                    components[key] = data ;
                }
                else if( line.find("list(") != std::string::npos)
                {
                    size_t start = line.find("(") ;
                    size_t end = line.find(")") ;
                    std::vector<std::string> list ;
                    std::string current ;
                    for(size_t i = start+1 ; i < end ; i++)
                    {
                        if( line[i] == ' ' ) { continue ; }
                        if( line[i] == ',' ) 
                        {
                            list.push_back( current ) ;
                            current = std::string() ;
                        }
                        else
                            current.append(1, line[i] ) ;
                    }
                    list.push_back(current) ;
                    lists[key] = list ;
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
    else
    {
        std::cout << "file " << path << " not found!" << std::endl ;
    }
}

bool Mineral::check()
{
    valid = false ;
    Matrix tmp = getStiffnessMatrix() ;
    if(tmp.array().max() > 0)
        valid = true ;
    return valid ;
}

void Mineral::set(std::string key, double v ) 
{
    if(components.find(key) == components.end())
    {
        std::vector<double> val ;
        val.push_back(v) ;
        components[key] = val ;
    }
    else
    {
        components[key].push_back(v) ;
    }
}

void Mineral::reset(std::string key, double v ) 
{
    std::vector<double> val ;
    val.push_back(v) ;
    components[key] = val ;
}

void Mineral::multiply(std::string key, double v, int index ) 
{
    if( components.find(key) == components.end() ) { return ; }
    if(index < 0 || (size_t) index >= components[key].size() )
    {
        for(size_t i = 0 ; i < components[key].size() ; i++)
            components[key][i] *= v ;
    }
    else
        components[key][index] *= v ;
}

double Mineral::get( std::string key, int index ) 
{
    if( components.find( key) == components.end() ) { return 0 ; }
    if(index < 0 || index >= (int) components[key].size() ) 
    {
        if(components[key].size() == 0) { return 0 ; }
        double sum = 0 ;
        for(size_t i = 0 ; i < components[key].size() ; i++)
            sum += components[key][i] ;
        return sum/components[key].size() ; 
    }
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
             c[0] = get("C11", index) ;
             c[1] = get("C44", index) ;
             c[2] = get("C12", index) ;
             break ;
        }
        case SYMMETRY_HEXAGONAL:
        {
             c.resize(5) ;
             c[0] = get("C11", index) ;
             c[1] = get("C33", index) ;
             c[2] = get("C44", index) ;
             c[3] = get("C12", index) ;
             c[4] = get("C13", index) ;
             break ;
        }
        case SYMMETRY_MONOCLINIC:
        {
             c.resize(13) ;
             c[0] = get("C11", index) ;
             c[1] = get("C22", index) ;
             c[2] = get("C33", index) ;
             c[3] = get("C44", index) ;
             c[4] = get("C55", index) ;
             c[5] = get("C66", index) ;
             c[6] = get("C12", index) ;
             c[7] = get("C13", index) ;
             c[8] = get("C16", index) ;
             c[9] = get("C23", index) ;
             c[10] = get("C26", index) ;
             c[11] = get("C36", index) ;
             c[12] = get("C45", index) ;
             break ;
        }
        case SYMMETRY_ORTHORHOMBIC:
        {
             c.resize(9) ;
             c[0] = get("C11", index) ;
             c[1] = get("C22", index) ;
             c[2] = get("C33", index) ;
             c[3] = get("C44", index) ;
             c[4] = get("C55", index) ;
             c[5] = get("C66", index) ;
             c[6] = get("C12", index) ;
             c[7] = get("C13", index) ;
             c[8] = get("C23", index) ;
             break ;
        }
        case SYMMETRY_TETRAGONAL:
        {
             c.resize(7) ;
             c[0] = get("C11", index) ;
             c[1] = get("C33", index) ;
             c[2] = get("C44", index) ;
             c[3] = get("C66", index) ;
             c[4] = get("C12", index) ;
             c[5] = get("C13", index) ;
             c[6] = get("C16", index) ;
             break ;
        }
        case SYMMETRY_TRIGONAL:
        {
             c.resize(7) ;
             c[0] = get("C11", index) ;
             c[1] = get("C33", index) ;
             c[2] = get("C44", index) ;
             c[3] = get("C66", index) ;
             c[4] = get("C12", index) ;
             c[5] = get("C13", index) ;
             c[6] = get("C14", index) ;
             break ;
        }
        case SYMMETRY_TRICLINIC:
        {
             c.resize(21) ;
             c[0] = get("C11", index) ;
             c[1] = get("C22", index) ;
             c[2] = get("C33", index) ;
             c[3] = get("C44", index) ;
             c[4] = get("C55", index) ;
             c[5] = get("C66", index) ;
             c[6] = get("C12", index) ;
             c[7] = get("C13", index) ;
             c[8] = get("C14", index) ;
             c[9] = get("C15", index) ;
             c[10] = get("C16", index) ;
             c[11] = get("C23", index) ;
             c[12] = get("C24", index) ;
             c[13] = get("C25", index) ;
             c[14] = get("C26", index) ;
             c[15] = get("C34", index) ;
             c[16] = get("C35", index) ;
             c[17] = get("C36", index) ;
             c[18] = get("C45", index) ;
             c[19] = get("C46", index) ;
             c[20] = get("C56", index) ;
             break ;
        }
   }
   return c ;
}

Matrix Mineral::getStiffnessMatrix( int index, bool force ) 
{
   return Tensor::orthotropicCauchyGreen( getElementaryComponents(index), symmetry, force ) ;
}

