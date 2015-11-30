// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "../../utilities/parser.h"
#include "../../utilities/itoa.h"
#include "../../utilities/matrixops.h"
#include "../../utilities/enumeration_translator.h"


#include <fstream>
#include <iostream>
#include <cmath>
#include <typeinfo>
#include <limits>
#include <time.h>
#include <sys/time.h>


using namespace Amie ;

struct TrgHeader
{
    size_t vertex ;
    size_t n_trg ;
    size_t n_fields ;
    bool ok ;

    TrgHeader( size_t v = 3, size_t t = 0, size_t f = 2) : vertex(v), n_trg(t), n_fields(f), ok(false) { } ; 

    bool fits(TrgHeader trg) 
    {
        return vertex == trg.vertex && n_trg == trg.n_trg && n_fields == trg.n_fields ;
    }
} ;

TrgHeader getValues( std::string file, std::vector<Vector> & val, TrgHeader trg = TrgHeader() ) 
{
    std::fstream reader ;
    reader.open(file.c_str(), std::ios::in) ;
    TrgHeader ret ;
    if(!reader.is_open())
    {
        std::cout << "cannot open file " << file << ", ignored" << std::endl ;
        return ret ;
    }

    std::cout << "parsing trg file " << file << std::flush ;
    std::string buffer ;
    size_t vertex = 0 ;
    size_t n_trg = 0 ;
    size_t n_fields = 0 ;
    reader >> buffer ; // TRIANGLES
    reader >> n_trg ;
    reader >> vertex ;
    reader >> n_fields ;

    ret.vertex = vertex ;
    ret.n_trg = n_trg ;
    ret.n_fields = n_fields ;

    if( trg.ok )
    {
        if(! ret.fits(trg) )
        {
            std::cout << std::endl ;
            std::cout << "incorrect header for file " << file << ", ignored" << std::endl ;
            return ret ;
        }
    }


    size_t line =  vertex*(n_fields+2) ;
    double tmp = 0 ;
    for(size_t i = 0 ; i < n_trg ; i++)
    {
        if(val.size() < i+1)
            val.push_back(Vector( line )) ;
        for(size_t j = 0 ; j < line ; j++)
        {
            reader >> tmp ;
            val[i][j] = tmp ;
        }
    }
    std::cout << "\rparsing trg file " << file << " ... done" << std::endl ;
    ret.ok = true ;
    return ret ;
}

std::vector< std::pair< size_t, size_t > > makeIndex(  BoundingBoxPosition position, std::vector<Vector> & val, size_t vertex, double tol)
{
    std::vector< std::pair< size_t, size_t > > index ;
    std::cout << "building index cache for position " << Enum::fromBoundingBoxPosition(position) << std::flush ;
    switch(position)
    {
        case TOP:
        {
            double maxy = val[0][1] ;
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                    maxy = std::max( maxy, val[i][2*j+1] ) ;
            }
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                {
                    if(val[i][2*j+1] > maxy-tol)
                        index.push_back(std::make_pair(i,j)) ;
                }
            }
            break ;
        }
        case BOTTOM:
        {
            double miny = val[0][1] ;
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                    miny = std::min( miny, val[i][2*j+1] ) ;
            }
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                {
                    if(val[i][2*j+1] < miny+tol)
                        index.push_back(std::make_pair(i,j)) ;
                }
            }
            break ;
        }
        case RIGHT:
        {
            double maxx = val[0][0] ;
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                    maxx = std::max( maxx, val[i][2*j] ) ;
            }
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                {
                    if(val[i][2*j] > maxx-tol)
                        index.push_back(std::make_pair(i,j)) ;
                }
            }
            break ;
        }
        case LEFT:
        {
            double minx = val[0][0] ;
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                    minx = std::min( minx, val[i][2*j] ) ;
            }
            for(size_t i = 0 ; i < val.size() ; i++)
            {
                for(size_t j = 0 ; j < vertex ; j++)
                {
                    if(val[i][2*j] < minx+tol)
                        index.push_back(std::make_pair(i,j)) ;
                }
            }
            break ;
        }
        default:
            index.clear() ;
    }
    std::cout << "\rbuilding index cache for position " << Enum::fromBoundingBoxPosition(position) << "... done" << std::endl ;
    return index ;
}

double getMinMaxValue( std::vector<Vector> & values, std::vector< std::pair< size_t, size_t > > & index, size_t field, bool min, size_t vertex )
{
    double ret = 0 ;
    for(size_t i = 0 ; i < index.size() ; i++)
    {
        double v = values[ index[i].first] [(2+field)*vertex + index[i].second] ;
        if( i == 0 || (v > ret && !min) ||  (v < ret && min))
            ret = v ;
    }
    return ret ;
}

std::vector<std::string> findFiles( std::string base, size_t max = 100 )
{
    std::vector<std::string> ret ;
    std::fstream stream ;
    if(base.find("%i") == std::string::npos)
    {
        stream.open( base.c_str(), std::ios::in ) ;
        if(stream.is_open())
            ret.push_back(base) ;
        stream.close() ;
        return ret ;
    }

    std::string before = base.substr( 0, base.find("%i") ) ;
    for(size_t i = 0 ; i < max+1 ; i++)
    {
        std::string test = before + itoa(i) ;
        stream.open( test.c_str(), std::ios::in ) ;
        if(stream.is_open())
            ret.push_back(test) ;
        stream.close() ;
    }
    std::cout << ret.size() << " files found matching pattern " << base << std::endl ;    

    return ret ;
}


int main( int argc, char *argv[] )
{
    CommandLineParser parser("Reads an AMIE trg file and obtains maximum values on a certain boundary") ;
    parser.addArgument("file_name","file","relative path to the trg file to read (use %i for multiple files)") ;
    parser.addValue("--tolerance",0.00001,"spatial tolerance on the edge-detection", "-t") ;
    parser.addValue("--steps",100,"maximum number of trg files read", "-s") ;
    parser.addValue("--x",0,"index of the field in the horizontal direction (default 0: displacements along XI)", "-x") ;
    parser.addValue("--y",1,"index of the field in the vertical direction (default 1: displacements along ETA)", "-y") ;
    parser.addString("--out","","file in which the results are stored (results are printed in the console if no file is defined)", "-o") ;
    parser.disableFeatureTreeArguments() ;
    parser.parseCommandLine(argc, argv) ;

    std::string file = parser.getStringArgument("file_name") ;
    std::string out = parser.getString("--out") ;
    double tol = parser.getValue("--tolerance") ;
    size_t x = parser.getValue("--x") ;
    size_t y = parser.getValue("--y") ;
    size_t steps = parser.getValue("--steps") ;

    std::vector<std::string> files = findFiles( file, steps) ;
    if(files.size() == 0)
    {
        std::cout << "no files found, exiting now..." << std::endl ;
        return 0 ;
    }

    Matrix ret( files.size(), 8 ) ;

    std::vector<Vector> values ;
    TrgHeader head = getValues( files[0], values ) ;
    if( x+3 > head.n_fields || y+3 > head.n_fields)
    {
        std::cout << "fields not found, exiting now..." << std::endl ;
        return 0 ;
    }

    std::vector< std::pair< size_t, size_t > > bottom = makeIndex( BOTTOM, values, head.vertex, tol ) ;
    std::vector< std::pair< size_t, size_t > > top = makeIndex( TOP, values, head.vertex, tol ) ;
    std::vector< std::pair< size_t, size_t > > left = makeIndex( LEFT, values, head.vertex, tol ) ;
    std::vector< std::pair< size_t, size_t > > right = makeIndex( RIGHT, values, head.vertex, tol ) ;
    ret[0][0] = getMinMaxValue( values, bottom, y, true, head.vertex )  ;
    ret[0][1] = getMinMaxValue( values, bottom, y, false, head.vertex )  ;
    ret[0][2] = getMinMaxValue( values, top, y, true, head.vertex )  ;
    ret[0][3] = getMinMaxValue( values, top, y, false, head.vertex )  ;
    ret[0][4] = getMinMaxValue( values, left, x, true, head.vertex )  ;
    ret[0][5] = getMinMaxValue( values, left, x, false, head.vertex )  ;
    ret[0][6] = getMinMaxValue( values, right, x, true, head.vertex )  ;
    ret[0][7] = getMinMaxValue( values, right, x, false, head.vertex )  ;

    for(size_t i = 1 ; i < files.size() ; i++)
    {
        TrgHeader next = getValues( files[i], values, head ) ;
        if(next.ok)
        {
            ret[i][0] = getMinMaxValue( values, bottom, y, true, head.vertex )  ;
            ret[i][1] = getMinMaxValue( values, bottom, y, false, head.vertex )  ;
            ret[i][2] = getMinMaxValue( values, top, y, true, head.vertex )  ;
            ret[i][3] = getMinMaxValue( values, top, y, false, head.vertex )  ;
            ret[i][4] = getMinMaxValue( values, left, x, true, head.vertex )  ;
            ret[i][5] = getMinMaxValue( values, left, x, false, head.vertex )  ;
            ret[i][6] = getMinMaxValue( values, right, x, true, head.vertex )  ;
            ret[i][7] = getMinMaxValue( values, right, x, false, head.vertex )  ;
        }
    }

    if(out.size() == 0)
        ret.print() ;
    else
    {
        std::fstream stream ;
        stream.open(out.c_str(), std::ios::out) ;
        for(size_t i = 0 ; i < ret.numRows() ; i++)
        {
            for(size_t j = 0 ; j < ret.numCols() ; j++)
                stream << ret[i][j] << "\t" ;
            stream << std::endl ;
        }
        stream.close() ;
        std::cout << "results exported in " << out << std::endl ;
    }


    return 0 ;
}
