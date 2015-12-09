//
// C++ Interface: voxel writer
//
// Description:
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "exodus_writer.h"
#include "../enumeration_translator.h"
#include "../../features/features.h"
#include "../../elements/integrable_entity.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>


namespace Amie
{

void ExodusTriangleWriter::setCache(std::vector<Feature *> feats)
{
    if(feats.size() > 0)
    {
        std::vector<Geometry *> geom ;
        for(size_t i = 0 ; i < feats.size() ; i++)
            geom.push_back( dynamic_cast<Geometry *>(feats[i]) ) ;
        index.push_back( source->get2DMesh()->generateCache(geom) ) ;
    }
}

void ExodusTriangleWriter::setCache(std::vector< std::vector<Feature *> > feats, bool last)
{
    for(size_t i = 0 ; i < feats.size() ; i++)
        setCache(feats[i]) ;
    if(last)
        finalizeCache() ;
}

void ExodusTriangleWriter::finalizeCache()
{
    unsigned int first = source->get2DMesh()->generateCacheOut( index ) ;
    index.insert( index.begin(), first ) ;

    sides.clear() ;
    Feature * box= source->getFeature(0) ;
    sides.push_back( std::make_pair( new Segment( box->getBoundingBox()[0], box->getBoundingBox()[1] ), 0 ) ) ;
    sides.push_back( std::make_pair( new Segment( box->getBoundingBox()[1], box->getBoundingBox()[2] ), 1 ) ) ;
    sides.push_back( std::make_pair( new Segment( box->getBoundingBox()[2], box->getBoundingBox()[3] ), 2 ) ) ;
    sides.push_back( std::make_pair( new Segment( box->getBoundingBox()[3], box->getBoundingBox()[0] ), 3 ) ) ;
    for(size_t i = 0 ; i < 4 ; i++)
        sides[i].second = source->get2DMesh()->generateCache( sides[i].first, POINT_TOLERANCE, true ) ; 


    elems.clear() ;
    unsigned int counter = 1 ;
    for(size_t i = 0 ; i < index.size() ; i++)
    {
        std::vector<int> cache = source->get2DMesh()->getCache( index[i] ) ;
        for(size_t j = 0 ; j < cache.size() ; j++)
        {
            elems[ cache[j] ] = counter ;
            counter++ ;
        }
    }
}

void ExodusTriangleWriter::write() 
{
    if(index.size() == 0)
        finalizeCache() ;

    std::ofstream out ;
    out.open(filename.c_str(), std::ios::out) ;

    Order order = source->getOrder() ;

    size_t planes = 1 ;
    if(order == LINEAR_TIME_LINEAR || order == QUADRATIC_TIME_LINEAR) { planes = 2 ; }
    if(order == LINEAR_TIME_QUADRATIC || order == QUADRATIC_TIME_QUADRATIC) { planes = 3 ; }

    size_t n = 3 ;
    if(order == QUADRATIC || order == QUADRATIC_TIME_LINEAR || order == QUADRATIC_TIME_QUADRATIC) { n = 6 ;}
    if(order == CUBIC) { n = 9 ; }
    if(order == QUADRIC) { n = 12 ; }

    size_t ps = n/3 ;

    std::vector<Point *> nodes = source->getNodes() ;
    Mesh<DelaunayTriangle, DelaunayTreeItem> * mesh = source->get2DMesh() ;

    out << "netcdf amie {" << std::endl ;
    out << "dimensions:" << std::endl ;
    out << "\tlen_string = 33 ;" << std::endl ;
    out << "\tlen_line = 81 ;" << std::endl ;
    out << "\tfour = 4 ;" << std::endl ;
    out << "\ttime_step = UNLIMITED ; // (0 currently)" << std::endl ;
    out << "\tnum_qa_rec = 1 ;" << std::endl ;
    out << "\tnum_dim = 2 ;" << std::endl ;
    out << "\tnum_nodes = " << nodes.size()/planes << " ;" << std::endl ;
    out << "\tnum_elem = " << elems.size() << " ;" << std::endl ;
    out << "\tnum_el_blk = " << index.size() << " ;" << std::endl ;
    for(size_t i = 0 ; i < index.size() ; i++)
    {
        out << "\tnum_el_in_blk" << itoa(i+1) << " = " << mesh->getCache(index[i]).size() << " ;" << std::endl ;
        out << "\tnum_nod_per_el" << itoa(i+1) << " = " << n << " ;" << std::endl ;
    }
    out << "\tnum_side_sets = " << sides.size() << " ;" << std::endl ;
    for(size_t i = 0 ; i < sides.size() ; i++)
    {
        out << "\tnum_side_ss" << itoa(i+1) << " = " << mesh->getCache(sides[i].second).size() << " ;" << std::endl ;
        out << "\tnum_df_ss" << itoa(i+1) << " = " << mesh->getCache(sides[i].second).size()*4 << " ;" << std::endl ;
    }
    
    out << std::endl ;

    out << "variables:" << std::endl ;
    out << "\tchar qa_records(num_qa_rec, four, len_string) ;" << std::endl ;
    out << "\tchar coor_names(num_dim, len_string) ;" << std::endl ;
    out << "\tint elem_map(num_elem) ;" << std::endl ;
    out << "\tint eb_prop1(num_el_blk) ;" << std::endl ;
    out << "\t\teb_prop1:name = \"ID\" ;" << std::endl ;
    out << "\tint eb_status(num_el_blk) ;" << std::endl ;
    for(size_t i = 0 ; i < index.size() ; i++)
    {
        std::string j = itoa(i+1) ;
        out << "\tint connect" + j+ "(num_el_in_blk" + j + ", num_nod_per_el" + j + ") ;" << std::endl ;
        out << "\t\tconnect" + j + ":elem_type = \"TRI\" ;" << std::endl ;
    }
    out << "\tint ss_status(num_side_sets) ;" << std::endl ;
    out << "\tint ss_prop1(num_side_sets) ;" << std::endl ;
    out << "\t\tss_prop1:name = \"ID\" ;" << std::endl ;
    out << "\tchar ss_names(num_side_sets, len_string) ;" << std::endl ;
    for(size_t i = 0 ; i < sides.size() ; i++)
    {
        std::string j = itoa(i+1) ;
        out << "\tint elem_ss"+j+"(num_side_ss"+j+") ;" << std::endl ;
        out << "\tint side_ss"+j+"(num_side_ss"+j+") ;" << std::endl ;
        out << "\tdouble dist_fact_ss"+j+"(num_df_ss"+j+") ;" << std::endl ;
    }
    out << "\tdouble coordx(num_nodes) ;" << std::endl ;
    out << "\tdouble coordy(num_nodes) ;" << std::endl ;
    out << std::endl ;

    out << "// global attributes:" << std::endl ;
    out << "  		:api_version = 4.98f ;" << std::endl ;
    out << "		:version = 4.98f ;" << std::endl ;
    out << "		:floating_point_word_size = 8 ;" << std::endl ;
    out << "		:file_size = 1 ;" << std::endl ;
    out << "		:title = \"amie\" ;" << std::endl ;
    out << std::endl ;

    out << "data:" << std::endl ;
    out << std::endl ;

    std::time_t time = std::time(NULL) ;
    std::tm* t = std::localtime( &time ) ;
    
    out << " qa_records = " << std::endl ;
    out << "   \"AMIE\" , " << std::endl ;
    out << "   \"development version\" , " << std::endl ;
    if(us == true)
        out << "   \"" << t->tm_mon+1 << "/" << t->tm_mday << "/" << t->tm_year+1900 << "\", " << std::endl ;
    else
        out << "   \"" << t->tm_mday << "/" << t->tm_mon+1 << "/" << t->tm_year+1900 << "\", " << std::endl ;
    out << "   \"" << t->tm_hour << ":" << t->tm_min << ":" << t->tm_sec << "\" ; " << std::endl ;
    out << std::endl ;

    out << " coor_names = "<< std::endl ;
    out << "  \"x\"," << std::endl ;
    out << "  \"y\" ;" << std::endl ;
    out << std::endl ;

    out << " elem_map = " << std::endl ;
    out << "  " ;
    for(size_t i = 0 ; i < elems.size() ; i++)
    {
        if(i%10 == 0 && i > 0)
        {
            out << " ," << std::endl ;
            out << "  " ;
        }
        else if(i>0)
            out << " , " ;
        out << i+1 ;
    }
    out << " ;" << std::endl ;
    out << std::endl ;

    out << " eb_prop1 = " ;
    for(size_t i = 0 ; i < index.size() ; i++)
    {
        out << i+1 ;
        if(i == index.size() -1)
            out << " ; " << std::endl ;
        else
            out << " , " ;
    }
    out << std::endl ;

    out << " eb_status = " ;
    for(size_t i = 0 ; i < index.size() ; i++)
    {
        out << 1 ;
        if(i == index.size() -1)
            out << " ; " << std::endl ;
        else
            out << " , " ;
    }
    out << std::endl ;

    out << " ss_prop1 = " ;
    for(size_t i = 0 ; i < sides.size() ; i++)
    {
        out << i+1 ;
        if(i == sides.size() -1)
            out << " ; " << std::endl ;
        else
            out << " , " ;
    }
    out << std::endl ;

    out << " ss_names = \"top\", \"right\", \"bottom\", \"left\" ;" << std::endl ;
    out << std::endl ;

    out << " ss_status = 1, 1, 1, 1 ;" << std::endl ;
    out << std::endl ;

    for(size_t i = 0 ; i < sides.size() ; i++)
    {
        std::vector<int> cache = mesh->getCache( sides[i].second ) ;

        out << " elem_ss" << i+1 << " = " ;
        for(size_t j = 0 ; j < cache.size() ; j++)
        {
            out << elems[ cache[j] ] ;
            if(j == cache.size()-1)
                out << " ;" << std::endl ;
            else
                out << ", "  ;
            if( j%10 == 0 && j > 0 && j != cache.size()-1)
            {
                out << std::endl ;
                out << "   " ;
            }
        }
        out << std::endl ;
        out << " side_ss" << i+1 << " = " ;
        for(size_t j = 0 ; j < cache.size() ; j++)
        {
            int s = 0 ;
            DelaunayTriangle * trg = dynamic_cast<DelaunayTriangle *>(mesh->getInTree( cache[j] )) ;
            if( !(sides[i].first->on( trg->getBoundingPoint( 0 )) ) ) { s = 2 ; }
            if( !(sides[i].first->on( trg->getBoundingPoint( ps )) ) ) { s = 3 ; }
            if( !(sides[i].first->on( trg->getBoundingPoint( ps*2 )) ) ) { s = 1 ; }
            out << s ;
            if(j == cache.size()-1)
                out << " ;" << std::endl ;
            else
                out << ", "  ;
            if( j%10 == 0 && j > 0 && j != cache.size()-1)
            {
                out << std::endl ;
                out << "   " ;
            }
        }
        out << std::endl ;
        out << " dist_fact_ss" << i+1 << " = " ;
        for(size_t j = 0 ; j < cache.size() ; j++)
        {
            out << "1, 1, 1, 1" ;
            if(j == cache.size()-1)
                out << " ;" << std::endl ;
            else
                out << ", "  ;
            if( j%3 == 0 && j > 0 && j != cache.size()-1)
            {
                out << std::endl ;
                out << "   " ;
            }
        }
        out << std::endl ;
    }

    for(size_t i = 0 ; i < index.size() ; i++)
    {
        std::string s = itoa(i+1) ;
        out << " connect" + s + " =" << std::endl ;
        std::vector<int> cache = mesh->getCache( index[i] ) ;
        for(size_t j = 0 ; j < cache.size() ; j++)
        {
            DelaunayTriangle * trg = dynamic_cast<DelaunayTriangle *>(mesh->getInTree( cache[j] )) ;
            out << "  " ;
            for(size_t p = 0 ; p < n ; p++)
            {
               out << trg->getBoundingPoint(p).getId()+1 ;
               if(p != n-1)
                   out << ", " ;
            }
            if(j == cache.size()-1)
                out << " ;" << std::endl ;
            else
                out << "," << std::endl ;
        }
        out << std::endl ;
    }

    out << " coordx = " ;
    out << "  " ;
    for(size_t i = 0 ; i < nodes.size()/planes ; i++)
    {
         if(i%10 == 0 && i > 0)
         {
             out << std::endl ;
             out << "  " ;
         }
         out << nodes[i]->getX() ;
         if(i == (nodes.size()/planes)-1)
             out << " ;" << std::endl ;
         else
             out << ", " ;
    }
    out << std::endl ;
    out << " coordy = " ;
    out << "  " ;
    for(size_t i = 0 ; i < nodes.size()/planes ; i++)
    {
         if(i%10 == 0 && i > 0)
         {
             out << std::endl ;
             out << "  " ;
         }
         out << nodes[i]->getY()  ;
         if(i == (nodes.size()/planes)-1)
             out << " ;" << std::endl ;
         else
             out << ", " ;
    }

    out << "}" << std::endl ;


    if(ncgen.size() > 0)
    {
        std::string command = ncgen + " -o " + filename + ".e " + filename ;
        std::system(command.c_str()) ;
    }

}

}



























