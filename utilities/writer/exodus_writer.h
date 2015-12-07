//
// C++ Interface: triangle writer
//
// Description: writer for 2D triangle file
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __EXODUS_WRITER_H__
#define __EXODUS_WRITER_H__

#include "triangle_writer.h"
#include "../../features/features.h"
#include "../inclusion_family.h"
#include "../itoa.h"

namespace Amie
{

class ExodusTriangleWriter : public TriangleWriter
{
    std::vector<unsigned int> index ;
    std::vector<std::pair<Segment *, unsigned int> > sides ;
    std::string ncgen ;
    std::map<unsigned int, unsigned int> elems ;
    bool us = false ;

public:
    ExodusTriangleWriter( std::string file, FeatureTree * F = nullptr, int t = 0) : TriangleWriter(file, F, t) { } ;

    void setCache( std::vector<unsigned int> & caches ) { index = caches ; }
    void setCache( std::vector<Feature *> feats ) ;
    void setCache( std::vector< std::vector< Feature *> > feats, bool last = true ) ;
    void finalizeCache() ;

    void setncgen(std::string n) { ncgen = n ; } 
    void setUSDateFormat(bool u) { us = u ; }

    virtual void write() ;


} ;



}

#endif
