//
// C++ Implementation: fracturecriterion
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "artificialfracture.h"
#include "../../mesher/delaunay.h"
#include "../../mesher/delaunay_3d.h"
namespace Mu {

ArtificialFracture::ArtificialFracture( MirrorState mirroring, double delta_x, double delta_y, double delta_z) : FractureCriterion(mirroring, delta_x, delta_y, delta_z)
{
}

ArtificialFracture::~ArtificialFracture()
{
}

}

