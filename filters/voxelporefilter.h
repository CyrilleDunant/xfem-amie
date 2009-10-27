//
// C++ Interface: voxelfilter
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUVOXELPOREFILTER_H
#define MUVOXELPOREFILTER_H

#include <vector>
#include <map>
#include <set>

#include "../elements/elements.h"
#include "../mesher/delaunay_3d.h"

namespace Mu {

/** \brief Read voxel information from µic. All voxels with value 0 will be used to form a pore network
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class VoxelPoreFilter{
protected:
	std::vector<Point *> points ;
	std::vector<DelaunayTetrahedron *> elems;
	
	bool existsPath(std::vector<std::vector<std::vector<int> > > & phase,
	                int isource, int jsource, int ksource,
	                int itarget, int jtarget, int ktarget,
			int istart, int jstart, int kstart,
	                int iend, int jend, int kend
	               ) const ;
	
	struct ConnectedNode
	{
		int i ; 
		int j ; 
		int k ;
		bool visited ;
		ConnectedNode(int i_, int j_, int k_) : i(i_), j(j_), k(k_), neighbour(0) { visited = false ;} ;
		
		std::vector<ConnectedNode *> neighbour ;
		
		bool isNeighbour(const ConnectedNode * n) const
		{
			return std::abs(i-n->i) + std::abs(j-n->j) + std::abs(k-n->k) < 2 ;
		}
	} ;
	
public:
	VoxelPoreFilter();

	LinearForm * behaviour ;
	virtual ~VoxelPoreFilter();
	
/** \brief read The file containing the voxel information*/
	void read(const char * filename) ;
	
	int poreIndex ;
	
/** \brief return the points forming the mesh*/
	std::vector<Point *> & getPoints() ;

/** \brief return the elements forming the mesh*/
	std::vector<DelaunayTetrahedron *> & getElements() ;
	
/** \brief return the points forming the mesh*/
	const std::vector<Point *> & getPoints() const ;

/** \brief return the elements forming the mesh*/
	const std::vector<DelaunayTetrahedron *> & getElements() const ;

};

}

#endif //MUVOXELFILTER_H
