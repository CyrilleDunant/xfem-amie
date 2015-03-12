//
// C++ Interface: voxelfilter
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2007-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef MUVOXELFILTER_H
#define MUVOXELFILTER_H

#include <vector>
#include <map>

#include "../elements/elements.h"
#include "../mesher/delaunay_3d.h"

namespace Amie {

/** \brief Voxel reader from µic output
	@author Cyrille Dunant <cyrille.dunant@epfl.ch>
*/
class VoxelFilter{
protected:
	std::vector<Point *> points ;
	std::vector<DelaunayTetrahedron *> elems;

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

/** \brief Constructor*/
	VoxelFilter();

	virtual ~VoxelFilter();
	
/** \brief read THe file containing the voxel information*/
	void read(const char * filename, Mesh<DelaunayTetrahedron, DelaunayTreeItem3D> * mesh = nullptr) ;
	
/** \brief assign voxel values to Behaviours*/
	std::map<unsigned char,Form *> behaviourMap ;
	
/** \brief return the points forming the mesh*/
	std::vector<Point *> & getPoints() ;

/** \brief return the elements forming the mesh*/
	std::vector<DelaunayTetrahedron *> & getElements() ;

	bool existsPath(std::vector<std::vector<std::vector<unsigned char> > > & phase,
	                int isource, int jsource, int ksource,
	                int itarget, int jtarget, int ktarget,
			int istart, int jstart, int kstart,
	                int iend, int jend, int kend
	               ) const ;
	
/** \brief return the points forming the mesh*/
	const std::vector<Point *> & getPoints() const ;

/** \brief return the elements forming the mesh*/
	const std::vector<DelaunayTetrahedron *> & getElements() const ;
    void update(std::vector< DelaunayTetrahedron* > & tree, const char* filename, Mesh <DelaunayTetrahedron, DelaunayTreeItem3D >* mesh = nullptr);

};

}

#endif //MUVOXELFILTER_H
