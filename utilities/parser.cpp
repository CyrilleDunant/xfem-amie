// Author: Cyrille Dunant <cyrille.dunant@gmail.com>, (C) 2005-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "parser.h"

void InclusionParser::readData()
{
	if (!file.fail())
	{
		circs->clear();
		size_t numberOfCircles ;
		size_t index ;
		file >> numberOfCircles ;
		
		for (size_t i = 0 ; i < numberOfCircles ; i++)
		{
			double centerX, centerY, radius ;
			file >> index >> radius >> centerX >> centerY ;
	#ifdef DEBUG
			std::cout << "index : " << index << ", r : " << radius << ", x : " << centerX << ", y : " << centerY << std::endl ;
	#endif
			circs->push_back(Circle(radius, centerX, centerY )) ;
		}
	}
}

void SegmentsParser::readData()
{
	if(!file.fail())
	{
		segments->clear() ;
		size_t numberOfSegments ;
		size_t index ;
		file >> numberOfSegments ;
		
		for (size_t i = 0 ; i < numberOfSegments ; i++)
		{
			if( file.eof())
				return ;
			
			size_t numberOfHeads, numberOfPoints ;
			file >> index >> numberOfHeads >> numberOfPoints ;
			std::valarray<Point *> m(numberOfPoints) ;
			
			for ( size_t j = 0 ; j < numberOfPoints ; j++)
			{
				double x, y ;
				file >> x >> y ;
				m[j] = new Point(x,y) ;

			}
			segments->push_back(SegmentedLine(m)) ;
		}
	}
}

void PointParser::readData()
{
	if(!file.fail())
	{
		pointSets->clear() ;
		
		size_t numberOfPoints, numberOfPointSets ;
		size_t index ;
		file >> numberOfPointSets ;
		
		for (size_t i = 0 ; i < numberOfPointSets ; i++)
		{
			file >> numberOfPoints ;
			PointSet * m = new PointSet(numberOfPoints) ;
			for (size_t j = 0 ;  j < numberOfPoints ; j++)
			{
				if( file.eof())
					return ;
				
				file >> index >> m->getPoint(j)->getX() >> m->getPoint(j)->getY() ;
	#ifdef DEBUG
			std::cout << "index : " << index << ", x : "<< m->getPoint(j)->getX() << ", y : " << m->getPoint(j)->getY() << std::endl ;
	#endif		
			}
			pointSets->push_back(m) ;
		}
	}
}



void PointParser3D::readData()
{
	if(!file.fail())
	{
		pointSets->clear() ;
		
		size_t numberOfPoints;
		size_t index ;
		size_t indexPeriodic;
		file >> numberOfPoints ;
		for (size_t i = 0 ;  i < numberOfPoints ; i++)
		{
			if( file.eof())
				return ;
				
			double x,y,z ;
//  			file >>  index>> x >> y >> z ;
			file >> index >> indexPeriodic >> x >> y >> z;
			Point * p = new Point(x,y,z) ; p->getId() = index-1 ;
			pointSets->push_back(p) ;
			this->PeriodicIds.push_back(indexPeriodic-1) ;

	#ifdef DEBUG
//   			std::cout << "index : " << index << ", x : "<< p->getX() << ", y : " << p->getY() << ", z : " << p->getZ() << std::endl ;
 	    	std::cout << "index : " << index << " indexPeriodic: " << indexPeriodic << ", x : "<< p->getX() << ", y : " << p->getY() << ", z : " << p->getZ() << std::endl ;
	#endif		
		}
	}

}

void PeriodParser::readData()
{
	if(!file.fail())
	{
		periodSets->clear() ;
		
		size_t numberOfindex;
		
		file >> numberOfindex ;
		 std::vector<Point *>  p;
		
		for (size_t i = 0 ;  i < numberOfindex ; i++)
		{
			if( file.eof())
				return ;

			std::valarray<int> index_id(1) ;	
			size_t index ;
			file >> index >> index_id[0] ;
			
			periodSets->push_back(std::pair<size_t,std::valarray<int> >(index,index_id )) ;
	#ifdef DEBUG
			std::cout << "index : " << index << "id "<< index_id[0]<< std::endl ;
	#endif
		}
	}

}
void TetrahedronParser::readData()
{
	if(!file.fail())
	{
		tets->clear() ;
		size_t numberOfTets ;
		file >> numberOfTets ;
		for (size_t i = 0 ;  i < numberOfTets ; i++)
		{
			if( file.eof())
				return ;
			
			std::valarray<int> vertex_id(4) ;
			size_t mat_index ;
			size_t index ;
			file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3] >>  mat_index;
			vertex_id[0]-- ; vertex_id[1]-- ; vertex_id[2]-- ; vertex_id[3]-- ;
			
			tets->push_back(std::pair<std::valarray<int>, size_t>(vertex_id, mat_index)) ;
	#ifdef DEBUG
		std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  << ", point 3 :" << vertex_id[3]  << ", material : "<< mat_index<< std::endl ;
	#endif		
		}
	}
}
void HexahedronParser::readData()
{
	if(!file.fail())
	{
		hex->clear() ;
		size_t numberOfHex ;
		file >> numberOfHex ;
		for (size_t i = 0 ;  i < numberOfHex ; i++)
		{
			if( file.eof())
				return ;
			
			std::valarray<int> vertex_id(8) ;
			size_t mat_index ;
			size_t index ;
			file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3] >> vertex_id[4] >> vertex_id[5] >> vertex_id[6]  >> vertex_id[7] >> mat_index;
			vertex_id[0]-- ; vertex_id[1]-- ; vertex_id[2]-- ; vertex_id[3]-- ; vertex_id[4]-- ;vertex_id[5]-- ;vertex_id[6]-- ;vertex_id[7]-- ;
			
			hex->push_back(std::pair<std::valarray<int>, size_t>(vertex_id, mat_index)) ;
	#ifdef DEBUG
		std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  << ", point 3 :" << vertex_id[3]  << ", point 4 :" << vertex_id[4]  << ", point 5 :" << vertex_id[5]  << ", point 6 :" << vertex_id[6]  << ", point 7 :" << vertex_id[7]  << ", material : "<< mat_index<< std::endl ;
	#endif		
		}
	}
}
void BoundaryParser::readData()
{
	if(!file.fail())
	{
		tri0->clear() ;
		tri1->clear() ;
		tri2->clear() ;
		tri3->clear() ;
		tri4->clear() ;
		tri5->clear() ;
		size_t numberOfTris ;
		file >> numberOfTris ;
		for (size_t i = 0 ;  i < numberOfTris ; i++)
		{
			if( file.eof())
				return ;
			
			std::valarray<int> vertex_id(3) ;
			size_t face_index ;
			size_t index ;
			size_t to_discard ;
			file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >>  face_index >> to_discard ;
			vertex_id[0]-- ; vertex_id[1]-- ; vertex_id[2]-- ;
			if(face_index == 1)
			{
				tri0->push_back(vertex_id) ;
			}
			else if(face_index == 2)
			{
				tri1->push_back(vertex_id) ;
			}
			else if(face_index == 3)
			{
				tri2->push_back(vertex_id) ;
			}
			else if(face_index == 4)
			{
				tri3->push_back(vertex_id) ;
			}
			else if(face_index == 5)
			{
				tri4->push_back(vertex_id) ;
			}
			else if(face_index == 6)
			{
				tri5->push_back(vertex_id) ;
			}
			else
			{
				assert(false) ;
			}
	#ifdef DEBUG
			std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2]  <<  ", face : "<< face_index<< std::endl ;
	#endif		
		}
	}
}

std::vector<std::valarray<int> > * BoundaryParser::getData(size_t face_index)
{
	if(face_index == 0)
	{
		return tri0 ;
	}
	else if(face_index == 1)
	{
		return tri1 ;
	}
	else if(face_index == 2)
	{
		return tri2 ;
	}
	else if(face_index == 3)
	{
		return tri3 ;
	}
	else if(face_index == 4)
	{
		return tri4 ;
	}
	else if(face_index == 5)
	{
		return tri5 ;
	}
	else
	{
		assert(false) ;
		return nullptr ; //shut up the compiler ;
	}
}

void HexahedronBoundaryParser::readData()
{
	if(!file.fail())
	{
		carre0->clear() ;
		carre1->clear() ;
		carre2->clear() ;
		carre3->clear() ;
		carre4->clear() ;
		carre5->clear() ;
		size_t numberOfCarres ;
		file >> numberOfCarres ;
		for (size_t i = 0 ;  i < numberOfCarres ; i++)
		{
			if( file.eof())
				return ;
			
			std::valarray<int> vertex_id(4) ;
			size_t face_index ;
			size_t index ;
			size_t to_discard ;
			file >> index >> vertex_id[0] >> vertex_id[1] >> vertex_id[2] >> vertex_id[3]>>  face_index >> to_discard ;
			vertex_id[0]-- ; vertex_id[1]-- ; vertex_id[2]-- ; vertex_id[3]-- ;
			if(face_index == 1)
			{
				carre0->push_back(vertex_id) ;
			}
			else if(face_index == 2)
			{
				carre1->push_back(vertex_id) ;
			}
			else if(face_index == 3)
			{
				carre2->push_back(vertex_id) ;
			}
			else if(face_index == 4)
			{
				carre3->push_back(vertex_id) ;
			}
			else if(face_index == 5)
			{
				carre4->push_back(vertex_id) ;
			}
			else if(face_index == 6)
			{
				carre5->push_back(vertex_id) ;
			}
			else
			{
				assert(false) ;
			}
	#ifdef DEBUG
			std::cout << "index : " << index << ", point 0 : "<< vertex_id[0]  << ", point 1 : " << vertex_id[1]  << ", point 2 : " << vertex_id[2] << ", point 3 : " << vertex_id[3] <<  ", face : "<< face_index<< std::endl ;
	#endif		
		}
	}
}

std::vector<std::valarray<int> > * HexahedronBoundaryParser::getData(size_t face_index)
{
	if(face_index == 0)
	{
		return carre0 ;
	}
	else if(face_index == 1)
	{
		return carre1 ;
	}
	else if(face_index == 2)
	{
		return carre2 ;
	}
	else if(face_index == 3)
	{
		return carre3 ;
	}
	else if(face_index == 4)
	{
		return carre4 ;
	}
	else if(face_index == 5)
	{
		return carre5 ;
	}
	else
	{
		assert(false) ;
		return nullptr ; //shut up the compiler ;
	}
}

int ConfigParser::getIndentLevel( std::string test ) 
{
	int i = 0 ;
	while(test[i] == '.' && i < test.size())
		i++ ;
	return i ;
}

void ConfigParser::readData()
{
	ConfigTreeItem * current = trunk ;
	int level = 1 ;
	std::string buffer ;
	if(!file.fail())
	{
		for( std::string line ; getline( file, line ) ; ) 
		{
			if(line[0] == '#')
				continue ;

			size_t comment = line.find("#") ;
			if(comment != std::string::npos)
			{
				line = line.substr(0, comment-1) ;
//				std::cout << line << std::endl ;
			}	

			size_t found = line.find(" ") ;
			size_t quote = line.find('"') ;
			while(found != std::string::npos && found < quote)
			{
				line = line.erase(found,1) ;
				found = line.find(" ") ;
			}
			quote = line.find('"') ;
			while(quote != std::string::npos)
			{
				line = line.erase(quote,1) ;
				quote = line.find('"') ;
			}

			size_t sep = line.find("=") ;
			int l = ConfigParser::getIndentLevel(line) ;
			if(l == level)
			{
				// current = current ;
			} 
			else if(l == level+1)
			{
				current = current->getLastChild() ;
				level++ ;
			}
			else if(l == 0 && line.length() == 0)
				break ;
			else if(l == 0 || l > level+1)
			{
				std::cout << "parsing error in file " << filename << std::endl ;
				current->print() ;
				exit(0) ;
			}
			else // (l between 1 and level-1)
			{
				while(level > l)
				{
					level-- ;
					current = current->getFather() ;
				}
			}
			if(sep == std::string::npos)
			{
//				std::cout << line.substr(level) << std::endl ;
				ConfigTreeItem * cnf = new ConfigTreeItem( current , line.substr(level) ) ;
//				current->printTree() ;
			}
			else
			{
				std::string right = line.substr( sep+1 ) ;
				bool isDouble = (right.find_first_not_of("0123456789.e-") == std::string::npos ) ;
				std::string label = line.substr(0, sep) ;
				label = label.substr(level) ;
				if(isDouble)
				{
					ConfigTreeItem * cnf = new ConfigTreeItem( current , label, atof(right.c_str()) ) ;
				}
				else
				{
/*					if(right.find("@") == 0 && authorizeIncludes)
					{
						right = right.substr(1) ;
						ConfigParser includes( right, false) ;
						includes.readData() ;
						ConfigTreeItem * cnf = includes.getData() ;
						ConfigTreeItem * child = cnf->getChild(label) ;
						current->addChild( child ) ;
						child->setFather( current ) ;
					} 
					else if( label == "include"  && authorizeIncludes)
					{
						ConfigParser includes( right, false) ;
						includes.readData() ;
						ConfigTreeItem * cnf = includes.getData() ;
						std::vector<ConfigTreeItem *> children = cnf->getAllChildren() ;
						for(size_t i = 0 ; i < children.size() ; i++)
						{
							current->addChild( children[i] ) ;
							children[i]->setFather( current ) ;
						}
					}
					else*/
						ConfigTreeItem * cnf = new ConfigTreeItem( current , label, right ) ;			
				}
			}

		}
		
		std::cout << filename << " parsed with success!" << std::endl ;

	}
	else
		std::cout << filename << " not found!" << std::endl ;


#ifdef _WIN32
	if(trunk)
		trunk->makeWindowsPath() ;
#endif

}

ConfigTreeItem * ConfigParser::readFile(std::string f, ConfigTreeItem * def, bool define) 
{
	ConfigParser parser(f) ;	
	parser.readData() ;
	ConfigTreeItem * ret = parser.getData() ;
	if(ret->hasChild("template"))
		ret = ret->getChild("template")->makeTemplate() ;
	if(define)
		ret->define(def, true) ;
	return ret ;
}



