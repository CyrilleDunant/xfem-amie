//
// C++ Interface: xml import/export
//
// Description: 
//
//
// Author: Alain Giorla
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "xml.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>


namespace Mu
{
XMLTree::XMLTree(XMLTree * item)
{
	this->descriptor = item->getDescriptor() ;
	for(size_t i = 0 ; i < item->nChildren() ; i++)
		this->addChild(new XMLTree(item->getChild(i))) ;
	for(size_t i = 0 ; i < item->nFields() ; i++)
		this->addField(item->getField(i)) ;
}

XMLTree::XMLTree(std::string d)
{
	this->descriptor = d ;
}

XMLTree::XMLTree(std::string d, int x)
{
	this->descriptor = d ;
	this->addField(x) ;
}

XMLTree::XMLTree(std::string d, double x)
{
	this->descriptor = d ;
	this->addField(x) ;
}

XMLTree::XMLTree(std::string d, Matrix m)
{
	this->descriptor = d ;

	XMLTree * mat = new XMLTree("matrix") ;
	mat->addField((int) m.numRows()) ;
	mat->addField((int) m.numCols()) ;
	mat->addField(m.array()) ;

	this->addChild(mat) ;
}

XMLTree::XMLTree(std::string d, Vector v)
{
	this->descriptor = d ;

	XMLTree * vec = new XMLTree("vector") ;
	vec->addField((int) v.size()) ;
	vec->addField(v) ;

	this->addChild(vec) ;
}

XMLTree::XMLTree(std::string d, std::vector<double> v)
{
	this->descriptor = d ;

	XMLTree * vec = new XMLTree("vector") ;
	vec->addField((int) v.size()) ;
	vec->addField(v) ;

	this->addChild(vec) ;
}

XMLTree::XMLTree(std::string d, std::vector<int> v)
{
	this->descriptor = d ;

	XMLTree * vec = new XMLTree("vector") ;
	vec->addField((int) v.size()) ;
	vec->addField(v) ;

	this->addChild(vec) ;
}

XMLTree::XMLTree(std::string d, std::pair<double, double> p)
{
	this->descriptor = d ;

	XMLTree * pr = new XMLTree("vector") ;
	pr->addField(2) ;
	pr->addField(p.first) ;
	pr->addField(p.second) ;

	this->addChild(pr) ;
}

void XMLTree::addChild(XMLTree * item)
{
	this->children.push_back(item) ;
}

void XMLTree::addChild(std::vector<XMLTree *> items)
{
	for(size_t i = 0 ; i < items.size() ; i++)
		this->children.push_back(items[i]) ;
}

void XMLTree::addChild(std::string d)
{
	this->children.push_back(new XMLTree(d)) ;
}

void XMLTree::addField(int val)
{
	this->fields.push_back((double) val) ;
}

void XMLTree::addField(std::vector<int> val)
{
	for(size_t i = 0 ; i < val.size() ; i++)
		this->addField(val[i]) ;
}

void XMLTree::addField(double val)
{
	this->fields.push_back(val) ;
}

void XMLTree::addField(std::vector<double> val)
{
	for(size_t i = 0 ; i < val.size() ; i++)
		this->addField(val[i]) ;
}

void XMLTree::addField(Vector val)
{
	for(size_t i = 0 ; i < val.size() ; i++)
		this->addField(val[i]) ;
}

void XMLTree::link(std::vector<std::pair<int, XMLTree * > > tree, int i)
{
	if(this->match(tree[i].second))
	{
		int level = tree[i].first ;
		int current = level ;
		size_t j = i + 1 ;
		while(j < tree.size())
		{
			current = tree[j].first ;
			if(current == level + 1)
			{
				tree[j].second->link(tree, j) ;
				this->addChild(tree[j].second) ;
			}
			if(current < level + 1)
				j = tree.size() ;
			j++ ;
		}
	}
}


std::string XMLTree::open()
{
	std::string xml = "<" ;
	xml.append(this->descriptor) ;
	xml.append(">") ;
	return xml ;
}

std::string XMLTree::close()
{
	std::string xml = "</" ;
	xml.append(this->descriptor) ;
	xml.append(">") ;
	return xml ;
}

std::string XMLTree::openclose()
{
	std::string o = this->open() ;
	std::string c = this->close() ;
	return o.append(c) ;
}

std::string XMLTree::line()
{
	std::string o = this->open() ;
	std::string m = " " ;
	for(size_t i = 0 ; i < nFields() ; i++)
	{
		m = m.append(tostring(fields[i])) ;
		m = m.append(" ") ;
	}
	std::string c = this->close() ;
	return o.append(m).append(c) ;
}

bool XMLTree::match(std::string balise)
{
	int test = strcmp(descriptor.c_str(),balise.c_str()) ;
	return (test == 0) ;
}

bool XMLTree::match(XMLTree * xml)
{
	std::string balise = xml->getDescriptor() ;
	return this->match(balise) ;
}

bool XMLTree::isChild(std::string balise)
{
	if(this->match(balise))
		return true ;
	size_t i = 0 ;
	while(i < children.size())
	{
		if(children[i]->isChild(balise))
			return true ;
		i++ ;
	}
	return false ;
}

bool XMLTree::isChild(XMLTree * xml)
{
	if(this->match(xml))
		return true ;
	size_t i = 0 ;
	while(i < children.size())
	{
		if(children[i]->isChild(xml))
			return true ;
		i++ ;
	}
	return false ;
}

std::vector<std::string> XMLTree::print(bool verbose)
{
	std::vector<std::string> xml = this->print(0) ;
	if(verbose)
	{
		for(size_t i = 0 ; i < xml.size() ; i++)
			std::cout << xml[i] << std::endl ;
	}
	return xml ;
}

std::vector<std::string> XMLTree::print(int n)
{
	std::vector<std::string> xml ;
	if(children.size() == 0 && fields.size() == 0)
	{
		xml.push_back(indent(n).append(this->openclose())) ;
		return xml ;
	}
	if(children.size() == 0)
	{
		xml.push_back(indent(n).append(this->line())) ;
		return xml ;
	}
	xml.push_back(indent(n).append(this->open())) ;	
	if(children.size() > 0)
	{
		for(size_t i = 0 ; i < children.size() ; i++)
			xml = appendStringVector(xml, children[i]->print(n+1)) ;
	}
	if(fields.size() > 0)
		xml.push_back(indent(n+1).append(this->line())) ;
	xml.push_back(indent(n).append(this->close())) ;
	return xml ;
}

std::vector<std::string> XMLTree::printInFile(std::string filename)
{
	std::cout << "exporting to xml file..." << std::endl ;
	std::vector<std::string> xml = this->print(false) ;
	std::fstream parser ;
	std::string fullfilename = filename.append(".xml") ;
	parser.open(fullfilename.c_str(), std::ios::out) ;
	for(size_t i = 0 ; i < xml.size() ; i++)
		parser << xml[i] << "\n" ;
	parser.close() ;
	std::cout << "...done" << std::endl ;
	return xml ;
}

void XMLTree::printBalise()
{
	std::cout << descriptor << std::endl ;
	if(children.size() > 0)
	{
		for(size_t i = 0 ; i < children.size() ; i++)
			children[i]->printBalise() ;
	}
}

std::pair<bool, double> XMLTree::buildDouble()
{
	bool found = false ;
	double db = 0 ;
	if(this->hasFields())
	{
		found = true ;
		db = fields[0] ;
	}
	return std::make_pair(found, db) ;
}

std::pair<bool, std::pair<double, double> > XMLTree::buildPair()
{
	bool found = false ;
	std::pair<double,double> pr = std::make_pair(0,0) ;
	std::pair<bool, double> vfirst = std::make_pair(false,0) ;
	std::pair<bool, double> vsecond = std::make_pair(false,0) ;
	if(this->match("vector"))
	{
		if(this->nFields() > 2)
		{
			found = true ;
			pr = std::make_pair(fields[1],fields[2]) ;
		}
	}
	return std::make_pair(found,pr) ;
}

std::pair<bool, std::vector<double> > XMLTree::buildVector()
{
	bool found = false ;
	std::vector<double> val ;
	if(this->match("vector"))
	{
		if(this->nFields() > 1)
		{
			size_t s = (size_t) fields[0] ;
			if(s == this->nFields() -1)
				found = true ;
			for(size_t i = 1 ; i < s+1 && i < this->nFields() ; i++)
				val.push_back(fields[i]) ;
		}
	}
	return std::make_pair(found,val) ;
}

std::pair<bool, Matrix> XMLTree::buildMatrix()
{
	bool found = false ;
	Matrix m(1,1) ;
	if(this->match("matrix"))
	{
		if(this->nFields() > 2)
		{
			size_t r = fields[0] ;
			size_t c = fields[1] ;
			if(r*c == this->nFields() -2)
				found = true ;
			m.resize(r,c) ;
			for(size_t i = 0 ; i < r ; i++)
			{
				for(size_t j = 0 ; j < c ; j++)
				{
					if(j*r + i < this->nFields() -2)
						m[i][j] = fields[2+j*r+i] ;
				}
			}
		}
	}
	return std::make_pair(found,m) ;
}




















std::string tostring(int a)
{
	std::stringstream buf ;
	buf << a ;
	return buf.str() ;	
}

std::string tostring(double d)
{
	std::stringstream buf ;
	buf << d ;
	return buf.str() ;	
}

std::vector<std::string> appendStringVector(std::vector<std::string> first, std::vector<std::string> second)
{
	std::vector<std::string> out ;
	for(size_t i = 0 ; i < first.size() ; i++)
		out.push_back(first[i]) ;
	for(size_t i = 0 ; i < second.size() ; i++)
		out.push_back(second[i]) ;
	return out ;
}

std::string indent(int n)
{
	std::string i = "" ;
	if(n < 1)
		return i ;
	i = "  " ;
	return i.append(indent(n-1)) ;
}







std::pair<int,std::string> opener(std::string balise)
{
	std::string seq = "<" ;
	std::string seqclose = "</" ;
	std::string geq = ">" ;
	std::string descriptor = balise ;
	int op = 0 ;
	int length = balise.size() ;
	if(balise.compare(0,1,seq)==0 && balise.compare(length-1,1,geq)==0)
	{
		op = 1 ;
		descriptor = balise.substr(1,length-2) ;
		if(balise.compare(0,2,seqclose) == 0)
		{
			op = -1 ;
			descriptor = balise.substr(2,length-3) ;
		}
	}
	return std::make_pair(op,descriptor) ;
}

XMLTree * importXMLFile(std::string filename)
{
	std::cout << "importing from xml file..." << std::endl ;
	std::string buffer ;

	std::vector<std::pair<int, XMLTree * > > tree ;
	int level = 0 ;
	size_t last_at_level = 0 ;
	std::vector<double> val ;
	bool previous_is_value = false ;
	std::pair<int, std::string> op ;

	std::fstream gway ;
	gway.open(filename.c_str(), std::ios::in) ;
	while(!gway.eof())
	{
		gway >> buffer ;
		op = opener(buffer.c_str()) ;
		switch(op.first)
		{
			case -1 :
			{
				level-- ;
				if(previous_is_value)
				{
					last_at_level = 0 ;
					for(size_t i = 0 ; i < tree.size() ; i++)
					{
						if(tree[i].first == level)
							last_at_level = i ;
					}
					tree[last_at_level].second->addField(val) ;
				}
				val.clear() ; 
				previous_is_value = false ;
				break ;
			}
			case 0 :
			{
				std::cout << op.second <<std::endl ;
				val.push_back(atof(op.second.c_str())) ;
				previous_is_value = true ;
				break ;
			}
			case 1 :
			{
				tree.push_back(std::make_pair(level,new XMLTree(op.second))) ;
				level++ ;
				val.clear() ;
				previous_is_value = false ;
				break ;
			}
		}
	}
	gway.close() ;

	tree[0].second->link(tree,0) ;

	return tree[0].second ;
}






























}
