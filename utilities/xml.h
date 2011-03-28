//
// C++ Interface: xml import/export
//
// Description: 
//
//
// Author: Alain Giorla <alain.giorla@epfl.ch>, (C) 2010-2011
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef __XML_H__
#define __XML_H__

#include "matrixops.h"
#include <string>
#include <sstream>

namespace Mu
{

class Point ;
  
class XMLTree
{
protected:
	/** \brief the XML balise. */
	std::string descriptor ;
	/** \brief the children in the tree. */
	std::vector<XMLTree *> children ;
	/** \brief the numerical values of the field. This tree works better if there is either children or fields, but not both. */
	std::vector<double> fields ;

public:
	/** \brief copy a XML tree item. */
	XMLTree(XMLTree * item) ;

	/** \brief void XML tree item. */
	XMLTree(std::string d) ;

	/** \brief XML tree item containing one value. */
	XMLTree(std::string d, int x) ;

	/** \brief XML tree item containing one value. */
	XMLTree(std::string d, double x) ;

	/** \brief XML tree item containing one matrix.
	* The data structure is: <d><matrix> (number of rows) (number of columns) (all values) </matrix></d>
	*/
	XMLTree(std::string d, Matrix m) ;

	/** \brief XML tree item containing one vector.
	* The data structure is: <d><vector> (size) (all values) </vector></d>
	*/
	XMLTree(std::string d, Vector v) ;

	/** \brief XML tree item containing one vector.
	* The data structure is: <d><vector> (size) (all values) </vector></d>
	*/
	XMLTree(std::string d, std::vector<int> v) ;

	/** \brief XML tree item containing one vector.
	* The data structure is: <d><vector> (size) (all values) </vector></d>
	*/
	XMLTree(std::string d, std::vector<double> v) ;

	/** \brief XML tree item containing a pair.
	* The data structure is: <d><vector> (size = 2) (all values) </vector></d>
	*/
	XMLTree(std::string d, std::pair<double, double> p) ;

	/** \brief add one child to the tree. */
	void addChild(XMLTree * item) ;

	/** \brief add children to the tree. */
	void addChild(std::vector<XMLTree * > items) ;

	/** \brief add one child (new tree item) to the tree. */
	void addChild(std::string d) ;

	/** \brief add a value in the field. */
	void addField(int val) ;

	/** \brief add values in the field. */
	void addField(std::vector<int> val) ;

	/** \brief add a value in the field. */
	void addField(double val) ;

	/** \brief add values in the field. */
	void addField(std::vector<double> val) ;

	/** \brief add values in the field. */
	void addField(Vector val) ;

	/** \brief links a tree together, knowing only the position of each item in the tree*/
	void link(std::vector<std::pair<int, XMLTree * > > tree, int i) ;

	/** \brief returns <d> */
	std::string open() ;

	/** \brief returns </d> */
	std::string close() ;

	/** \brief returns <d></d> */
	std::string openclose() ;

	/** \brief returns <d>(all values)</d> */
	std::string line() ;

	std::string getDescriptor() {return descriptor ; } ;
	std::vector<XMLTree *> getChildren() {return children ; } ;
	std::vector<double> getFields() {return fields ; } ;
	XMLTree * getChild(int i) {return children[i] ; } ;
	double getField(int i) {return fields[i] ; } ;
	size_t nChildren() {return children.size() ; } ;
	size_t nFields() {return fields.size() ; } ;
	bool hasChildren() {return children.size() > 0 ; } ;
	bool hasFields() {return fields.size() > 0 ; } ;

	/** \brief checks if this tree item has (balise) as descriptor */
	bool match(std::string balise) ;

	/** \brief checks if thwo tree items have the same descriptor */
	bool match(XMLTree * xml) ;

	/** \brief checks if (balise) is a descriptor of one of the children */
	bool isChild(std::string balise) ;

	/** \brief checks if a tree item is a children */
	bool isChild(XMLTree * xml) ;


	/** \brief prints as string (use verbose = true to actually display results in console) */
	std::vector<std::string> print(bool verbose) ;

	/** \brief prints as string, and indent at level n */
	std::vector<std::string> print(int n) ;

	/** \brief prints into a file (the ".xml" extension is automatically added) */
	std::vector<std::string> printInFile(std::string filename) ;

	/** \brief displays the descriptor */
	void printBalise() ;

	/** \brief returns the first field value (if exists) */
	std::pair<bool, double> buildDouble() ;

	/** \brief returns second and third field values as a pair (if exists) ;  */
	std::pair<bool, std::pair<double, double> > buildPair() ;

	/** \brief returns second and plus field values as a vector (if exists) ;  */
	std::pair<bool, std::vector<double> > buildVector() ;

	/** \brief returns third and plus field values as a matrix (if exists) ;  */
	std::pair<bool, Matrix> buildMatrix() ;



} ;











std::string tostring(int i) ;
std::string tostring(double d) ;
std::vector<std::string> appendStringVector(std::vector<std::string> first, std::vector<std::string> second) ;
std::string indent(int n) ;

XMLTree * importXMLFile(std::string filename) ;

}

#endif
