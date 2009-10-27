
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef __CONFIGURATION_H_
#define __CONFIGURATION_H_

#include <vector>
#include <assert.h>
#include <string>
#include "../geometry/geometry_2D.h"

enum  Type
{
	TRUNK,
	DOUBLE,
	INCLUSION,
	CRACK,
	TYPE_PROBLEME
} ;


/** Configuration atom for the configuration tree.
 * 
 */
struct ConfigAtom
{
	/**  Parent coniguration atom. */
	ConfigAtom * parent ;
	
	/**  Children coniguration atom. */
	std::vector< ConfigAtom * > children ;
	
	/** Label for the configuration entry. Serves also as an identifyier.*/
	std::string label ;
	
	/** Returns a human-readable form of the type of object that is configured.*/
	std::string translate(Type t) ;
	
	/** Type of configurable entity described.*/
	Type type ;
	ConfigAtom() { } ;
	virtual ~ConfigAtom() { } ;
	
	/** Accessor. Is coded for consistency, despite the fact that parent is a public entity.*/
	ConfigAtom * getParent();
	
	/** Return the label of this entity */
	std::string getLabel() ;
	
	/** Return the label of this entity with all its ascendency.
	  * It is written in the form grandparent.parent.self .
	  */
	std::string getFullLabel() ;
	
	/** Set the label for this atom. */
	void setLabel(std::string label) ;
} ;


/** \c Configlet is a a member of the configuration tree. 
 * It contains data and is identified by its labal and position in the tree.
 */
template <typename T>
class Configlet : public ConfigAtom
{
protected:
	/** The actual data payload.*/
	T * data ;
	
public:
	
	/** Construct a configlet for a given object.
	 * 
	 * @param p parent.
	 * @param data data payload of the configlet.
	 * @param type type of the configured object.
	 */
	Configlet(ConfigAtom * p, T * data, Type type) ;
	
	/** Construct a configlet for a given object.
	 * 
	 * @param p parent.
	 * @param data data payload of the configlet.
	 * @param type type of the configured object.
	 * @param label set the label of the configlet.
	 */
	Configlet(ConfigAtom * p, T * data, Type type, std::string label);
	~Configlet() ;
	
	template <typename T2>
		/** Add to the datatree from this configlet.
		 * 
		 * @param data new data.
		 * @param type type of the sub-object.
		 * @return a pointer to the newly created Configlet.
		 */
		Configlet<T2> * putData( T2 * data, Type type) ;
	template <typename T2>
		/** Add to the datatree from this configlet.
		 * 
		 * @param data new data.
		 * @param type type of the sub-object.
		 * @param label set the label of the new Configlet.
		 * @return a pointer to the newly created Configlet.
		 */
		Configlet<T2> * putData( T2 * data, Type type, std::string label) ;
	
	void setData(T * data) ;
	
	T * getData() ;
} ;


template <typename T>
Configlet<T>::Configlet(ConfigAtom * _parent, T * _data, Type _type)
{
	if ((_parent == NULL) && (_type != TRUNK))
	{
		this->type = TRUNK ;
		this->data = NULL ;
		putData(_data, _type) ;
		parent = NULL ;
		label = std::string("trunk") ;
	}
	else
	{
		type = _type ;
		data = _data ;
		label = translate(_type) ;
		parent = _parent ;
	}
}

template <typename T>
Configlet<T>::Configlet(ConfigAtom * _parent, T * _data, Type _type, std::string _label)
{
	if ((_parent == NULL) && (_type != TRUNK))
	{
		this->type = TRUNK ;
		this->data = NULL ;
		parent = NULL ;
		putData(_data, _type, _label) ;
		label = "trunk" ;
	}
	else
	{
		type = _type ;
		data = _data ;
		label = _label ;
		parent = _parent ;
	}
}


template <typename T>
Configlet<T>::~Configlet()
{
	typedef typename std::vector<ConfigAtom *>::iterator VI ;
	
	for(VI i = children.begin() ; i < children.end() ; i++)
		delete &i ;
	
	//delete data ;
}

template <typename T>
template <typename T2>
Configlet<T2> * Configlet<T>::putData( T2 * _data, Type _type)
{	
	Configlet<T2> *c = new Configlet<T2>(this, _data, _type) ;
	return c ;
}

template <typename T>
template <typename T2>
Configlet<T2> * Configlet<T>::putData( T2 * _data, Type _type, std::string _label)
{
	Configlet<T2> *c = new Configlet<T2>( this, _data, _type, _label) ;
	return c ;
}

template <typename T>
void Configlet<T>::setData(T * _data)
{
	if(type == TRUNK)
		return ;
	
	assert(data !=NULL) ;
	data = _data ;
}


template <typename T>
T * Configlet<T>::getData()
{
	return data ;
}



#endif // __CONFIGURATION_H_
