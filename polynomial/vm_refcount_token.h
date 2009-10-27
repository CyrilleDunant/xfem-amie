//
// C++ Interface: vm_refcount_token
//
// Description: 
//
//
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef REF_COUNTED_TOKEN_H
#define REF_COUNTED_TOKEN_H

#include<iostream>
#include<valarray>

#include "vm_token.h"

namespace Mu
{

/** \brief Internal storage of a reference counted token.
 * This class contains the pointer to the refered token and does the reference counting.
*/
class RefCountedTokenInternal
{
private:
	const Token * t ;
	size_t count ;
	
	/** \brief Copy-contructor. Private to prevent copy
	 * 
	 * @param reference
	 */
	RefCountedTokenInternal( const RefCountedTokenInternal & ) ;

	/** \brief Copy-contructor. Private to prevent copy
	 * 
	 * @param  reference
	 * @return self
	 */
	RefCountedTokenInternal &operator=( const RefCountedTokenInternal & )
	{
		return *this ;
	}
	
public:

	/** Evaluate the token given a Context
	 * 
	 * @param context 
	 */
	void eval(Context &context) const
	{
		t->eval(context) ;
	}
	
	/** \brief Constructor, increment reference count.
	 * 
	 * @param tok Token to point to
	 */
	RefCountedTokenInternal(const Token * tok) 
	{
		t = tok ;
		count = 1 ;
	}
	
	/** \brief Destructor. Also destroys the pointed Token
	 * 
	 */
	~RefCountedTokenInternal()
	{
		delete t ;
	}
	
	/** \brief Cast operator so the class behaves exactly as the Token pointed to
	 * 
	 * @return Token * pointed to
	 */
	const Token * token() const 
	{
		return t ;
	}
	
	/** \brief decrement reference count 
	 * 
	 * @return 0 if count == 0
	 */
	inline bool dereferenced()
	{
		count-- ;
		if (count == 0)
		{
			return true ;
		}
		else
		{
			return false ;
		}
	}
	
	/** \brief increment reference count
	 * 
	 */
	inline void referenced()
	{
		count++ ;
	}
	
	/** \brief accessor, give number of time the Token is referenced
	 * 
	 * @return count
	 */
	size_t getCount() const
	{
		return count ;
	}
} ;

/** \brief Public reference counted token. Behaves like a Token *
*/
class RefCountedToken
{
protected:
	RefCountedTokenInternal * t ; 

	/** \brief Add a reference
	 * 
	 */
	void referenced() 
	{
		if(t!=NULL)
			t->referenced() ;
	}
	
	/** \brief Remove a reference, cleanup if it was the last.
	 * 
	 */
	void dereferenced() 
	{
		if ( t!=NULL && t->dereferenced() )
		{
			delete t;
			t = NULL ;
		}
	}
	
public:
	
	/** \brief Default constructor. Behaves like NULL
	 * 
	 */
	RefCountedToken() : t(NULL) { } ; 
	
	
	/** \brief Evaluate the token given a Context
	 * 
	 * @param context context
	 */
	void eval(Context &context) const
	{
		t->eval(context) ;
	}
	
	/** \brief decrement reference count
	 * 
	 */
	~RefCountedToken()
	{
		dereferenced() ;
	}
	
	/** \brief Constructor. Behaves like the Token given as an argument
	 * 
	 * @param token 
	 */
	RefCountedToken(const Token * token)
	{
		t = new RefCountedTokenInternal(token) ;
	}
	
	/** \brief Copy-constructor. Behaves like the Token given as an argument, and increments reference counts
	 * 
	 * @param token 
	 */
	RefCountedToken(const RefCountedToken & token)
	{
		t = token.t ;
		referenced() ;
	}
	
	/** \brief Copy-constructor. Behaves like the Token given as an argument, and increments reference counts
	 * 
	 * @param token 
	 * @return 
	 */
	RefCountedToken & operator = ( const RefCountedToken & token)
	{
		if(t == NULL)
		{
			t = token.t ;
			referenced() ;
		}
		else
		{
			dereferenced();
			t = token.t ;
			referenced() ;
		}
		return *this ;
	}
	
	/** \brief Copy-constructor. Behaves like the Token given as an argument, and increments reference counts
	 * 
	 * @param token 
	 * @return 
	 */
	RefCountedToken & operator = (const Token * token)
	{
		if ( t && t->getCount() == 1 ) {
			referenced() ;
		}
		else 
		{
			dereferenced();
			t = new RefCountedTokenInternal( token );
		}
		return *this;
	}
	
	/** \brief Compare the Token pointers 
	 * 
	 * @param p to compare
	 * @return true if the pointers are equal
	 */
	bool operator== ( const RefCountedToken& p ) const
	{
		return (const Token *)(*this) == (const Token*) p;
	}

	/** \brief Compare the Token pointers 
	 * 
	 * @param p to compare
	 * @return true if the pointers are different
	 */
	bool operator!= ( const RefCountedToken & p ) const
	{
		return !( *this == p );
	}
	
	/** \brief Behave like the Token pointed to
	 * 
	 * @return The pointer to the Token 
	 */
	const Token * operator-> () const
	{
		return (const Token *)(t?t->token():0) ; 
	} 
	
	/** \brief Behave like the Token pointed to. Cast operator
	 * 
	 */
	operator const Token * () const
	{
		return (const Token *)(t?t->token():0) ; 
	}

} ;

typedef std::valarray<RefCountedToken> ByteCode  ;

} ;

#endif
