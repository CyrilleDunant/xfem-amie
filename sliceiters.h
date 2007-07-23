
// Author: Cyrille Dunant <cyrille.dunant@epfl.ch>, (C) 2005-2007
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef __SLICEITERS_H_
#define __SLICEITERS_H_

#include<valarray>

namespace Mu
{
template<class T>
class Slice_iter
{
	std::valarray<T> *v ;
	std::slice s ;
	size_t curr;
	
	T& ref(size_t i) const {return (*v)[s.start()+i*s.stride()]; }
	
public:
	
Slice_iter(std::valarray<T> *vv, std::slice ss) : v(vv), s(ss), curr(0) { }
	
	Slice_iter end() const
	{
		Slice_iter t = *this ;
		t.curr = s.size();
		return t ;
	}
	
	Slice_iter& operator ++() {curr++; return *this; }
	Slice_iter operator ++(int) {Slice_iter t = *this ; curr++ ; return t; }
	
	T& operator[](size_t i){ return ref(i) ; }
	T& operator[](size_t i) const { return ref(i) ; }
	T& operator()(size_t i){ return ref(i) ; }
	T& operator()(size_t i) const { return ref(i) ; }
	T& operator*(){ return ref(curr) ;}
	
	template < class C >
	friend bool operator ==(const Slice_iter &p, const Slice_iter &q) ;
	template < class C >
	friend bool operator !=(const Slice_iter &p, const Slice_iter &q) ;
	template < class C >
	friend bool operator < (const Slice_iter &p, const Slice_iter &q) ;
} ;


template < class T > bool operator ==(const Slice_iter<T> &p, const Slice_iter<T> &q)
{
	return p.curr==q.curr && p.s.stride()==q.s.stride() && p.s.start() == q.s.start() ;
} ;

template < class T > bool operator !=(const Slice_iter<T> &p, const Slice_iter<T> &q)
{
	return !(p==q) ;
} ;

template < class T > bool operator <(const Slice_iter<T> &p, const Slice_iter<T> &q)
{
	return p.curr < q.curr && p.s.stride() == q.s.stride() && p.s.start() == q.s.start() ;
} ;

template<class T>
class Cslice_iter
{
public:
	std::valarray<T> *v ;
	std::slice s ;
	size_t curr;
	
	T& ref(size_t i) const {return (*v)[s.start()+i*s.stride()]; }
	
public:
	
Cslice_iter(std::valarray<T> *vv, std::slice ss) : v(vv), s(ss), curr(0) { }
	
	Cslice_iter end() const
	{
		Cslice_iter t = *this ;
		t.curr = s.size();
		return t ;
	}
	
	Cslice_iter& operator ++() {curr++; return *this; }
	Cslice_iter operator ++(int) {Cslice_iter t = *this ; curr++ ; return t; }
	
	T& operator[](size_t i) const { return ref(i) ; }
	T& operator()(size_t i) const { return ref(i) ; }
	T& operator*() const { return ref(curr) ;}
	
	template < class C >
	friend bool operator==(const Cslice_iter &p, const Cslice_iter &q) ;
	template < class C >
	friend bool operator!=(const Cslice_iter &p, const Cslice_iter &q) ;
	template < class C >
	friend bool operator< (const Cslice_iter &p, const Cslice_iter &q) ;
	
} ;

template < class T > bool operator==(const Cslice_iter<T> &p, const Cslice_iter<T> &q)
{
	return p.curr==q.curr && p.s.stride()==q.s.stride() && p.s.start() == q.s.start() ;
}

template < class T > bool operator!=(const Cslice_iter<T> &p, const Cslice_iter<T> &q)
{
	return !(p==q) ;
}

template < class T > bool operator<(const Cslice_iter<T> &p, const Cslice_iter<T> &q) 
{
	return p.curr < q.curr && p.s.stride() == q.s.stride() && p.s.start() == q.s.start() ;
}

}

#endif //__SLICEITERS_H_
