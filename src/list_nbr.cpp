// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//-----------------------------------------------------------------------------
// Created on:	12 Apr 2004
// Last update:	25 Feb 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NBRLIST_H
#define CLASS_NBRLIST_H

#include <list>

namespace NIFS{    

template < typename T>
class list_nbr : public std::list< T* >
{
public:
    list_nbr() {}
	~list_nbr() {}
    T* operator[]( int index );	
	bool has_obj( T* op );
	void remove( T* op );
	void ret( int index, T* op );
};


template < typename T>
T* 
list_nbr::operator[]( int index ) 
{
    if ( empty() ) return NULL; 
    
	int n = index % size();
	iterator iter = begin();
	for ( int i=0; i<n; i++ ) 
 		iter++;
	return (*iter);
}


template < typename T>
bool 
list_nbr::has( T* op )
{
	for ( iterator iter = begin(); iter != end(); iter++) 
		if ( *iter == op)
			return true;
	return false;
}

template < typename T>
void
list_nbr::remove( T* op)
{
    if ( empty() ) return; 
    
    iterator iter;
	for ( iter = begin(); ( iter != end() ) && ( *iter != op ); iter++ );
	if ( iter != end() )
		erase( iter );
}


template < typename T>
void
list_nbr::set( int index, T* op )
{
    if ( empty() ) return;
	 
	int n = index % size();
	iterator iter = begin();
	for ( int i=0; i<n; i++ ) 
 		iter++;
	(*iter) = op;
}



} // namespace NIFS

#endif // CLASS_NBRLIST_H

