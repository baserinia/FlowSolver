// Incompressible Flow Solver
// Copyright (C) 2005  Amir R. Baserinia
// Created on:    12 Apr 2004
// Last update:    29 Jun 2017
//-----------------------------------------------------------------------------

#ifndef TEMPLATE_CLASS_LIST_NBR_H
#define TEMPLATE_CLASS_LIST_NBR_H

#include <list>

namespace NIFS{    

//--- INTERFACE ---

template < typename T>
class list_nbr : public std::list< T* >
{
    public:
        list_nbr() {}
        ~list_nbr() {}
        T* operator[]( unsigned index );    
        bool has_obj( T* op );
        void remove( T* op );
        void set( unsigned index, T* op );
};

//--- IMPLEMENTATION ---

template < typename T>
T* 
list_nbr<T>::operator[]( unsigned index ) 
{
    if ( this->size() == 0 ) return 0; 
    
    unsigned n = index % this->size();
    typename list_nbr<T>::iterator iter = this->begin();
    for ( unsigned i=0; i<n; i++ ) 
         ++iter;
    return (*iter);
}

template < typename T>
bool 
list_nbr<T>::has_obj( T* op )
{
    typename list_nbr<T>::iterator iter;
    for ( iter = this->begin(); iter != this->end(); ++iter ) 
        if ( *iter == op)
            return true;
    return false;
}

template < typename T>
void
list_nbr<T>::remove( T* op)
{
    if ( this->size() == 0 ) return; 
    
    typename list_nbr<T>::iterator iter;
    for ( iter = this->begin(); ( iter != this->end() ) && ( *iter != op ); iter++ );
    if ( iter != this->end() )
        this->erase( iter );
}

template < typename T>
void
list_nbr<T>::set( unsigned index, T* op )
{
    if ( this->size() == 0 ) return; 
     
    unsigned n = index % this->size();
    typename list_nbr<T>::iterator iter = this->begin();
    for ( unsigned i=0; i<n; i++ ) 
         iter++;
    (*iter) = op;
}

} // namespace NIFS

#endif //--- TEMPLATE_CLASS_LIST_NBR_H

