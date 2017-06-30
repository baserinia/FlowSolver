// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on: 	09 Mar 2005
// Last update: 11 Mar 2005
//-----------------------------------------------------------------------------

#ifndef TEMPLATE_CLASS_NIFS_MATRIX_N_H
#define TEMPLATE_CLASS_NIFS_MATRIX_N_H

#include <cmath>
#include "vector_n.h"

namespace NIFS{

//--- INTERFACE ---

template < typename T, unsigned n >
class matrix_n {
	public:
		matrix_n( void );
		matrix_n( const T& r );
		matrix_n( const matrix_n& mat );
		virtual ~matrix_n( void ) { }
		void operator=( const matrix_n& mat );
		void operator=( const T& r );
		void operator+=( const matrix_n& mat ); 
		void operator-=( const matrix_n& mat ); 
		void operator*=( T r ); 
		void operator/=( T r ); 
		matrix_n operator+(const matrix_n& mat) const; 
		matrix_n operator-( void ) const;
		matrix_n operator-(const matrix_n& mat) const; 	
		matrix_n operator*( T r ) const; 
		matrix_n operator*( const matrix_n& mat ) const;
		matrix_n operator/( T r ) const; 	
		vector_n<T,n> operator*( const vector_n<T,n>& vec ) const;
		vector_n<T,n>& operator[]( unsigned i );
		vector_n<T,n> operator()( unsigned i ) const;	
		void set( unsigned i, unsigned j, T r );
		T get( unsigned i, unsigned j );
		matrix_n transpose( void ) const;
	private:
		vector_n<T,n>	t[ n ]; //--- n*n square matrix---
};

//--- IMPLEMENTATION ---

template < typename T, unsigned n >
matrix_n<T,n>::matrix_n( void )
{
	for ( unsigned i = 0; i < n; ++i )
		t[i] = static_cast<T>(0.0);
}  


template < typename T, unsigned n >
matrix_n<T,n>::matrix_n( const T& r )
{
	for ( unsigned i = 0; i < n; ++i )
		t[i] = r;
}  


template < typename T, unsigned n >
matrix_n<T,n>::matrix_n( const matrix_n<T,n>& mat ) 
{
    *this = mat;
}  


template < typename T, unsigned n >
void 
matrix_n<T,n>::operator=( const matrix_n<T,n>& mat ) 
{
	for ( unsigned i = 0; i < n; ++i )
    	t[i] = mat.t[i];
}  

template < typename T, unsigned n >
void 
matrix_n<T,n>::operator=( const T& r )
{
	for ( unsigned i = 0; i < n; ++i )
		t[i] = r;
}  

template < typename T, unsigned n >
void 
matrix_n<T,n>::operator+=( const matrix_n<T,n>& mat ) 
{
	for ( unsigned i = 0; i < n; ++i )
		t[i] += mat.t[i];  
}  


template < typename T, unsigned n >
void 
matrix_n<T,n>::operator-=( const matrix_n<T,n>& mat ) 
{
	for ( unsigned i = 0; i < n; ++i )
		t[i] -= mat.t[i];  
}  


template < typename T, unsigned n >
void 
matrix_n<T,n>::operator*=( T r ) 
{
	for ( unsigned i = 0; i < n; ++i )
		t[i] *= r;  
}  


template < typename T, unsigned n >
void 
matrix_n<T,n>::operator/=( T r ) 
{
	for ( unsigned i = 0; i < n; ++i )
		t[i] /= r;  
}  


template < typename T, unsigned n >
matrix_n<T,n>
matrix_n<T,n>::operator+( const matrix_n<T,n>& mat ) const
{
	matrix_n<T,n> mat_tmp;
	for ( unsigned i = 0; i < n; ++i )
		mat_tmp.t[i] = t[i] + mat.t[i];  
    return mat_tmp;
}  


template < typename T, unsigned n >
matrix_n<T,n>
matrix_n<T,n>::operator-( void ) const
{
	matrix_n<T,n> mat_tmp;
	for ( unsigned i = 0; i < n; ++i )
		mat_tmp.t[i] -= t[i];  
    return mat_tmp;
}  


template < typename T, unsigned n >
matrix_n<T,n>
matrix_n<T,n>::operator-( const matrix_n& mat ) const 	
{
	matrix_n<T,n> mat_tmp;
	for ( unsigned i = 0; i < n; ++i )
		mat_tmp.t[i] = t[i] - mat.t[i];  
    return mat_tmp;
}  


template < typename T, unsigned n >
matrix_n<T,n>
matrix_n<T,n>::operator*( T r ) const
{
	matrix_n<T,n> mat_tmp;
	for ( unsigned i = 0; i < n; ++i )
		mat_tmp.t[i] = t[i] * r;  
    return mat_tmp;
}  


template < typename T, unsigned n >
matrix_n<T,n>
matrix_n<T,n>::operator*( const matrix_n<T,n>& mat ) const
{
	matrix_n<T,n> mat_trn = mat.transpose();
	matrix_n<T,n> mat_tmp;
	for ( unsigned i = 0; i < n; ++i )
		for ( unsigned j = 0; j < n; ++j )
			mat_tmp[i][j] = t[i] * mat_trn.t[j];
	return mat_tmp;
}  


template < typename T, unsigned n >
vector_n<T,n>
matrix_n<T,n>::operator*( const vector_n<T,n>& vec ) const
{
	vector_n<T,n> vec_tmp;
	for ( unsigned i = 0; i < n; ++i )
		vec_tmp[i] = vec * t[i];
	return vec_tmp;
}  


template < typename T, unsigned n >
matrix_n<T,n>
matrix_n<T,n>::operator/( T r ) const
{
	matrix_n<T,n> mat_tmp;
	for ( unsigned i = 0; i < n; ++i )
		mat_tmp[i] = t[i] / r;
	return mat_tmp;
}  


template < typename T, unsigned n >
vector_n<T,n>&
matrix_n<T,n>::operator[]( unsigned i )
{
	return t[i]; //--- CAUTION: no range checking
}  


template < typename T, unsigned n >
vector_n<T,n>
matrix_n<T,n>::operator()( unsigned i ) const
{
    return t[i];
}  


template < typename T, unsigned n >
void 
matrix_n<T,n>::set( unsigned i, unsigned j, T r )
{
	t[i][j] = r; //--- CAUTION: no range checking
}  


template < typename T, unsigned n >
T
matrix_n<T,n>::get( unsigned i, unsigned j )
{
	return t[i](j);
}  


template < typename T, unsigned n >
matrix_n<T,n>
matrix_n<T,n>::transpose( void ) const
{
	matrix_n<T,n> mat_tmp;
	for ( unsigned i = 0; i < n; ++i )
		for ( unsigned j = 0; j < n; ++j )
			mat_tmp.t[i][j] = t[j](i);
	return mat_tmp;
}

/*
template < typename T, unsigned n >
std::ostream& 
operator<<( std::ostream& os,  matrix_n<T,n>& mat )
{
	for ( unsigned i = 0; i < n; ++i )
	    os << i << ": " << mat[i];
    return os;
}
*/

} // namespace NIFS

#endif // TEMPLATE_CLASS_NIFS_MATRIX_N_H
