// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//-----------------------------------------------------------------------------
// Created on:	09 Mar 2005
// Last update:	11 Mar 2005
//-----------------------------------------------------------------------------

#ifndef TEMPLATE_CLASS_NIFS_VECTOR_N_H
#define TEMPLATE_CLASS_NIFS_VECTOR_N_H

#include <cmath>
//#include <iostream>
//#include <iomanip>


namespace NIFS {    //--- Nummy Incompressible Flow Solver ---

//--- INTERFACE ---

template < typename T, unsigned n >
class vector_n {
	public:
		vector_n( void );
		vector_n( const T& r );
		vector_n( const vector_n& vec );
		virtual ~vector_n() { }
		void operator=( const vector_n& vec );
		void operator=( const T& r );
		void operator+=( const vector_n& vec );
		void operator-=( const vector_n& vec );
		void operator*=( const T& r );
		void operator/=( const T& r );
		vector_n operator-() const;
		vector_n operator+( const vector_n& vec ) const;
		vector_n operator-( const vector_n& vec ) const;
		vector_n operator*( const T& r ) const;
		vector_n operator/( const T& r ) const;
		T operator*( const vector_n& vec ) const;
		T& operator[]( unsigned i );
		T operator()( unsigned i ) const;
		T get_magnit_sqr() const;
		T get_magnit() const;
		vector_n get_unit() const;
		vector_n mult_elem( const vector_n& vec ); // element-wise mutiplication
	private:
		T	v[ n ];	//--- n-dimensional vector ---
};

//--- IMPLEMENTATION ---

template < typename T, unsigned n >
vector_n<T,n>::vector_n( void )
{
	for ( unsigned i = 0; i < n; ++i )
		v[i] = static_cast<T>(0.0);
}  


template < typename T, unsigned n >
vector_n<T,n>::vector_n( const T& r )
{
	for ( unsigned i = 0; i < n; ++i )
		v[i] = r;
}  


template < typename T, unsigned n >
vector_n<T,n>::vector_n( const vector_n& vec ) 
{
    *this = vec;
}  


template < typename T, unsigned n >
void 
vector_n<T,n>::operator=( const vector_n<T,n>& vec ) 
{
	for ( unsigned i = 0; i < n; ++i )
		v[i] = vec.v[i];
}  


template < typename T, unsigned n >
void 
vector_n<T,n>::operator=( const T& r )
{
	for ( unsigned i = 0 ; i < n; ++i )
		v[i] = r;
}


template < typename T, unsigned n >
void 
vector_n<T,n>::operator+=( const vector_n<T,n>& vec ) 
{
	for ( unsigned i = 0; i < n; ++i )
		v[i] += vec.v[i];  
}  


template < typename T, unsigned n >
void 
vector_n<T,n>::operator-=( const vector_n<T,n>& vec ) 
{
	for ( unsigned i = 0; i < n; ++i )
		v[i] -= vec.v[i];  
}  


template < typename T, unsigned n >
void 
vector_n<T,n>::operator*=( const T& r ) 
{
	for ( unsigned i = 0; i < n; ++i )
		v[i] *= r;  
}  


template < typename T, unsigned n >
void 
vector_n<T,n>::operator/=( const T& r ) 
{
	for ( unsigned i = 0; i < n; ++i )
		v[i] /= r;  
}  


template < typename T, unsigned n >
vector_n<T,n> 
vector_n<T,n>::operator-( ) const
{
	vector_n<T,n> vec_tmp;
	for ( unsigned i = 0; i < n; ++i )
		vec_tmp.v[i] -= v[i];  
    return vec_tmp;
}  


template < typename T, unsigned n >
vector_n<T,n> 
vector_n<T,n>::operator+( const vector_n& vec ) const
{
	vector_n<T,n> vec_tmp;
	for ( unsigned i = 0; i < n; ++i )
		vec_tmp[i] = (*this)(i) + vec(i);  
    return vec_tmp;
}  


template < typename T, unsigned n >
vector_n<T,n> 
vector_n<T,n>::operator-( const vector_n<T,n>& vec ) const
{
	vector_n<T,n> vec_tmp;
	for ( unsigned i = 0; i < n; ++i )
		vec_tmp.v[i] = v[i] - vec.v[i];  
    return vec_tmp;

}  


template < typename T, unsigned n >
vector_n<T,n> 
vector_n<T,n>::operator*( const T& r ) const
{
	vector_n<T,n> vec_tmp;
	for ( unsigned i = 0; i < n; ++i )
		vec_tmp.v[i] = v[i] * r;  
    return vec_tmp;
}  


template < typename T, unsigned n >
vector_n<T,n> 
vector_n<T,n>::operator/( const T& r ) const
{
	vector_n<T,n> vec_tmp;
	for ( unsigned i = 0; i < n; ++i )
		vec_tmp.v[i] = v[i] / r;  
    return vec_tmp;
}


template < typename T, unsigned n >
T
vector_n<T,n>::operator*( const vector_n<T,n>& vec ) const
{
	T r = static_cast<T>( 0.0 );
	for ( unsigned i = 0; i < n; ++i )
		r += v[i] * vec.v[i];  
    return r;
}  


template < typename T, unsigned n >
T&
vector_n<T,n>::operator[]( unsigned i )
{
    return v[i];
}  


template < typename T, unsigned n >
T
vector_n<T,n>::operator()( unsigned i ) const
{
    return v[i];
}  


template < typename T, unsigned n >
T
vector_n<T,n>::get_magnit_sqr() const
{
	T r = static_cast<T>( 0.0 );
	for ( unsigned i = 0; i < n; ++i )
		r += v[i] * v[i];  
    return r;
}  


template < typename T, unsigned n >
T
vector_n<T,n>::get_magnit() const
{
	T r = static_cast<T>( 0.0 );
	for ( unsigned i = 0; i < n; ++i )
		r += v[i] * v[i];  
    return std::sqrt( r );
}  


template < typename T, unsigned n >
vector_n<T,n>
vector_n<T,n>::get_unit() const
{
    return (*this) / get_magnit( ); 
}  


template< typename T , unsigned n >
vector_n< T , n > 
vector_n< T , n  >::mult_elem( const vector_n<T,n>& vec )
{
	vector_n< T , n > vec_temp;
	for ( unsigned i = 0 ; i < n ; ++i ) {
		vec_temp[i] = (*this)(i) * vec(i);
	}
	return vec_temp;
}


/*
template < typename T, unsigned n >
std::ostream& 
operator<<( std::ostream& os, const vector_n<T,n>& vec )
{
	for (unsigned i=0; i< n; ++i )
		os << std::setw(10) << vec(i);
	os << std::endl;
	return os;
}
*/

} //--- namespace NIFS ---

#endif // TEMPLATE_CLASS_NIFS_VECTOR_N_H
