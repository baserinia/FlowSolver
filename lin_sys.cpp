// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	18 Apr 2004
// Last update:	12 Apr 2005
//------------------------------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <iomanip>
#include "lin_sys.h"
#include "umfpack.h" 
#include "my_except.h"
#include "fstream"

void 
NIFS::c_lin_sys::init_lin_sys( int eqn_num, int entry_num )
{
	m_eqn_num = eqn_num;
	m_entry_num = entry_num;
	
	m_elem.reserve( m_entry_num );

	m_A.reserve( m_entry_num );
	m_c.reserve( m_eqn_num+1 ); 
	m_r.reserve( m_entry_num );
	m_b.reserve( m_eqn_num );

	m_x.reserve( m_eqn_num ); 
	m_x.assign ( m_eqn_num , 0.0 );
	
}

void 
NIFS::c_lin_sys::init_lin_sys_block( int eqn_num_block, 
									 int entry_num_block, 
									 unsigned index_beg, 
									 unsigned index_end   )
{
	m_index_beg = index_beg;
	m_index_end = index_end;
	
	m_block_size = m_index_end - m_index_beg;
	
	m_eqn_num = eqn_num_block * m_block_size;
	m_entry_num = entry_num_block * m_block_size * m_block_size;

	m_elem.reserve( m_entry_num );

	m_A.reserve( m_entry_num );
	m_c.reserve( m_eqn_num+1 ); 
	m_r.reserve( m_entry_num );
	m_b.reserve( m_eqn_num );

	m_x.reserve( m_eqn_num ); 
	m_x.assign ( m_eqn_num , 0.0 );
}


void 
NIFS::c_lin_sys::convert_compressed_column( void )
{
	std::sort( m_elem.begin(), m_elem.end() );
	
	unsigned k = 0;
	unsigned nz_entry_num = 0; // nonzero entry num
	
	m_c.push_back( 0 );
	for ( unsigned i = 0 ; i < m_entry_num ; ++i ) {
		if ( m_elem[i].get_val() != 0.0 ) {
			m_A.push_back( m_elem[i].get_val() );
			m_r.push_back( m_elem[i].get_row() );
			if ( m_elem[i].get_col() != k ) {
				k = m_elem[i].get_col();
				m_c.push_back( nz_entry_num ); 
			}
			++nz_entry_num;
		}
	}
	m_c.push_back( nz_entry_num );
	m_elem.clear();
}


void 
NIFS::c_lin_sys::reset()
{
	m_b.clear();
	m_c.clear();
	m_A.clear();
	m_r.clear();
	m_elem.clear();
}


void 
NIFS::c_lin_sys::solve_lin_sys( void )
{
	double* null = ( double* ) NULL;
	void* symbolic;
	void* numeric;
	double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
	Control[UMFPACK_PRL] = 6;
	
	umfpack_di_symbolic( get_eqn_num() , get_eqn_num() , 
						 &m_c[0] , &m_r[0] , &m_A[0] , 
						 &symbolic , null , null);

	umfpack_di_numeric( &m_c[0] , &m_r[0] , &m_A[0] , 
						symbolic , &numeric , Control , Info);
//	umfpack_di_report_info(Control,Info);

	umfpack_di_free_symbolic( &symbolic );
	
	umfpack_di_solve( UMFPACK_A , &m_c[0] , &m_r[0] , &m_A[0] , 
					  &m_x[0] , &m_b[0] , numeric , null , null);
					  
	umfpack_di_free_numeric( &numeric );
}


void 
NIFS::c_lin_sys::push_entry( double val, int row, int col )
{
	c_elem elem;
	elem.set_val( val );
	elem.set_row( row );
	elem.set_col( col );
	m_elem.push_back( elem );	
}

void 
NIFS::c_lin_sys::push_rhs( double val, int row )
{
	m_b[row] = val;
}

void 
NIFS::c_lin_sys::push_entry_block( const matrix_block< double >& mat, int row_blk, int col_blk )
{
	c_elem elem;
	
	for ( unsigned i = 0 ; i < m_block_size ; ++i )
		for ( unsigned j = 0 ; j < m_block_size ; ++j ) {
			elem.set_val( mat( i + m_index_beg )( j + m_index_beg ) );
			elem.set_row( row_blk * m_block_size + i );
			elem.set_col( col_blk * m_block_size + j );
	
			m_elem.push_back( elem );	
		}
}

void 
NIFS::c_lin_sys::push_rhs_block( const vector_block< double >& vec, int row_blk )
{
	for ( unsigned i = 0 ; i < m_block_size ; ++i ) {
		m_b[ row_blk * m_block_size + i ] = vec( i + m_index_beg );
	}
}

NIFS::vector_block<double>
NIFS::c_lin_sys::get_x_block( int i )
{
	vector_block<double> state( 0.0 );
	for ( unsigned j = 0; j < m_block_size ; ++j )
		state[ j + m_index_beg ] = m_x[ i * m_block_size + j ];
	return state;
}

double
NIFS::c_lin_sys::get_x_block_elem( int i, unsigned j )
{
	return m_x[ i * m_block_size + j - m_index_beg];
}

void 
NIFS::c_lin_sys::print( void )
{
//	std::cout << m_elem.size() <<  std::endl;
	for ( unsigned i = 0; i < m_entry_num; ++i ) {
		std::cout  
		<<	m_elem[i].get_val() << ' ' 
		<<	m_elem[i].get_row() << ' ' 
		<<	m_elem[i].get_col() << std::endl;
	}

}

void 
NIFS::c_lin_sys::print_coef_mat( std::string cof_file_name )
{
	std::ofstream ofs;
	ofs.open( cof_file_name.c_str() , std::ios::out );
	if ( !ofs ) 
		throw NIFS::my_except(	"Coef matrix file cannot be created." );

	ofs.setf( std::ios::scientific | std::ios::showpoint  );
	ofs.precision( 10 ); 
	for ( unsigned i = 0; i < m_entry_num; ++i ) {
		ofs <<
		std::setw( 20  ) << m_elem[i].get_val() <<
		std::setw( 8 ) << m_elem[i].get_row() <<
		std::setw( 8 ) << m_elem[i].get_col() << std::endl;
	}
}

void 
NIFS::c_lin_sys::print_rhs_vec( std::string rhs_file_name )
{
	std::ofstream ofs;
	ofs.open( rhs_file_name.c_str() , std::ios::out );
	if ( !ofs ) 
		throw NIFS::my_except(	"RHS vector file cannot be created." );

	ofs.setf( std::ios::scientific | std::ios::showpoint  );
	ofs.precision( 10 ); 
	for ( unsigned i = 0; i < m_eqn_num; ++i ) {
		ofs << std::setw( 20 ) << m_b[i] << std::endl;
	}
	
}
