// Incompressible Flow lin_sys

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	16 Apr 2004
// Last update:	12 Apr 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_LIN_SYS_H
#define CLASS_NIFS_C_LIN_SYS_H

#include <vector>
#include "elem.h"
#include "matrix_block.h"
#include "vector_block.h"

namespace NIFS {    

// linear system manipulator and lin_sys
class c_lin_sys {
	public:
		c_lin_sys() {}
		~c_lin_sys() {}

	void init_lin_sys( int eqn_num, int entry_num);
	void push_entry( double val, int row, int col );
	void push_rhs( double val, int row );
	double get_x( int i ) { return m_x[i]; }

	void init_lin_sys_block( int eqn_num_block, int entry_num_block, unsigned index_beg, unsigned index_end  );
	void push_entry_block( const matrix_block< double >& mat, int row_blk, int col_blk );
	void push_rhs_block( const vector_block< double >& vec, int row_blk );
	vector_block<double> get_x_block( int i );
	double get_x_block_elem( int i, unsigned j );

	int  get_eqn_num() { return m_eqn_num; }
	int  get_entry_num() { return m_entry_num; }
	void reset();
	void solve_lin_sys();
	void convert_compressed_column( void );
	void print( void );	// test
	void print_coef_mat( std::string cof_file_name );
	void print_rhs_vec( std::string rhs_file_name );
	

	private:
		std::vector<double>	m_b;	// right hand side vector
		std::vector<double>	m_x;	// solution vector
		std::vector<double>	m_A;	// element value array
		std::vector<int>	m_r;	// row array
		std::vector<int>	m_c;	// compressed column array
		std::vector<c_elem>	m_elem;	// full elements vector
		unsigned m_eqn_num;
		unsigned m_entry_num;
		unsigned m_index_beg;
		unsigned m_index_end;
		unsigned m_block_size;
};

} // namespace NBCFD



#endif //--- CLASS_NIFS_C_LIN_SYS_H
