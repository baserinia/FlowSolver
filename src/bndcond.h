// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	16 Mar 2005
// Last update:	16 Mar 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_BNDCOND_H
#define CLASS_NIFS_BNDCOND_H

#include <string>
#include <vector>
#include "parser.h"
#include "vector_block.h"

namespace NIFS{ // Naamey Incompressible Flow Solver

// INTERFACE

class Cbndcond {
	public:
		Cbndcond() {}
		virtual ~Cbndcond() {}
		void load_bc_file( std::string file_name );
		unsigned get_bnd_num( void );
		double get_pres( unsigned bnd_code, double x, double y );
		double get_velu( unsigned bnd_code, double x, double y );
		double get_velv( unsigned bnd_code, double x, double y );
		double get_temp( unsigned bnd_code, double x, double y );
		unsigned get_bnd_type( unsigned bnd_code );
	private:
		unsigned extract_index( const std::string& str, 
								const std::string& sub_str );
		void set_bc_string( const std::string& bc_string );
		
	private:
		Cparser m_parser;
		unsigned m_bnd_num;
		std::vector< unsigned > m_bc_type;
		std::vector< std::string > m_bc_pres_exp;
		std::vector< std::string > m_bc_velu_exp;
		std::vector< std::string > m_bc_velv_exp;
		std::vector< std::string > m_bc_temp_exp;
};

} // namespace NIFS

#endif //--- CLASS_NIFS_BNDCOND_H ---
