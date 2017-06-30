// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	16 Mar 2005
// Last update:	16 Mar 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_LOG_FILE_H
#define CLASS_NIFS_C_LOG_FILE_H

#include <string>
#include <ctime>
#include "fstream"
#include "vector_block.h"

namespace NIFS{ 

//--- INTERFACE ---

class c_log_file
{
	public:
		c_log_file() {}
		c_log_file( std::string log_file_name );
		virtual ~c_log_file() {}
		void open( std::string log_file_name );
		void close();
		void print_str( const char* str );
		void print_line( const char* str );
		void print_cal_time( void );
		void print_elapsed_time( void );
		void print_separator( void );
		void blank_line( );
		void residual_header( void );
		void print_residual( int, vector_block<double>, vector_block<double> );
		void residual_footer( void );
		std::time_t get_elapsed_time( void );
		
	private:
		std::ofstream	m_ofs;
		std::time_t		m_start_time;
};

} //--- namespace NIFS ---

#endif // CLASS_NIFS_C_LOG_FILE_H
