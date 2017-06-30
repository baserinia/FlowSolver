// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	16 Mar 2005
// Last update:	16 Mar 2005
//------------------------------------------------------------------------------

#ifndef CLASS_NIFS_LOGFILE_H
#define CLASS_NIFS_LOGFILE_H

#include <string>
#include <ctime>
#include "fstream"

namespace NIFS{ //---  

class Clogfile
{
	public:
		Clogfile() {}
		Clogfile( std::string logfile_name );
		virtual ~Clogfile() {}
		void open_logfile( std::string logfile_name );
		void close_logfile();
		void print_str( const char* str );
		void print_line( const char* str );
		void print_cal_time( void );
		void print_elapsed_time( void );
		void print_separator( void );
	private:
		std::ofstream m_ofs;
		std::time_t	m_start_time;
};

} //--- namespace NIFS ---

#endif // CLASS_NIFS_LOGFILE_H
