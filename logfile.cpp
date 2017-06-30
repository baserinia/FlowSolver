// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:	16 Mar 2005
// Last update:	16 Mar 2005
//-----------------------------------------------------------------------------

#include "logfile.h"
#include "my_except.h"


NIFS::Clogfile::Clogfile( std::string logfile_name )
{
	open_logfile( logfile_name );
}

void 
NIFS::Clogfile::open_logfile( std::string logfile_name )
{
	m_ofs.open( logfile_name.c_str(), std::ios::out );	
	if ( !m_ofs ) {
		throw NIFS::my_except(	"logfile \"" + 
								logfile_name + 
								"\" cannot be created." );
	}
	print_line( "--- NIFS LOG FILE ---"  );
	print_line( "--- COPYRIGHT (C) A. R. BASERINIA ---" );
	m_ofs << "--- RUN DATE: ";	
	print_cal_time( );
	print_separator();
	m_start_time = time( NULL );
}


void 
NIFS::Clogfile::close_logfile()
{
	m_ofs.close();
}

void 
NIFS::Clogfile::print_str( const char* str )
{
	std::string std_str( str );
	m_ofs << std_str;
}	


void 
NIFS::Clogfile::print_line( const char* str )
{
	std::string std_str( str );
	m_ofs << std_str << std::endl;
}	


void 
NIFS::Clogfile::print_cal_time( void )
{
	std::time_t the_time;
	the_time = std::time( NULL );
	m_ofs << std::ctime( &the_time );
}

void 
NIFS::Clogfile::print_elapsed_time( void )
{
	std::time_t the_time = std::time( NULL );
	std::time_t elapsed_time = the_time - m_start_time;
	m_ofs << "elapsed time: " << elapsed_time << std::endl;
}

void 
NIFS::Clogfile::print_separator( void )
{
	m_ofs 	<< std::endl 
			<< "------------------------------------------------------------" 
			<< std::endl
			<< std::endl;	
}