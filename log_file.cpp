// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:	16 Mar 2005
// Last update:	19 Apr 2005
//-----------------------------------------------------------------------------

#include <iomanip>
#include "log_file.h"
#include "my_except.h"



NIFS::c_log_file::c_log_file( std::string log_file_name )
{
	open( log_file_name );
}

void 
NIFS::c_log_file::open( std::string log_file_name )
{
	m_ofs.open( log_file_name.c_str(), std::ios::out );	
	if ( !m_ofs ) {
		throw NIFS::my_except(	"log_file \"" + 
								log_file_name + 
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
NIFS::c_log_file::close()
{
	m_ofs.close();
}

void 
NIFS::c_log_file::print_str( const char* str )
{
	std::string std_str( str );
	m_ofs << std_str;
}	


void 
NIFS::c_log_file::print_line( const char* str )
{
	std::string std_str( str );
	m_ofs << std_str << std::endl;
}	


void 
NIFS::c_log_file::print_cal_time( void )
{
	std::time_t the_time;
	the_time = std::time( NULL );
	m_ofs << std::ctime( &the_time );
}

void 
NIFS::c_log_file::print_elapsed_time( void )
{
	m_ofs << "elapsed time: " << get_elapsed_time() << std::endl;
}

void 
NIFS::c_log_file::print_separator( void )
{
	blank_line();
	m_ofs << std::setw(80) << std::setfill('-') << "";
	blank_line();
	blank_line();
	m_ofs.fill( ' ' );
}

void 
NIFS::c_log_file::blank_line( )
{
	m_ofs << std::endl;
}



void 
NIFS::c_log_file::residual_header( void )
{
	blank_line();
	print_line( "--- Convergence History ---" );
	blank_line();
	m_ofs 	
		<< "-----+-------+----" 
		<< "------+----------+----------+----------+----" 
		<< "------+----------+----------+----------+"
		<< std::endl;
	m_ofs 	
		<< " itr |  time |    " 
		<< "p_rms |    p_max |    u_rms |    u_max |    " 
		<< "v_rms |    v_max |    T_rms |    T_max |"
		<< std::endl;
	m_ofs 	
		<< "-----+-------+----" 
		<< "------+----------+----------+----------+----" 
		<< "------+----------+----------+----------+"
		<< std::endl;
}

void 
NIFS::c_log_file::print_residual( int itr, 
								  vector_block<double> res_rms, 
								  vector_block<double> res_max)
{
	m_ofs << std::setw( 4 ) << itr << " | ";
	m_ofs << std::setw( 5 ) << get_elapsed_time() << " | ";
	
	m_ofs.setf( std::ios::scientific | std::ios::showpoint  );
	m_ofs.precision( 2 ); 
	m_ofs <<  res_rms.get_pres() << " | ";
	m_ofs <<  res_max.get_pres() << " | ";
	m_ofs <<  res_rms.get_velu() << " | ";
	m_ofs <<  res_max.get_velu() << " | ";
	m_ofs <<  res_rms.get_velv() << " | ";
	m_ofs <<  res_max.get_velv() << " | ";
	m_ofs <<  res_rms.get_temp() << " | ";
	m_ofs <<  res_max.get_temp() << " | ";
	m_ofs << std::endl;
}

void 
NIFS::c_log_file::residual_footer( void )
{
	m_ofs 	
		<< "-----+-------+----"
		<< "------+----------+----------+----------+----"
		<< "------+----------+----------+----------+"
		<< std::endl;
	blank_line();
}


std::time_t 
NIFS::c_log_file::get_elapsed_time( void )
{
	std::time_t the_time = std::time( NULL );
	return ( the_time - m_start_time );
}


