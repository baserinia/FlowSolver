// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on: 	16 Mar  2004
// Last update: 24 Mar 2005
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include "bndcond.h"
#include "my_except.h"

//--- IMPLEMENTATION ---

void 
NIFS::c_bndcond::load_bc_file( std::string file_name )
{
	std::string line_str;
	unsigned line_num = 0;
    std::string::size_type index;
	
    std::ifstream ifs( file_name.c_str() , std::ios::in );	
   	if ( !ifs ) {
		throw NIFS::my_except( 
			"unable to open boundary conditions file \"" + file_name + "\"" );
	}
	try {
		m_parser.set_param( "x", 0.0 );
		m_parser.set_param( "y", 0.0 );
		m_parser.set_param( "bnd_num", 0 );
		std::getline( ifs, line_str );
		while ( !ifs.eof( ) ) {
			line_num++;
			index = line_str.find( '%' );
			if ( index != std::string::npos )	
				line_str = line_str.erase ( index ); // erase all after '%'
			m_parser.eval_postfix( line_str );
			set_bc_string( line_str );	// might raise exception
			std::getline( ifs, line_str );
		}
	}
	catch ( NIFS::my_except& e ) {
		e.set_msg( 	"error in \"" + 
					file_name +
					"\" on line " +
					e.ui_to_str( line_num ) +
					": " +
					e.what( ) );
		throw;
	}
	return;
}


unsigned 
NIFS::c_bndcond::get_bnd_num( void )
{ 
	return static_cast<unsigned>( m_parser.get_param( "bnd_num" ) ); 
}


void 
NIFS::c_bndcond::set_bc_string( const std::string& bc_string )
{
	std::string::size_type index;
	std::string first_token;
	unsigned n = 0;
	
	first_token = m_parser.get_token( bc_string, 0 );
	
	if ( first_token == "" )
		return;
	else if ( ( index = first_token.find( "bnd_num" ) ) != std::string::npos ) {
		m_bnd_num = static_cast<unsigned>( m_parser.get_param( first_token ) );
		m_bc_type.reserve( m_bnd_num + 1 );
		m_bc_pres_exp.resize( m_bnd_num + 1 );
		m_bc_velu_exp.resize( m_bnd_num + 1 );
		m_bc_velv_exp.resize( m_bnd_num + 1 );
		m_bc_temp_exp.resize( m_bnd_num + 1 );
	}
	else if ( ( index = first_token.find( "bnd_type" ) ) != std::string::npos ) {
		n = extract_index( first_token, "bnd_type" );
		m_bc_type[n] = static_cast<unsigned>( m_parser.get_param( first_token ) );
	}
	else if ( ( index = first_token.find( "pres" ) ) != std::string::npos ) {
		n = extract_index( first_token, "pres" );
		m_bc_pres_exp[n] = bc_string ;
	}
	else if ( ( index = first_token.find( "velu" ) ) != std::string::npos ) {
		n = extract_index( first_token, "velu" );
		m_bc_velu_exp[0] = bc_string;
	}
	else if ( ( index = first_token.find( "velv" ) ) != std::string::npos ) {
		n = extract_index( first_token, "velv" );
		m_bc_velv_exp[n] = bc_string;
	}
	else if ( ( index = first_token.find( "temp" ) ) != std::string::npos ) {
		n = extract_index( first_token, "temp" );
		m_bc_temp_exp[n] = bc_string;
	} 
	else {
		throw NIFS::my_except( 
			"unidentified boundary condition parameter " + first_token );
	}
	return;
}


double
NIFS::c_bndcond::get_pres( unsigned bnd_code, double x, double y )
{
	m_parser.set_param( "x", x );
	m_parser.set_param( "y", y );
	m_parser.eval_postfix( m_bc_pres_exp[ bnd_code ] );
	return m_parser.get_answer( );
}


double
NIFS::c_bndcond::get_velu( unsigned bnd_code, double x, double y )
{
	m_parser.set_param( "x", x );
	m_parser.set_param( "y", y );
	m_parser.eval_postfix( m_bc_velu_exp[ bnd_code ] );
	return m_parser.get_answer( );	
}


double
NIFS::c_bndcond::get_velv( unsigned bnd_code, double x, double y )
{
	m_parser.set_param( "x", x );
	m_parser.set_param( "y", y );
	m_parser.eval_postfix( m_bc_velv_exp[ bnd_code ] );
	return m_parser.get_answer( );
}


double
NIFS::c_bndcond::get_temp( unsigned bnd_code, double x, double y )
{
	m_parser.set_param( "x", x );
	m_parser.set_param( "y", y );
	m_parser.eval_postfix( m_bc_temp_exp[ bnd_code ] );
	return m_parser.get_answer( );
}


unsigned 
NIFS::c_bndcond::extract_index(	const std::string& str, 
								const std::string& sub_str )
{
	std::string temp_str = str;
	std::string temp_str2;
	std::string::size_type index;
	unsigned n = 0;
	
	index = str.find( sub_str );
	if (  index != std::string::npos )
	{
		temp_str2 = temp_str.erase( index, sub_str.length() );
	}
	
	if ( temp_str2.find_first_not_of( "0123456789" ) != std::string::npos )
		throw NIFS::my_except( 
			"boundary condition identifier \"" + str + "\" is invalid" );
	else  {
		n = atoi( temp_str2.c_str() );
		if ( ( n > get_bnd_num() ) || ( n == 0 ) )
			throw NIFS::my_except( 
				"boundary number " + temp_str + " does not exist" );
	}
	return n;	
}

unsigned 
NIFS::c_bndcond::get_bc_type( unsigned bnd_code )
	m_parser.eval_postfix( m_bc_type[ bnd_code ] );
	return static_cast<unsigned>( m_parser.get_answer( ) );
}