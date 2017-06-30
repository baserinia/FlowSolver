// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on:	21 Mar 2005
// Last update:	21 Mar 2005
//------------------------------------------------------------------------------

#include "my_except.h"



std::string 
NIFS::my_except::ui_to_str( unsigned n )
{
	char str[10];
	sprintf( str, "%d", n );
	return static_cast< std::string >( str );
}

NIFS::param_except::param_except( const std::string& msg ) :
					my_except( msg )
{
	
}


const char* 
NIFS::param_except::what() const throw()
{
	std::string new_msg;
	char		line_num_str[5];
	sprintf( line_num_str, "%d", m_line_num );
	new_msg = 	"error in \"" + 
				m_file_name + 
				"\" on line " + 
				line_num_str +
				": " + 
				m_msg;
	return new_msg.c_str();
}


