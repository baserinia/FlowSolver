// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//------------------------------------------------------------------------------
// Created on: 	11 Apr  2004
// Last update: 13 Mar 2005
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include "my_except.h"
#include "params.h"
//#include "info.h"


//--- IMPLEMENTATION ---

void 
NIFS::c_params::load_params_file( std::string file_name)
{
	std::string line_str;
	unsigned line_num = 0;
    std::string::size_type index;
	
    std::ifstream ifs( file_name.c_str() , std::ios::in );	
   	if ( !ifs ) {
		throw NIFS::my_except( 
			"unable to open parameters file \"" + file_name + "\"" );
	}
	
	try {	
		std::getline( ifs, line_str );
		while ( !ifs.eof( ) ) {
			line_num++;
			index = line_str.find( '%' );
			if ( index != std::string::npos )	
				line_str = line_str.erase ( index );
			m_parser.eval_postfix( line_str );
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

double	
NIFS::c_params::get_pres_scale()
{
	return 0.5 * get_rho() * sqr( get_vel_scale() );
}

double	
NIFS::c_params::get_temp_scale()
{
	return sqr( get_vel_scale() ) / get_c_m();
}

bool
NIFS::c_params::is_grad_out() 
{
	int out_flags = static_cast< int >( m_parser.get_param( "out_flags" ) );
	return ( out_flags & 1 );
}
/*
void 
NIFS::c_params::export_params( void )
{
	char str[128];
	c_info info;
	
	sprintf( str, " order ofaccuracy:      alpha = %8.4f" ,  get_acc_order() );
	info.get_log_file()->print_line( str );
	sprintf( str, " density:               rho   = %8.4f [kg/m^3]" ,  get_rho() );
	info.get_log_file()->print_line( str );
	sprintf( str, " viscosity:             mu    = %8.4f [kg/m.s]" ,  get_mu() );
	info.get_log_file()->print_line( str );
	sprintf( str, " diffusivity:           gamma = %8.4f [W/K.m] " ,  get_gamma() );
	info.get_log_file()->print_line( str );
	sprintf( str, " heat capacity:         c_p   = %8.4f [J/K.kg]" ,  get_c_m() );
	info.get_log_file()->print_line( str );
	sprintf( str, " expansion coef:        beta  = %8.4f [1/K]   " ,  get_beta() );
	info.get_log_file()->print_line( str );
	sprintf( str, " gravity-x:             g_x   = %8.4f [m/s^2] " ,  get_g_x() );
	info.get_log_file()->print_line( str );
	sprintf( str, " gravity-y:             g_y   = %8.4f [m/s^2] " ,  get_g_y() );
	info.get_log_file()->print_line( str );
	sprintf( str, " reference temperature: T_0   = %8.4f [K]     " ,  get_ref_temp() );
	info.get_log_file()->print_line( str );
	sprintf( str, " pressure scale:        p_ref = %8.4f [Pa]    " ,  get_pres_scale() );
	info.get_log_file()->print_line( str );
	sprintf( str, " velocity scale:        v_ref = %8.4f [m/s]   " ,  get_vel_scale() );
	info.get_log_file()->print_line( str );
	sprintf( str, " temperature scale:     T_ref = %8.4f [K]     " ,  get_temp_scale() );
	info.get_log_file()->print_line( str );
	sprintf( str, " gradient relaxation:   omega = %8.4f " ,  get_grad_rlx() );
	info.get_log_file()->print_line( str );
	sprintf( str, " convergnece criterion: eps   = %8.4f " ,  get_res_trsh() );
	info.get_log_file()->print_line( str );
	sprintf( str, " max iteration:         it    = %8d   " ,  get_max_iter() );
	info.get_log_file()->print_line( str );
}
*/
