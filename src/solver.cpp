//	Copyright (C) 2006  Amir R. Baserinia
//
// 	This file is part of SMAIF (Simple Mesh Adaptor for Incompressible Flow)
//
//	SMAIF is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	SMAIF is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with SMAIF; if not, write to the Free Software
//	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
//------------------------------------------------------------------------------
// Created on:	16 Apr 2004
// Last update:	20 Jun 2006
//------------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include "solver.h"
#include "cell.h"
#include "face.h"
#include "vector_block.h"

//--- IMPLEMENTATION ---


void 
NIFS::c_solver::initialize( c_mesh* mesh_ptr )
{
	m_mesh_ptr = mesh_ptr;
	m_mesh_ptr->update_geometry();
  m_mesh_ptr->load_initial_sol( m_info.get_args()->get_init_file_name( ) );
  m_mesh_ptr->update_solution();
}

//---( solve_lin_sys )---------------------------------------
// solve linear system of equations
void 
NIFS::c_solver::solve_lin_sys()
{
	c_cell* cell_ptr;
	c_cell* cell_nbr_ptr;
	unsigned i;
	unsigned index_beg = 0;
	unsigned index_end = 4;
	
	if ( m_info.get_params()->get_solve_puv() == 0.0 )
		index_beg = 3;

	if ( m_info.get_params()->get_solve_t() == 0.0 )
		index_end = 3;
	
	if ( index_beg == index_end )
		return;	

	//--- setting up the system
	m_lin_sys.init_lin_sys_block( m_mesh_ptr->get_cell_num(), 
						 		  m_mesh_ptr->get_coef_num(),
						 		  index_beg,
						 		  index_end );

						 		  
	for ( cell_ptr = m_mesh_ptr->get_first_cell(); 
		  cell_ptr != NULL; 
		  cell_ptr = m_mesh_ptr->get_next_cell() ) {

		m_lin_sys.push_entry_block( cell_ptr->get_coef_slf(), 
									cell_ptr->get_index(),
									cell_ptr->get_index() );

		for ( i = 0; i < cell_ptr->get_cell_num(); ++i ) {
			cell_nbr_ptr = static_cast<c_cell*>( cell_ptr->get_cell( i ) );
			m_lin_sys.push_entry_block( cell_ptr->get_coef_nbr( i ), 
										cell_ptr->get_index(),
										cell_nbr_ptr->get_index() );
		}

		m_lin_sys.push_rhs_block( cell_ptr->get_coef_rhs(), cell_ptr->get_index() );
	}
	m_lin_sys.convert_compressed_column();
	if ( m_info.get_params()->get_soln_rlx() != 0 )
	  m_lin_sys.solve_lin_sys();

	for ( cell_ptr = m_mesh_ptr->get_first_cell(), i = 0; 
		  cell_ptr != NULL; 
		  cell_ptr = m_mesh_ptr->get_next_cell(), ++i ) {
//		if ((i+1)%16!=0)
		  cell_ptr->set_state( m_lin_sys.get_x_block( i ), index_beg, index_end );
//		else
//		  cell_ptr->set_state( cell_ptr->get_state() );
	}
	m_lin_sys.reset();
}


//---( solution_loop )---
// solution iteration loop
void 
NIFS::c_solver::solution_loop()
{
	int i = 0;
	double	crit;
	vector_block< double > res_rms;
	vector_block< double > res_max;

	m_info.get_log_file()->residual_header();

	//--- main solution loop
	do {
		m_mesh_ptr->update_solution();
		solve_lin_sys();
		
		i++;
		res_rms = m_mesh_ptr->get_residual_rms();
		res_max = m_mesh_ptr->get_residual_max();
		
		crit = 0.0;
		for ( int j = 0 ; j < 4 ; ++j ) {
			if ( crit < res_max( j ) )
				crit = res_max( j );
		}
		
		m_info.get_log_file()->print_residual( i , res_rms , res_max );
		
	} while ( ( crit > m_info.get_params()->get_res_trsh() ) && 
			  (   i  < m_info.get_params()->get_max_iter() ) ); 

//	m_mesh_adap_ptr = new c_mesh_adap( *m_mesh_ptr );
//	m_mesh_adap_ptr->update_geometry();	
	
//	m_mesh_adap_ptr->set_bg_mesh( m_mesh_ptr );
//	m_mesh_adap_ptr->adapt();
//	m_mesh_adap_ptr->update_geometry();
	
//	m_mesh_ptr->release_mesh(); 
//	m_mesh_ptr = m_mesh_adap_ptr;
//	m_mesh_ptr->update_geometry();
//	std::cout << "mesh copy complete" << std::endl;

//	m_mesh_adap_ptr->save_mesh();
			  
	m_info.get_log_file()->residual_footer();
}


//---( export_results )---
// export results to the corresponding files
void 
NIFS::c_solver::export_results( void )
{
	c_cell* cell_ptr;


	if ( m_info.get_params()->get_soln_out() != 0.0 ) 
	  m_info.get_res_file()->open_soln_file( m_info.get_args()->get_soln_file_name( ) ); 		

	if ( m_info.get_params()->get_grad_out() != 0.0 ) 
	  m_info.get_res_file()->open_grad_file( m_info.get_args()->get_grad_file_name( ) ); 		

	if ( m_info.get_params()->get_hess_out() != 0.0 ) 
	  m_info.get_res_file()->open_hess_file( m_info.get_args()->get_hess_file_name( ) ); 	

	if ( m_info.get_params()->get_solb_out() != 0.0 ) 
	  m_info.get_res_file()->open_solb_file( m_info.get_args()->get_solb_file_name( ) ); 		

	if ( m_info.get_params()->get_geom_out() != 0.0 ) 
	  m_info.get_res_file()->open_geom_file( m_info.get_args()->get_geom_file_name( ) ); 		

	if ( m_info.get_params()->get_rsdl_out() != 0.0 ) 
	  m_info.get_res_file()->open_rsdl_file( m_info.get_args()->get_rsdl_file_name( ) ); 		
	
	for ( cell_ptr = m_mesh_ptr->get_first_cell(); 
		  cell_ptr != NULL; 
		  cell_ptr = m_mesh_ptr->get_next_cell() ) {
		cell_ptr->print_result();
	}
	m_info.get_log_file()->print_line( "results dumped to file" );
	m_info.get_log_file()->print_elapsed_time();

  if ( m_info.get_params()->do_out_coef_mat() != 0.0 ) {
    m_lin_sys.print_coef_mat( m_info.get_args()->get_cmat_file_name() );
  }
  
  if ( m_info.get_params()->do_out_rhs_vec() != 0.0 ) {
    m_lin_sys.print_rhs_vec( m_info.get_args()->get_rhsv_file_name() );
  }
}


//---( export_params )---
// exports physical parametes to the log file
void 
NIFS::c_solver::export_params( void )
{
	char str[128];
	
	sprintf( str, " order ofaccuracy:      alpha = %8.4f" ,  m_info.get_params()->get_acc_order() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " density:               rho   = %8.4f [kg/m^3]" ,  m_info.get_params()->get_rho() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " viscosity:             mu    = %8.4f [kg/m.s]" ,  m_info.get_params()->get_mu() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " diffusivity:           gamma = %8.4f [W/K.m] " ,  m_info.get_params()->get_gamma() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " heat capacity:         c_p   = %8.4f [J/K.kg]" , m_info.get_params()->get_c_m() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " expansion coef:        beta  = %8.4f [1/K]   " ,  m_info.get_params()->get_beta() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " gravity-x:             g_x   = %8.4f [m/s^2] " ,  m_info.get_params()->get_g_x() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " gravity-y:             g_y   = %8.4f [m/s^2] " ,  m_info.get_params()->get_g_y() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " reference temperature: T_0   = %8.4f [K]     " ,  m_info.get_params()->get_ref_temp() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " pressure scale:        p_ref = %8.4f [Pa]    " ,  m_info.get_params()->get_pres_scale() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " velocity scale:        v_ref = %8.4f [m/s]   " ,  m_info.get_params()->get_vel_scale() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " temperature scale:     T_ref = %8.4f [K]     " ,  m_info.get_params()->get_temp_scale() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " gradient relaxation:   omega = %8.4f " ,  m_info.get_params()->get_grad_rlx() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " convergnece criterion: eps   = %8.4f " ,  m_info.get_params()->get_res_trsh() );
	m_info.get_log_file()->print_line( str );
	sprintf( str, " max iteration:         it    = %8d   " ,  m_info.get_params()->get_max_iter() );
	m_info.get_log_file()->print_line( str );

	m_info.get_log_file()->print_separator();
}
