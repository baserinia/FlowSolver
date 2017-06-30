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
//----------------------------------------------------------------------------
// Created on: 	10 Apr 2004
// Last update: 27 Jun 2006
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_PARAMS_H
#define CLASS_NIFS_C_PARAMS_H

#include <string>
#include "parser.h"

namespace NIFS{ 

// INTERFACE

class c_params {
	public:
		c_params() {}
		virtual ~c_params() {}
		void	load_params_file( std::string file_name);
		
		double	get_solve_puv()  { return m_parser.get_param( "solve_puv" ); }
		double	get_solve_t()    { return m_parser.get_param( "solve_t" ); }
		double	get_rho()        { return m_parser.get_param( "rho" ); }
		double	get_mu()         { return m_parser.get_param( "mu" ); }
		double	get_gamma()      { return m_parser.get_param( "gamma" ); }
		double	get_c_m()        { return m_parser.get_param( "c_m" ); }
		double	get_beta()       { return m_parser.get_param( "beta" ); }
		double	get_g_x()        { return m_parser.get_param( "g_x" ); }
		double	get_g_y()        { return m_parser.get_param( "g_y" ); }
		double	get_ref_temp()   { return m_parser.get_param( "ref_temp" ); }
		double	get_grad_rlx()   { return m_parser.get_param( "grad_rlx" ); }
		double	get_hess_rlx()   { return m_parser.get_param( "hess_rlx" ); }
		double	get_soln_rlx()   { return m_parser.get_param( "soln_rlx" ); }
		double	get_grad_pow()   { return m_parser.get_param( "grad_pow" ); }
		double	get_hess_pow()   { return m_parser.get_param( "hess_pow" ); }
		double	get_res_trsh()   { return m_parser.get_param( "res_trsh" ); }
		double	get_acc_order()  { return m_parser.get_param( "acc_order" ); }
		double	get_grad_swc()   { return m_parser.get_param( "gradient" ); }
		double	get_hess_swc()   { return m_parser.get_param( "hessian" ); }
		double	get_grad_vrt()   { return m_parser.get_param( "grad_vrt" ); }
		double	get_hess_vrt()   { return m_parser.get_param( "hess_vrt" ); }
		double	get_hess_acc()   { return m_parser.get_param( "hess_acc" ); }
		
		double	get_vel_scale()  { return m_parser.get_param( "vel_scale" ); }
		
    bool	  is_grad_out(); 
    double  get_initialize() { return m_parser.get_param( "initialize" ); }
    double  get_bnd_out()    { return m_parser.get_param( "bnd_out" ); } 
    double  get_soln_out()   { return m_parser.get_param( "soln_out" ); } 
    double  get_grad_out()   { return m_parser.get_param( "grad_out" ); } 
    double  get_hess_out()   { return m_parser.get_param( "hess_out" ); } 
    double  get_solb_out()   { return m_parser.get_param( "solb_out" ); }
    double  get_rsdl_out()   { return m_parser.get_param( "rsdl_out" ); }
    double  get_geom_out()   { return m_parser.get_param( "geom_out" ); }
    double  do_out_coef_mat(){ return m_parser.get_param( "coef_mat" ); }
    double  do_out_rhs_vec() { return m_parser.get_param( "rhs_vec" ); }
    double	get_pres_scale();
    double	get_temp_scale();
    double  get_smt_trsh()   { return m_parser.get_param( "smt_trsh" ); }
    double  get_adap_trsh()  { return m_parser.get_param( "adap_trsh" ); }

    int get_max_iter()  { return static_cast<int>( m_parser.get_param( "max_iter" ) ); }
		int get_gh_iter()   { return static_cast<int>( m_parser.get_param( "gh_iter" ) ); }
		int get_smt_iter()  { return static_cast<int>( m_parser.get_param( "smt_iter" ) ); }
		int get_adap_iter() { return static_cast<int>( m_parser.get_param( "adap_iter" ) ); }
				
  private:
    c_parser	m_parser;
    double sqr( double x ) { return x * x; }
};

} // namespace NIFS

#endif // CLASS_NIFS_C_PARAMS_H

