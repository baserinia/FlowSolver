// Incompressible Flow Solver

// Copyright (C) 2005  Amir R. Baserinia

//----------------------------------------------------------------------------
// Created on:	16 Mar Feb 2005
// Last update:	16 Mar 2005
//-----------------------------------------------------------------------------

#include <string>
#include "info.h"

//-- IMPLEMENTATION ---

// initialize static variables
NIFS::c_params*   NIFS::c_info::m_params_ptr   = NULL;
NIFS::c_log_file* NIFS::c_info::m_log_file_ptr = NULL;
NIFS::c_res_file* NIFS::c_info::m_res_file_ptr = NULL;
NIFS::c_bnd_cond* NIFS::c_info::m_bnd_cond_ptr = NULL;
NIFS::c_args*     NIFS::c_info::m_args_ptr     = NULL;

NIFS::c_info::c_info()
{

}


//---( ~c_info )---
// default constructor
NIFS::c_info::~c_info()
{
//  release_log_file( );
//  release_params( );
//  release_bnd_cond( );
//  release_res_file( );
//  release_args( );
}


NIFS::c_params* 
NIFS::c_info::get_params( void ) const
{
	if ( c_info::m_params_ptr == NULL) 
		return c_info::m_params_ptr = new c_params();
	else
		return c_info::m_params_ptr;
}

NIFS::c_log_file* 
NIFS::c_info::get_log_file( void ) const
{
	if ( c_info::m_log_file_ptr == NULL) 
		return ( c_info::m_log_file_ptr = new c_log_file() );
	else
		return c_info::m_log_file_ptr;
}


NIFS::c_res_file* 
NIFS::c_info::get_res_file( void ) const
{
	if ( c_info::m_res_file_ptr == NULL) 
		return ( c_info::m_res_file_ptr = new c_res_file() );
	else
		return c_info::m_res_file_ptr;
}

NIFS::c_bnd_cond* 
NIFS::c_info::get_bnd_cond( void ) const
{
	if ( c_info::m_bnd_cond_ptr == NULL) 
		return c_info::m_bnd_cond_ptr = new c_bnd_cond();
	else
		return c_info::m_bnd_cond_ptr;
}

NIFS::c_args* 
NIFS::c_info::get_args( void ) const
{
	if ( c_info::m_args_ptr == NULL) 
		return c_info::m_args_ptr = new c_args();
	else
		return c_info::m_args_ptr;
}

void 
NIFS::c_info::release_params( void )
{
	delete c_info::m_params_ptr;	
	c_info::m_params_ptr = NULL;
}

void 
NIFS::c_info::release_log_file( void )
{
	if ( c_info::m_log_file_ptr != NULL ) {
		m_log_file_ptr->close();
		delete c_info::m_log_file_ptr;	
		c_info::m_log_file_ptr = NULL;
	}
}

void 
NIFS::c_info::release_res_file( void )
{
	if ( c_info::m_res_file_ptr != NULL ) {
		m_res_file_ptr->close_all();
		delete c_info::m_res_file_ptr;	
		c_info::m_res_file_ptr = NULL;
	}
}

void 
NIFS::c_info::release_bnd_cond( void )
{
	delete c_info::m_bnd_cond_ptr;	
	c_info::m_bnd_cond_ptr = NULL;
}


void 
NIFS::c_info::release_all( void )
{
	release_log_file( );
	release_params( );
	release_bnd_cond( );
	release_res_file( );
	release_args( );
}

void 
NIFS::c_info::release_args( void )
{
	delete c_info::m_args_ptr;	
	c_info::m_args_ptr = NULL;
}
