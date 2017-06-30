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
// Created on:	28 Feb 2005
// Last update:	16 Mar 2005
//-----------------------------------------------------------------------------

#ifndef CLASS_NIFS_C_INFO_H
#define CLASS_NIFS_C_INFO_H

#include <string>
#include "params.h"
#include "log_file.h"
#include "res_file.h"
#include "bnd_cond.h"
#include "args.h"


namespace NIFS{ 

//--- INTERFACE ---

class c_info {
  public:
    c_info();
    virtual ~c_info();

    c_params* get_params( void ) const; 
    c_log_file* get_log_file( void ) const;
    c_res_file* get_res_file( void ) const;
    c_bnd_cond* get_bnd_cond( void ) const;
    c_args*     get_args( void ) const;

    void release_params( void );
    void release_log_file( void );
    void release_res_file( void );
    void release_bnd_cond( void );
    void release_args( void );
		
    void release_all( void );
		
  private:
    static c_params*   m_params_ptr;
    static c_log_file* m_log_file_ptr;
    static c_res_file* m_res_file_ptr;
    static c_bnd_cond* m_bnd_cond_ptr;
    static c_args*     m_args_ptr;
};

} // namespace NIFS

#endif // CLASS_NIFS_C_INFO_H
