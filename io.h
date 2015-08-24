//
//  io.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__io__
#define __TetrahedralParticlesInConfinement__io__

#include <stdio.h>
#include "struct_def.h"
#include <iostream>
#include <string>
#include "molecule_list.h"


//add in LAMMPS functions from smac

namespace TetrahedralParticlesInConfinement{
    
    void load(const char*, coord_list_t&, arg_t=0x0);
    void loadxyz(const char* filename, coord_list_t& x, arg_t=0x0);
    void loadxyz(const char* filename, MoleculeList& system, Box& box);
    
    struct xyz_info{
        xyz_info() : outstream(0x0),instream(0x0){
            xres = coord_t(3, 0.0);
        }
        
        Box box;
        std::vector<std::string> type;
        std::ostream* outstream;
        std::istream* instream;
        std::map<std::string, int> reservoir;
        coord_t xres;
    };
    
    struct xyzfile{
        unsigned int n;
        coord_list_t x;
        std::string commentstr;
        std::vector<std::string> type;
    };
    
    void loadxyz(const char*, xyzfile&);
    void loadxyz(std::istream&, xyzfile&);
    
    void savexyz(const char* filename, coord_list_t& x, xyz_info&);
    void savevarxyz(const char* filename, coord_list_t& x, xyz_info&);
    
    std::ostream& operator << (std::ostream&, Box&);
    std::ostream& operator << (std::ostream&, function1d_t&);
    std::ostream& operator << (std::ostream&, coord_list_t&);
    std::ostream& operator << (std::ostream&, coord_t&);
    
    //std::ostream& operator << (std::ostream&, shpdesc_t&);
    
    std::istream& operator >> (std::istream&, coord_list_t&);
    std::istream& operator >> (std::istream&, coord_t&);
    //std::istream& operator >> (std::istream&, shpdesc_t&);
    
    void deltype(coord_list_t& x, std::vector<std::string>& types, std::string);
}

#endif /* defined(__TetrahedralParticlesInConfinement__io__) */
