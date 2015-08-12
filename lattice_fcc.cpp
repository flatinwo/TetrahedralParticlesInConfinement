//
//  lattice_fcc.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "lattice_fcc.h"
#include <cmath>
#include <iostream>

namespace TetrahedralParticlesInConfinement {
    
    LatticeFCC::LatticeFCC():Lattice(){
        lattice_type = "FCC";
        density = 1.0;
    }
    
    LatticeFCC::~LatticeFCC(){
        
    }
    
    void LatticeFCC::generateLattice(){
        // lattice density
        
        int N = number_of_lattice_points;
        box_length = std::pow((double) N/density,1.0/3.0);
        
        //Find an M large enough to fit N atoms on a fcc lattice
        int M = 1;
        while (4*M*M*M < N){//abstract number of atoms in a unit cell
            ++M;
        }
        
        double a = box_length/((double) M); //lattice constant of conventional cell
        
        //shifted unit cell of lattice, update later to vector or coord_t
        double xCell[4] = {0.25,0.75,0.75,0.25};
        double yCell[4] = {0.25,0.75,0.25,0.75};
        double zCell[4] = {0.25,0.25,0.75,0.75};
        
        //place atoms on lattice
        int n=0; // number of atoms placed so far
        for (int x=0;x<M;x++){
            for (int y=0;y<M;y++){
                for (int z=0;z<M;z++){
                    for (int k=0;k<4;k++){
                        if (n<N){
                            coord_t r;
                            r.push_back(((double) x + xCell[k])*a);
                            r.push_back(((double) y + yCell[k])*a);
                            r.push_back(((double) z + zCell[k])*a);
                            ++n;
                            points.push_back(r);
                        }
                    }
                }
            }
        }
    }
    
    
}
