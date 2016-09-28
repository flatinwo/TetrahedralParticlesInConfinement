//
//  spatial.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__spatial__
#define __TetrahedralParticlesInConfinement__spatial__

#include <stdio.h>
#include "struct_def.h"

namespace TetrahedralParticlesInConfinement{
    
    //the squared distance between two particles
    double distancesq(coord_t& x1, coord_t& x2);
    
    //the squared distance between two particles
    double distancesq(coord_t& x1, coord_t& x2, coord_t&, bool_list_t&);
    
    //the squared distance between two particles in a box
    double distancesq(coord_t& x1, coord_t& x2, Box&);
    
    //the squared distance and vector distance between two particles in a box
    double_coord_t distancesqandvec(coord_t& x1, coord_t& x2, Box&);
    
    //the squared distance and vector distance between two particles in a box
    void distancesqandvec(coord_t& x1, coord_t& x2, Box&, double_coord_t& temp);
    
    //the distance between two particles
    double distance(coord_t& x1, coord_t&x2);
    
    //the distance between two particles
    double distance(coord_t& x1, coord_t& x2, Box&);
    
    //vector separation between two particles
    coord_t distancevec(coord_t& x1, coord_t& x2, Box&);
    
    //the r, angular theta and phi between two particles
    coord_t spherical_orientation(coord_t& x1, coord_t& x2, Box&);

    
    //are particles in the box
    bool are_particles_in_box(coord_t& x, Box& box);
    
    //are particles in the box
    bool are_particles_in_box(coord_list_t& x, Box& box);

    void pbc(coord_t& x, const coord_t& period, const bool_list_t&);
    
    void pbc(double& x, double period, bool periodic=true);
    void pbcwithfloor(double& x, double period, bool periodic=true);
}

#endif /* defined(__TetrahedralParticlesInConfinement__spatial__) */
