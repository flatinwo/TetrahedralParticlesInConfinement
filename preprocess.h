//
//  preprocess.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__preprocess__
#define __TetrahedralParticlesInConfinement__preprocess__

#include <stdio.h>
#include "struct_def.h"

namespace TetrahedralParticlesInConfinement {
    
    void perturb(coord_list_t& x, double max);
    
    void blur(coord_list_t& x, double max, int npts);
    
    
    void unmap(coord_list_t& x, coord_t& x_c, Box& box);
    void unmapcentroid(coord_list_t& x, Box& box);
    void unmapwithfloor(coord_list_t&x, coord_t& x_c, Box& box); // this is for view, and here we will use the floor command
    
    coord_t centroid(coord_list_t& x, Box& box);
    
    coord_t centroid(coord_list_t& x);
    
    void translate(coord_list_t& x, const coord_t& x_c);
    void translate(coord_t& x, const coord_t& x_c);
    
    void untranslate(coord_list_t& x, const coord_t& x_c);
    void untranslate(coord_t& x, const coord_t& x_c);
    
    /**
     \brief rotate a coordinate list by a specified quaternion
     */
    void rotateq(coord_list_t& x, std::vector<double> Q);
    
    /**
     \brief compute quaternion of rotation between two vectors
     */
    coord_t quaternion(const coord_t& v1, const coord_t& v2);
}

#endif /* defined(__TetrahedralParticlesInConfinement__preprocess__) */
