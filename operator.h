//
//  operator.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/2/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__operator__
#define __TetrahedralParticlesInConfinement__operator__

#include <stdio.h>
#include "struct_def.h"

namespace TetrahedralParticlesInConfinement{
    
    void normalize(coord_t&);
    coord_t hamilton_product(coord_t& q1, coord_t& q2);
    coord_t multiply_quaternions(coord_t& q1, coord_t& q2);
    void matrix_vector_product(coord_list_t& A, coord_t& x);
    
    double cosine_angle(const coord_t&, const coord_t&);
    
}
#endif /* defined(__TetrahedralParticlesInConfinement__operator__) */
