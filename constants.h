//
//  constants.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef TetrahedralParticlesInConfinement_constants_h
#define TetrahedralParticlesInConfinement_constants_h

#include "struct_def.h"
#include <vector>

namespace TetrahedralParticlesInConfinement {
    
    //tetrahedron
    struct tetrahedron_t{
        coord_list_t orientation_list;
        tetrahedron_t()
        {
            double x0 = sqrt(2./3.);
            double x1 = sqrt(1./3.);
            
            orientation_list.push_back({  x0,  0.,  -x1});
            orientation_list.push_back({ -x0,  0.,  -x1});
            orientation_list.push_back({   0., x0,   x1});
            orientation_list.push_back({   0.,-x0,   x1});
        }
        
    };
    
}

#endif
