//
//  lattice_fcc.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__lattice_fcc__
#define __TetrahedralParticlesInConfinement__lattice_fcc__

#include <stdio.h>
#include "struct_def.h"

namespace TetrahedralParticlesInConfinement{
    
    class LatticeFCC : public Lattice{
    public:
        LatticeFCC();
        ~LatticeFCC();
        void generateLattice();
    };

}

#endif /* defined(__TetrahedralParticlesInConfinement__lattice_fcc__) */
