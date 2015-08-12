//
//  random_number_generator_48.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__random_number_generator_48__
#define __TetrahedralParticlesInConfinement__random_number_generator_48__

#include <stdio.h>
#include "random_number_generator.h"

namespace TetrahedralParticlesInConfinement{
    
    class RandomNumberGenerator48 : public RandomNumberGenerator {
    public:
        RandomNumberGenerator48(unsigned short=11, unsigned short=7, unsigned short=19);
        RandomNumberGenerator48(unsigned short seed[3]);
        virtual ~RandomNumberGenerator48();
        
        void setSeed(unsigned short seed[3]);
        double randDouble();
        int randInt();
        
    private:
        unsigned short _seed[3];
    };
}

#endif /* defined(__TetrahedralParticlesInConfinement__random_number_generator_48__) */
