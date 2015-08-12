//
//  random_number_generator.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__random_number_generator__
#define __TetrahedralParticlesInConfinement__random_number_generator__

#include <stdio.h>
namespace TetrahedralParticlesInConfinement{
    class RandomNumberGenerator{
    public:
        RandomNumberGenerator();
        virtual ~RandomNumberGenerator();
        
        double randUniform(double min=0.0,double max=1.0);
        double randNormal(double mean=0.0, double stdev=1.0);
        
        virtual double randDouble();
        virtual int randInt();
        
    protected:
        bool _save;
        double _second;
    };
}

#endif /* defined(__TetrahedralParticlesInConfinement__random_number_generator__) */
