//
//  random_number_generator_48.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "random_number_generator_48.h"
#include <cstdlib>

namespace TetrahedralParticlesInConfinement {
    
    RandomNumberGenerator48::RandomNumberGenerator48(unsigned short s0, unsigned short s1, unsigned short s2){
        _seed[0] = s0;
        _seed[1] = s1;
        _seed[2] = s2;
    }
    
    RandomNumberGenerator48::RandomNumberGenerator48(unsigned short seed[3]){
        setSeed(seed);
    }
    
    RandomNumberGenerator48::~RandomNumberGenerator48(){
        
    }
    
    void RandomNumberGenerator48::setSeed(unsigned short seed[3]){
        _seed[0] = seed[0];
        _seed[1] = seed[1];
        _seed[2] = seed[2];
    }
    
    double RandomNumberGenerator48::randDouble(){
        return erand48(_seed);
    }
    
    int RandomNumberGenerator48::randInt(){
        return int(nrand48(_seed));
    }
}