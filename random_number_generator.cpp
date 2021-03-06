//
//  random_number_generator.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "random_number_generator.h"
#include <cassert>
#include <cmath>

namespace TetrahedralParticlesInConfinement{
    
    RandomNumberGenerator::RandomNumberGenerator(): _save(false){
        
    }
    
    RandomNumberGenerator::~RandomNumberGenerator(){
        
    }
    
    double RandomNumberGenerator::randUniform(double min, double max){
        double range = max-min;
        return min + range*randDouble();
    }
    
    double RandomNumberGenerator::randNormal(double mean, double stdev){
        double first,v1,v2,rsq,fac;
        
        v1=v2=rsq=0.0;
        
        if (!_save){
            int again = 1;
            while (again){
                v1 = 2.0*randDouble()-1.0;
                v2 = 2.0*randDouble()-1.0;
                rsq = v1*v1 + v2*v2;
                if (rsq < 1.0 && rsq!=0.0) again = 0;
            }
            fac = sqrt(-2.0*log(rsq)/rsq);
            _second = v1*fac;
            first = v2*fac;
            _save = true;
        }
        else {
            first = _second;
            _save = false;
        }
        return (mean+first*stdev);
    }
    
    double RandomNumberGenerator::randDouble(){
        int virtual_function_was_overriden = 0;
        assert(virtual_function_was_overriden);
        return 0.0;
    }
    
    int RandomNumberGenerator::randInt(){
        int virtual_function_was_overriden = 0;
        assert(virtual_function_was_overriden);
        return 0;
    }
}
