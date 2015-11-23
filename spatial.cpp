//
//  spatial.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "spatial.h"
#include <cassert>
#include <cmath>

namespace TetrahedralParticlesInConfinement{
    
    double distancesq(coord_t& x1, coord_t& x2, Box& box){
        unsigned long dim = x1.size();

        if (dim != x2.size())
            assert(0);
        
        double sum = 0.0, dx = 0.;
        
        for (unsigned int j=0; j<dim; j++){
            dx = x1[j]-x2[j];
            pbc(dx,box.box_period[j],box.periodic[j]);
            sum += dx*dx;
        }
        return sum;
    }
    
    
    double_coord_t distancesqandvec(coord_t& x1, coord_t& x2, Box& box){
        unsigned long dim = x1.size();
        
        double_coord_t temp;
        
        if (dim != x2.size())
            assert(0);
        
        double sum = 0.0,dx=0.;
        
        for (unsigned int j=0; j<dim; j++){
            dx = x1[j]-x2[j];
            pbc(dx,box.box_period[j],box.periodic[j]);
            sum += dx*dx;
            temp.second.push_back(dx);
        }
        temp.first = sum;
        return temp;
    }
    
    
    double distancesq(coord_t& x1, coord_t& x2, coord_t& box_period, bool_list_t& periodic){
        
        unsigned long dim = x1.size();
        
        if (dim != x2.size())
            assert(0);
        
        double sum = 0.0, dx = 0.;
        
        for (unsigned int j=0; j<dim; j++){
            dx = x1[j]-x2[j];
            pbc(dx,box_period[j],periodic[j]);
            sum += dx*dx;
        }
        return sum;
        
    }
    
    double distancesq(coord_t& x1, coord_t& x2){
        
        unsigned long dim = x1.size();
        
        if (dim != x2.size())
            assert(0);
        
        double sum = 0.0,dx=0.;
        
        for (unsigned int j=0; j<dim; j++){
            dx = x1[j]-x2[j];
            sum += dx*dx;
        }
        return sum;
        
    }
    
    
    coord_t distancevec(coord_t& x1, coord_t& x2, Box& box){
        
        unsigned long dim = x1.size();
        coord_t dx(dim);
        
        if (dim != x2.size())
            assert(0);
        
        for (unsigned int j=0; j<dim; j++){
            dx[j] = x1[j]-x2[j];
            pbc(dx[j],box.box_period[j],box.periodic[j]);
        }
        return dx;
        
    }
    
    double distance(coord_t& x1, coord_t& x2, Box& box){
        return sqrt(distancesq(x1,x2,box));
    }
    
    double distance(coord_t& x1, coord_t& x2){
        return sqrt(distancesq(x1,x2));
    }
    
    void pbc(coord_t& x, const coord_t& period, const bool_list_t& periodic){
        for (unsigned int j=0; j<x.size();j++){
            pbc(x[j],period[j],periodic[j]);
        }
    }
    
    
    inline double anint(double x){
        if (x > 0.5){
            return 1.0;
        }
        else if (x < -0.5){
            return -1.0;
        }
        else{
            return 0.0;
        }
    }
    
    
    void pbc(double& x, double period, bool periodic){
        
        if (periodic){
            //x -= period*anint(x/period);
            x -= period*round(x/period);
        }
        
    }
    
    void pbcwithfloor(double& x, double period, bool periodic){
        if (periodic) {
            x -= period*floor(x/period);
        }
    }
    
    bool are_particles_in_box(coord_t& x, Box& box){
        
        if (x.size() != box.box_hi.size())
            assert(0);
        
        for (unsigned int i=0; i<x.size();i++){
            if (x[i] < box.box_lo[i]){
                return false;
            }
            if (x[i] >= box.box_hi[i]) {
                return false;
            }
        }
        return true;
        
    }
    
    
    bool are_particles_in_box(coord_list_t& x, Box& box){
        for (unsigned int i=0; i<x.size();i++){
            if (!are_particles_in_box(x[i],box)){
                return false;
            }
        }
        return true;
    }
    
    
    
}
