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
    
    //the rsq, angular theta and phi between two particles
    coord_t spherical_orientation(coord_t& x1, coord_t& x2, Box& box){
        coord_t temp(3,0.);
        double_coord_t result = distancesqandvec(x1, x2, box);
        
        double r = sqrt(result.first);
        double theta = acos(result.second[2]/r); //assumes three-d vector
        double phi = atan2(result.second[1], result.second[0]);
        if (phi < 0) phi += 2.0*M_PI;
        
        temp[0] = r; temp[1] = theta; temp[2] = phi;
        
        return  temp;
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
