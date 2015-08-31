//
//  preprocess.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "preprocess.h"
#include <cmath>
#include "spatial.h"
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "operator.h"

namespace TetrahedralParticlesInConfinement {
    
    /*
     void perturb(coord_list_t& x, double max){
     for (unsigned int i=0; i<x.size(); i++) {
     x[i][0] += randgauss()*max;
     x[i][1] += randgauss()*max;
     x[i][2] += randgauss()*max;
     }
     }
     
     void blur(coord_list_t& x, double max, int npts){
     if (x.size() == 0) {
     return;
     }
     
     coord_list_t x_rot;
     
     int dim = (int) x[0].size();
     for (unsigned int i=0; i<x.size(); i++) {
     for (int j=0; j<npts; j++) {
     coord_t xnew = x[i];
     for (int k=0; k<dim; k++) {
     xnew[k] += randgauss()*max;
     }
     x_rot.push_back(xnew);
     }
     }
     x.insert(x.end(), x_rot.begin(), x_rot.end());
     }
     */
    
    
    /**
     \param x is the coordinates to unmap
     \param box is the simulation box
     
     Unmaps the periodic boundary conditions and places the object
     at the centroid
     */
    
    void unmapcentroid(coord_list_t& x, Box& box){
        coord_t x_c = centroid(x, box);
        unmap(x, x_c, box);
        x_c[0] = -x_c[0]; x_c[1] = -x_c[1]; x_c[2] = -x_c[2];
        translate(x, x_c);
    }
    
    /**
     \param x is the coordinates to unmap
     \param box is the simulation box
     
     Unmaps the periodic boundary conditions and places the object
     at a specified coordinate
     */
    void unmap(coord_list_t& x, coord_t& x_c, Box& box){
        if (x.size() == 0) {
            return;
        }
        int dim = (int) x[0].size();
        for (unsigned int i=0; i<x.size(); i++) {
            for (int k=0; k<dim; k++) {
                double dx = x[i][k] - x_c[k];
                pbc(dx, box.box_period[k], box.periodic[k]);
                x[i][k] = dx;
            }
        }
    }
    
    void unmapwithfloor(coord_list_t&x, coord_t& x_c, Box& box){
        if (x.size() == 0) {
            return;
        }
        
        int dim = (int) x[0].size();
        for (unsigned int i=0; i<x.size(); i++) {
            for (int k=0; k<dim; k++) {
                double dx = x[i][k] - x_c[k];
                pbcwithfloor(dx, box.box_period[k], box.periodic[k]);
                x[i][k] = dx;
            }
        }
    }
    
    
    coord_t centroid(coord_list_t& x){
        int n = (int) x.size();
        
        if (n==0) {
            return coord_t(3,0.0);
        }
        
        int dim = (int) x[0].size();
        
        coord_t xc(dim, 0.0);
        
        for (int i=0; i<n; i++) {
            for (int j=0; j<dim; j++) {
                xc[j] += x[i][j];
            }
        }
        
        for (int j=0; j<dim; j++) {
            xc[j] /= (double) n;
        }
        return xc;
    }
    
    coord_t centroid(coord_list_t& x, Box& box){
        assert(x.size()>0);
        coord_t rc(3,0.0);
        coord_t x_ref = x[0];
        
        for (unsigned int i=0; i<x.size(); i++) {
            coord_t xi(3);
            xi[0] = x[i][0] - x_ref[0];
            xi[1] = x[i][1] - x_ref[1];
            xi[2] = x[i][2] - x_ref[2];
            pbc(xi, box.box_period, box.periodic);
            rc[0] += xi[0];
            rc[1] += xi[1];
            rc[2] += xi[2];
        }
        
        rc[0] /= (double)x.size();
        rc[1] /= (double)x.size();
        rc[2] /= (double)x.size();
        
        rc[0] += x_ref[0];
        rc[1] += x_ref[1];
        rc[2] += x_ref[2];
        
        return rc;
        
    }
    
    void translate(coord_list_t& x, const coord_t& x_c)
    {
        for (unsigned int i=0; i<x.size(); i++) {
            x[i][0] += x_c[0]; x[i][1] += x_c[1]; x[i][2] += x_c[2];
        }
    }
    
    void untranslate(coord_list_t& x, const coord_t& x_c)
    {
        for (unsigned int i=0; i<x.size(); i++) {
            x[i][0] -= x_c[0]; x[i][1] -= x_c[1]; x[i][2] -= x_c[2];
        }
    }
    
    void translate(coord_t& x, const coord_t& x_c){
        for (unsigned int i=0; i<x.size(); i++) {
            x[i] += x_c[i];
        }
    }
    
    void untranslate(coord_t& x, const coord_t& x_c){
        for (unsigned int i=0; i<x.size(); i++) {
            x[i] -= x_c[i];
        }
    }
    
    
    coord_t quaternion(const coord_t& v1, const coord_t& v2){
        assert(v1.size()==3);
        coord_t q(4), v(3);
        v = cross_product(v1,v2);
        
        q[0] = sqrt(dot_product(v1, v1)*dot_product(v2, v2)) + dot_product(v1, v2);
        for (unsigned int i=0; i<v.size(); i++) q[i+1] = v[i];
        return q;
    }
    
    
    void rotateq(coord_list_t& x, std::vector<double> Q)
    {
        assert(Q.size()==4);
        double norm = sqrt(Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]);
        Q[0] /= norm; Q[1] /= norm; Q[2] /= norm; Q[3] /= norm;
        
        double AXX = Q[0] * Q[0] + Q[1] * Q[1] - Q[2] * Q[2] - Q[3] * Q[3];  //modified by FL to match moves.cpp
        double AXY = 2.0 * ( Q[1] * Q[2] - Q[0] * Q[3] ); //modified FL
        double AXZ = 2.0 * ( Q[1] * Q[3] + Q[0] * Q[2] ); //modified FL
        double AYX = 2.0 * ( Q[1] * Q[2] + Q[0] * Q[3] ); //FL
        double AYY = Q[0] * Q[0] - Q[1] * Q[1] + Q[2] * Q[2] - Q[3] * Q[3];
        double AYZ = 2.0 * ( Q[2] * Q[3] - Q[0] * Q[1] ); //FL
        double AZX = 2.0 * ( Q[1] * Q[3] - Q[0] * Q[2] ); //FL
        double AZY = 2.0 * ( Q[2] * Q[3] + Q[0] * Q[1] ); //FL
        double AZZ = Q[0] * Q[0] - Q[1] * Q[1] - Q[2] * Q[2] + Q[3] * Q[3];
        
        int n = (int)x.size();
        for (int i=0; i<n; i++) {
            //create vector between a vertex and the center
            
            norm = sqrt(
                        x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2]);
            
            double unit_vector[3];
            if (norm == 0.0) {
                unit_vector[0] = unit_vector[i] = unit_vector[2] = 0.0;
            }
            else {
                unit_vector[0] = x[i][0] / norm;
                unit_vector[1] = x[i][1] / norm;
                unit_vector[2] = x[i][2] / norm;
            }
            
            double temp[3];
            temp[0] = AXX * unit_vector[0] +
            AXY * unit_vector[1] + AXZ * unit_vector[2];
            temp[1] = AYX * unit_vector[0] +
            AYY * unit_vector[1] + AYZ * unit_vector[2];
            temp[2] = AZX * unit_vector[0] +
            AZY * unit_vector[1] + AZZ * unit_vector[2];
            
            x[i][0] = temp[0] * norm;
            x[i][1] = temp[1] * norm;
            x[i][2] = temp[2] * norm;
            
        }
    }
    
}