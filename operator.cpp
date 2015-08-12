//
//  operator.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/2/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "operator.h"
#include <cmath>
#include <cassert>

namespace TetrahedralParticlesInConfinement {
    
    void normalize(coord_t& x){
        
        int dim = (int)x.size();
        double sum = 0;
        
        for (int i=0; i<dim; i++) {
            sum += x[i]*x[i];
        }
        
        sum = sqrt(sum);
        for (int i=0; i<dim; i++) {
            x[i] /= sum;
        }
    }
    
    coord_t hamilton_product(coord_t& q1, coord_t& q2){
        coord_t product(4,0);
        
        assert(q1.size() == 4);
        assert(q1.size() == q2.size());
        
        product[0] = q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3];
        product[1] = q1[0]*q2[1] + q1[1]*q2[0] + q1[2]*q2[3] - q1[3]*q2[2];
        product[2] = q1[0]*q2[2] - q1[1]*q2[3] + q1[2]*q2[0] + q1[3]*q2[1];
        product[3] = q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1] + q1[3]*q2[0];
        
        return product;
    }
    
    coord_t multiply_quaternions(coord_t& q1, coord_t& q2){
        return hamilton_product(q1, q2);
    }
    
    double cosine_angle(coord_t& x1, coord_t& x2){
        
        double cosphi = 0.;
        coord_t x1t = x1;
        coord_t x2t = x2;
        
        normalize(x1t);
        normalize(x2t);
        
        assert(x1t.size() == x2t.size());
        
        for (unsigned int i=0; i<x1t.size(); i++) {
            cosphi += x1t[i]*x2t[i];
        }
        
        return cosphi;
    }
    
    void matrix_vector_product(coord_list_t& A, coord_t& x){
        const int dim = (int) x.size();
        assert(A.size() == dim);
        assert(A[0].size() == dim);
        
        coord_t b(dim, 0.);
        
        for (int i=0; i<dim; i++){
            for (int j=0; j<dim; j++) b[i] += A[i][j]*x[j];
        }
        
        x = b;
    }
}
