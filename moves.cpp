//
//  moves.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "moves.h"
#include "spatial.h"
#include <cassert>
#include "operator.h"

//do not copy over translate
//make return type void
//and use pass by reference to update values

namespace TetrahedralParticlesInConfinement {
    
    //io operator for move_info
    std::ostream& operator << (std::ostream& _os, const move_info& info){
        _os << "\ndelta_move\t" << info.delta_move << std::endl;
        _os << "move_probability\t" << info.move_probability << std::endl;
        return _os;
    }
    
    
    //update trajectory of individual particle
    void translate(coord_t& x, Box& box, move_info& translate_info){
        
        for (unsigned int i=0;i<x.size();i++){
            x[i] += (translate_info.delta_move)*gasdev();
            pbc(x[i], box.box_period[i], box.periodic[i]);
        }
    }
    
    
    //update trajectories of total system
    void translate(coord_list_t& x, Box& box, move_info& translate_info){
        
        for (unsigned int i=0;i<x.size();i++){
            translate(x[i],box,translate_info);
        }
    }
    
    void rotate(coord_t& orientation, move_info& rotate_info){
        
        assert(orientation.size()==3);
        generateQuaternionPair(rotate_info, rotate_info.result);
        
        for (int i=0; i<3; i++) {
            double Rdx = 0;
            for (int j=0; j<3; j++) {
                Rdx += rotate_info.result.second[i][j]*orientation[j];
            }
            orientation[i] = Rdx;
        }
    }
    
    //test
    void rotate(Colloid& colloid1, move_info& rotate_info){
        assert(colloid1.orientation.size() == 3);
        assert(colloid1.quaternion.size() == 4);
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        colloid1.quaternion = hamilton_product(rotate_info.result.first, colloid1.quaternion); //order matters
        matrix_vector_product(rotate_info.result.second, colloid1.orientation);
        

        
    }
    
    void rotate(Colloid& colloid1, Colloid& colloid_ref, move_info& rotate_info){
        
        assert(colloid1.orientation.size() == 3);
        assert(colloid1.quaternion.size() == 4);
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        colloid1.quaternion = hamilton_product(rotate_info.result.first, colloid1.quaternion); //order matters
        matrix_vector_product(rotate_info.result.second, colloid1.orientation);
        
        for (int i=0; i<3; i++)
            colloid1._center_of_mass[i] = colloid_ref._center_of_mass[i] + colloid1.orientation[i];
        
    }
    
    void rotate(Colloid& colloid1, Colloid& colloid_ref, double bond_length, move_info& rotate_info){
        assert(colloid1.orientation.size() == 3);
        assert(colloid1.quaternion.size() == 4);
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        colloid1.quaternion = hamilton_product(rotate_info.result.first, colloid1.quaternion); //note order of the product matters
        matrix_vector_product(rotate_info.result.second, colloid1.orientation);
        
        for (int i=0; i<3; i++)
            colloid1._center_of_mass[i] = colloid_ref._center_of_mass[i] + bond_length*colloid1.orientation[i];
        
    }
    
    //rotates "patches" on core
    void rotate(Colloid& colloid1, coord_list_t& orientation_list, move_info& rotate_info){
        generateQuaternionPair(rotate_info, rotate_info.result);
        colloid1.quaternion = hamilton_product(rotate_info.result.first, colloid1.quaternion); //note the order of the product matters
        matrix_vector_product(rotate_info.result.second, colloid1.orientation);
        for (unsigned int i=0; i<orientation_list.size(); i++)
            matrix_vector_product(rotate_info.result.second, orientation_list[i]);
        
    }

    
    void rotate(Colloid& colloid1, coord_pair& QR){
        colloid1.quaternion = hamilton_product(QR.first, colloid1.quaternion);
        matrix_vector_product(QR.second, colloid1.orientation);
        matrix_vector_product(QR.second, colloid1._center_of_mass);
        
    }
    
    
    void rotate(TetramerPatchyColloid& molecule, move_info& rotate_info){
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        coord_t molecule_com = molecule.colloid_list[0]._center_of_mass;
        
        int dim = (int) molecule.colloid_list.size();
        
        for (int i=0; i< dim; i++) {
            rotate(molecule.colloid_list[i], rotate_info.result);
            if (i<dim-1) matrix_vector_product(rotate_info.result.second, molecule.orientation_list[i]);
        }
        molecule.setCenterOfMass(molecule_com);
    }
    
#pragma mark with BOX INFO
    
    void rotate(Colloid& colloid1, Box& box, coord_pair& QR){
        colloid1.quaternion = hamilton_product(QR.first, colloid1.quaternion);
        matrix_vector_product(QR.second, colloid1.orientation);
        matrix_vector_product(QR.second, colloid1._center_of_mass); //not needed because this relative to arbitrary zero
        
        //this is not needed either
        for (unsigned int i=0; i<3; i++) {
            pbc(colloid1._center_of_mass[i], box.box_period[i], box.periodic[i]);
        }
        
    }
    
    
    void rotate(TetramerPatchyColloid& molecule, Box& box, move_info& rotate_info){
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        coord_t molecule_com = molecule.colloid_list[0]._center_of_mass;
        
        int dim = (int) molecule.colloid_list.size();
        
        for (int i=0; i< dim; i++) {
            rotate(molecule.colloid_list[i], box, rotate_info.result);
            if (i<dim-1) matrix_vector_product(rotate_info.result.second, molecule.orientation_list[i]);
        }
        molecule.setCenterOfMass(molecule_com, box);
    }
    
    void rotate(Colloid& colloid1, Colloid& colloid_ref, Box& box, move_info& rotate_info){
        
        assert(colloid1.orientation.size() == 3);
        assert(colloid1.quaternion.size() == 4);
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        colloid1.quaternion = hamilton_product(rotate_info.result.first, colloid1.quaternion); //order matters
        matrix_vector_product(rotate_info.result.second, colloid1.orientation);
        
        for (unsigned int i=0; i<3; i++){
            colloid1._center_of_mass[i] = colloid_ref._center_of_mass[i] + colloid1.orientation[i];
            pbc(colloid1._center_of_mass[i], box.box_period[i], box.periodic[i]);
        }
        
    }
    
    void rotate(Colloid& colloid1, Colloid& colloid_ref, double bond_length, Box& box, move_info& rotate_info){
        assert(colloid1.orientation.size() == 3);
        assert(colloid1.quaternion.size() == 4);
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        colloid1.quaternion = hamilton_product(rotate_info.result.first, colloid1.quaternion); //note order of the product matters
        matrix_vector_product(rotate_info.result.second, colloid1.orientation);
        
        for (unsigned int i=0; i<3; i++){
            colloid1._center_of_mass[i] = colloid_ref._center_of_mass[i] + bond_length*colloid1.orientation[i];
            pbc(colloid1._center_of_mass[i], box.box_period[i], box.periodic[i]);
        }
        
    }
    
    
    void rotate_sites_per_molecule(coord_list_t& orientation_list, move_info& rotate_info){
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        double rx, ry, rz;
        for (unsigned int n=0; n<orientation_list.size(); n++) {
            rx = orientation_list[n][0];
            ry = orientation_list[n][1];
            rz = orientation_list[n][2];
            for (int i=0; i<3; i++) {
                orientation_list[n][i] = rotate_info.result.second[i][0]*rx
                                        + rotate_info.result.second[i][1]*ry
                                        + rotate_info.result.second[i][2]*rz;
            }
            
        }
        
    }
    
    
    coord_t rotate_sites(coord_list_t& orientation_list, move_info& rotate_info){
        
        generateQuaternionPair(rotate_info, rotate_info.result);
        double rx, ry, rz;
        
        for (unsigned int n=0; n<orientation_list.size(); n++) {
            rx = orientation_list[n][0];
            ry = orientation_list[n][1];
            rz = orientation_list[n][2];
            for (int i=0; i<3; i++) {
                orientation_list[n][i] = rotate_info.result.second[i][0]*rx
                + rotate_info.result.second[i][1]*ry
                + rotate_info.result.second[i][2]*rz;
            }
        }
        
        return rotate_info.result.first;
        
    }
    
    
    // will this function ever be called?
    void rotate(coord_list_t& orientation_list, move_info& rotate_info){
        
        for (unsigned int i=0;i<orientation_list.size();i++){
            rotate(orientation_list[i],rotate_info);
        }
    }
    
    
    
    coord_pair generateQuaternionPair(move_info& rotate_info){
        
        double s1,s2,s3,p0,p1,p2,p3;
        double v1,v2,v3,theta_quat=0.;
        
        double delt_rot = rotate_info.delta_move;
        
        coord_pair x;
        coord_list_t R;
        
        if (delt_rot == 0.) {
            R.push_back({1.,0.,0.});
            R.push_back({0.,1.,0.});
            R.push_back({0.,0.,1.});
            coord_t temp(4); //perhaps modify to return same q or unit quaternion
            x.first = temp;
            x.second = R;
            return x;
        }
        else{
            do{
                s1 = delt_rot*gasdev();
                s2 = delt_rot*gasdev();
                s3 = delt_rot*gasdev();
                theta_quat = sqrt(s1*s1 + s2*s2 + s3*s3); // can optimize by not taking square root
            }while (theta_quat > M_PI);
            
            v1 = s1/theta_quat;
            v2 = s2/theta_quat;
            v3 = s3/theta_quat;
            
            p0 = cos(0.5*theta_quat);
            double sinthetaquat = sin(0.5*theta_quat);
            p1 = sinthetaquat*v1;
            p2 = sinthetaquat*v2;
            p3 = sinthetaquat*v3;
            
            x.first.push_back(p0);
            x.first.push_back(p1);
            x.first.push_back(p2);
            x.first.push_back(p3);
            
            coord_t temp1;
            temp1.push_back(p0*p0 + p1*p1 - p2*p2 - p3*p3);
            temp1.push_back(2.*(p1*p2 - p0*p3));
            temp1.push_back(2.*(p1*p3 + p0*p2));
            
            R.push_back(temp1);
            
            temp1.clear();
            temp1.push_back(2.*(p1*p2 + p0*p3));
            temp1.push_back(p0*p0 - p1*p1 + p2*p2 - p3*p3);
            temp1.push_back(2.*(p2*p3 - p0*p1));
            
            R.push_back(temp1);
            
            temp1.clear();
            temp1.push_back(2.*(p1*p3 - p0*p2));
            temp1.push_back(2.*(p2*p3 + p0*p1));
            temp1.push_back(p0*p0 - p1*p1 - p2*p2 + p3*p3);
            
            R.push_back(temp1);
            
            x.second = R;
            
            return x;
        }
        
    }
    
    
    void move_info::update_delta_move(){
        
        compute_move_probability();
        
        if (delta_move > delta_move_max) {
            exit(1);
        }
        
        if (move_probability < 0.30) {
            delta_move *= 0.95;
        }
        else if (move_probability > 0.40){
            delta_move *= 1.05;
        }
    }
    
    
    void move_info::compute_move_probability(){
        //assert(total_moves>0);
        if (total_moves == 0) move_probability = 1.00;
        else{
        	move_probability = (double) accepted_moves / (double) total_moves;
            if (move_probability < 0.0){
                std::cout << "Move probability is\t" << move_probability << std::endl;
                std::cout << "Number of accepted moves is\t" << accepted_moves << std::endl;
                std::cout << "Number of total moves is\t" << total_moves << std::endl;
                exit(0);
            }
	}
    }
    
    
    double move_info::get_move_probability(){
        compute_move_probability();
        return move_probability;
    }
    
    void move_info::reset(){
        accepted_moves = 0.0;
        rejected_moves = 0.0;
        total_moves = 0.0;
    }
    
    
#pragma mark RESCALES
    void rescale(coord_t& x, double s){
        int n = (int) x.size();
        assert(n>0);
        for (int i=0; i<n; i++) {
            x[i] *= s;
        }
    }
    
    void rescale(coord_list_t& x, double s){
        int n = (int) x.size();
        if (n == 0) {
            return;
        }
        
        for (int i=0; i<n; i++) {
            rescale(x[i], s);
        }
        
    }
    
    void rescale(Box& box, double s){
        assert(box.box_lo.size()>0);
        rescale(box.box_lo, s);
        rescale(box.box_hi, s);
        rescale(box.box_period, s);
        //Monte Carlo moves on periodicity
    }
    
    
    void rescale(MoleculeList& System, double s){
        int dim = (int) System.molecule_list.size();
        assert(dim>0);
        
        for (int i=0; i<dim; i++) {
            coord_t x = System.molecule_list[i].colloid_list[0]._center_of_mass;
            rescale(x, s);
            System.molecule_list[i].setCenterOfMass(x);
        }
        
    }
    
    void rescale(MoleculeList& System, Box& box, double s){
        rescale(box, s);
        
        int dim = (int) System.molecule_list.size();
        
        assert(dim>0);
        
        for (int i=0; i<dim; i++) {
            coord_t x = System.molecule_list[i].colloid_list[0]._center_of_mass;
            rescale(x, s);
            System.molecule_list[i].setCenterOfMass(x,box);
        }
    }
    
    void rescale2D(Box& box, double s){
        unsigned int dim = (unsigned int) box.box_lo.size();
        
        //Box must be non-periodic in 2D
        assert(!box.periodic[2]); //z-direction
        assert(dim>0);
        
        for (unsigned int i=0; i<dim-1; i++) {
            box.box_lo[i] *= s;
            box.box_hi[i] *= s;
            box.box_period[i] *= s;
        }

    }
    
    void rescale2D(MoleculeList& System, Box& box, double s){
        rescale2D(box, s);
        
        int dim = (int) System.molecule_list.size();
        assert(dim>0);
        for (int i=0; i<dim; i++) {
            coord_t x = System.molecule_list[i].colloid_list[0]._center_of_mass;
            x[0] *= s; x[1] *= s;
            System.molecule_list[i].setCenterOfMass(x,box);
        }
    }
    
    // Generate gaussian random number, mean zero and variance 1
    inline double gasdev(){
        static bool available = false;
        static double gset;
        double fac, rsq, v1, v2;
        
        if (!available) {
            do{
                v1 = 2.0*rand() / double(RAND_MAX) - 1.0;
                v2 = 2.0*rand() / double(RAND_MAX) - 1.0;
                rsq = v1*v1 + v2*v2;
            } while (rsq >= 1.0 || rsq == 0.0);
            fac = sqrt(-2.0*log(rsq)/rsq);
            gset = v1*fac;
            available = true;
            return v2*fac;
        }else{
            available = false;
            return gset;
        }
    }
    
    
    void generateQuaternionPair(move_info& rotate_info, coord_pair& rotate){
        
        double s1,s2,s3,p0,p1,p2,p3;
        double v1,v2,v3,theta_quat=0.;
        
        double delt_rot = rotate_info.delta_move;
        
        assert(rotate_info.result.first.size()==4);
        assert(rotate.second.size()==3);
        
        
        if (delt_rot == 0.) {
            rotate.second[0] = {1.,0.,0.};
            rotate.second[1] = {0.,1.,0.};
            rotate.second[2] = {0.,0.,1.};
            rotate.first = {1.,0.,0.,0.};
            
        }
        else{
            do{
                s1 = delt_rot*gasdev();
                s2 = delt_rot*gasdev();
                s3 = delt_rot*gasdev();
                theta_quat = sqrt(s1*s1 + s2*s2 + s3*s3); // can optimize by not taking square root
            }while (theta_quat > M_PI);
            
            v1 = s1/theta_quat;
            v2 = s2/theta_quat;
            v3 = s3/theta_quat;
            
            p0 = cos(0.5*theta_quat);
            double sinthetaquat = sin(0.5*theta_quat);
            p1 = sinthetaquat*v1;
            p2 = sinthetaquat*v2;
            p3 = sinthetaquat*v3;
            
            rotate.first = {p0,p1,p2,p3};
            
            rotate.second[0]={p0*p0 + p1*p1 - p2*p2 - p3*p3,
                2.*(p1*p2 - p0*p3),2.*(p1*p3 + p0*p2)};
            
            rotate.second[1]={2.*(p1*p2 + p0*p3),
                p0*p0 - p1*p1 + p2*p2 - p3*p3,2.*(p2*p3 - p0*p1)};
            
            rotate.second[3]={2.*(p1*p3 - p0*p2),2.*(p2*p3 + p0*p1),
                p0*p0 - p1*p1 - p2*p2 + p3*p3};
            
        }
        
    }

    
}
