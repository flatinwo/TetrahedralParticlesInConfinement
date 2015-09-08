//
//  tetramer_patchy_colloid.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "tetramer_patchy_colloid.h"
#include "io.h"
#include "spatial.h"
#include "operator.h"
#include "moves.h"
#include "preprocess.h"
#include <cassert>

namespace TetrahedralParticlesInConfinement {
    
#pragma mark XTORS
    /*TetramerPatchyColloid::TetramerPatchyColloid():_number_of_particles(5){
        buildMolecule();
        //_center_of_mass = &colloid_list[0]._center_of_mass;
        
    }*/
    
    TetramerPatchyColloid::TetramerPatchyColloid(double bond_length):_number_of_particles(5), _bond_length(bond_length){
        buildMolecule();
    }
    
    TetramerPatchyColloid::~TetramerPatchyColloid(){
        
    }
    
#pragma mark SETS
    void TetramerPatchyColloid::setCenterOfMass(const coord_t& x){
        
        assert(x.size() == colloid_list[0]._center_of_mass.size());
        
        coord_t temp = colloid_list[0]._center_of_mass;
        
        colloid_list[0]._center_of_mass = x;
        
        for (unsigned int i=1; i<colloid_list.size(); i++) {
            coord_t dx(x.size(),0);
            
            for (unsigned int j=0; j<x.size(); j++) {
                dx[j] = (colloid_list[i]._center_of_mass[j] - temp[j]);
            }
            
            normalize(dx);
            rescale(dx,_bond_length);
            
            for (unsigned int j=0; j<x.size(); j++) {
                colloid_list[i]._center_of_mass[j] = x[j] + dx[j]; //pbc?
            }
            
        }
        
    }
    
    void TetramerPatchyColloid::setCenterOfMass(const coord_t& x, const Box& box){
        assert(x.size() == colloid_list[0]._center_of_mass.size());
        
        coord_t temp = colloid_list[0]._center_of_mass;
        
        colloid_list[0]._center_of_mass = x;
        for (unsigned int i=1; i<colloid_list.size(); i++) {
            coord_t dx(x.size(),0);
            
            for (unsigned int j=0; j<x.size(); j++) {
                dx[j] = (colloid_list[i]._center_of_mass[j] - temp[j]);
            }
            
            normalize(dx);
            rescale(dx,_bond_length);
            
            for (unsigned int j=0; j<x.size(); j++) {
                pbc(dx[j],box.box_period[j],box.periodic[j]);
                colloid_list[i]._center_of_mass[j] = x[j] + dx[j];
            }

        }
        
    }
    
    
    void TetramerPatchyColloid::setCenterOfMasses(const coord_list_t& x, const coord_list_t& orientation, const Box& box){
        assert((int) x.size() == _number_of_particles);
        
        for (unsigned int i=0; i<x.size(); i++) {
            colloid_list[i]._center_of_mass = x[i];
            if (i > 0) colloid_list[i].orientation = orientation[i];
        }
        
        orientation_list = Template.orientation_list;
        coord_t q = orientation[0];
        
        assert(q.size()==4);
        rotateq(orientation_list, q);
        colloid_list[0].orientation = orientation_list[0];
        colloid_list[0].quaternion = orientation[0];
    }
    
    void TetramerPatchyColloid::setMoleculeID(int index){
        for (unsigned int i=0; i<colloid_list.size(); i++) {
            colloid_list[i].molecule_id = index;
        }
    }
    
    void TetramerPatchyColloid::setBondLength(double length){
        _bond_length = length;
    }
    
#pragma mark BUILDS
    void TetramerPatchyColloid::buildMolecule(){
        
        assert(colloid_list.size()==0);
        
        coord_t unit_quaternion(4,0.);
        unit_quaternion[0] = 1.0;
        int j = 0;
        
        for (int i=0; i<_number_of_particles; i++) {
            Colloid temp;
            temp.quaternion = unit_quaternion;
            if (i!=0) {
                coord_t orientation(Template.orientation_list[j]);
                
                for (unsigned int m=0; m<orientation.size(); m++) {
                    orientation[m] *= _bond_length;
                }
                
                temp.setCenterOfMass(orientation);
                temp.setOrientationVector(Template.orientation_list[j]);
                temp.core = false;
                j++;
            }
            else{
                temp.core = true;
                temp.setOrientationVector(Template.orientation_list[j]);//treat this system as a four patch
            }

            colloid_list.push_back(temp);
        }
        
        orientation_list.resize(_number_of_particles-1);
        assert(orientation_list.size() == Template.orientation_list.size());
        orientation_list = Template.orientation_list; 
        assert(colloid_list.size()==_number_of_particles);
    }
    
    
#pragma mark WRITES
    std::ostream& operator << (std::ostream& _os, const TetramerPatchyColloid& tetramer){
        assert(tetramer.colloid_list.size() == 5);
        //unmapwithfloor(x, _simulation._box.box_lo, _simulation._box);
        
        //_os << tetramer.colloid_list.size() << "\n";
        //_os << "Atoms. MC Step: " << "" << "\n";
        
        for (unsigned int i=0; i<tetramer.colloid_list.size(); i++){
            coord_t x = tetramer.colloid_list[i]._center_of_mass;
            if (i==0) {
                _os << "\tN\t" << x << std::endl;
            }
            else
                _os << "\tO\t" << x << std::endl;
            
        }
        return _os;
    }
}