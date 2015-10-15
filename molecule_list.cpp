//
//  molecule_list.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "molecule_list.h"
#include "io.h"
#include <cassert>
#include "preprocess.h"

namespace TetrahedralParticlesInConfinement {
    MoleculeList::MoleculeList(): _is_good(false), nparticle_per_molecule(5), _bond_length(1.0){
        //setColloidListOrientation();
    }
    
    
    MoleculeList::~MoleculeList(){
    }
    
    
    /*
    //copy constructor
    MoleculeList::MoleculeList(const MoleculeList& other){
        molecule_list = other.molecule_list;
        nparticle_per_molecule = other.nparticle_per_molecule;
        _bond_length = other._bond_length;
        _is_good = other._is_good;
        _is_bad = other._is_bad;
        _box = other._box;
        _center_of_mass_list = other._center_of_mass_list;
        _center_of_mass_colloid_list = other._center_of_mass_colloid_list;
        _orientation_list = other._orientation_list;
        
        buildFullColloidListPointer(); //check to see if this works
        
    }
    
    //assignment operator
    MoleculeList MoleculeList::operator=(MoleculeList rhs){
        std::swap(full_colloid_list, rhs.full_colloid_list);
        return *this;
    }
    */
    
    
    void MoleculeList::buildMoleculeList(Lattice& lattice_input){
        molecule_list.clear();
        coord_list_t points = lattice_input.getLattice();
        for (unsigned int i=0;i<points.size();i++){
            TetramerPatchyColloid A(_bond_length);
            A.setMoleculeID(i);
            A.setCenterOfMass(points[i]);
            molecule_list.push_back(A);
        }
        buildFullColloidListPointer();
    }
    
    //uses a random orientation, for detailed balance pick random orientation
    void MoleculeList::addToMoleculeList(const coord_t& x_c){
        int n = (int) molecule_list.size();
        TetramerPatchyColloid A(_bond_length);
        A.setCenterOfMass(x_c);
        A.setMoleculeID(n);
        //A.setOrientationRandomly();
        molecule_list.push_back(A);
        buildFullColloidListPointer();
    }
    
    //uses a random orientation, for detailed balance pick random orientation
    void MoleculeList::addToMoleculeList(const coord_t& x_c, const Box& box){
        int n = (int) molecule_list.size();
        TetramerPatchyColloid A(_bond_length);
        A.setCenterOfMass(x_c,box);
        A.setMoleculeID(n);
        //A.setOrientationRandomly();
        molecule_list.push_back(A);
        buildFullColloidListPointer();
    }
    
    
    void MoleculeList::addToMoleculeList(const coord_list_t& x_c, const coord_list_t& orientation, const Box& box){
        int n = (int) molecule_list.size();
        TetramerPatchyColloid A(_bond_length);
        A.setCenterOfMasses(x_c, orientation, box);
        A.setMoleculeID(n);
        molecule_list.push_back(A);
        buildFullColloidListPointer();
    }
    

    
    void MoleculeList::buildMoleculeListAndBox(Lattice& lattice_input, Box& box){
        
        buildMoleculeList(lattice_input);
        
        box.clear();
        double box_length = lattice_input.box_length;
        double offset = 2.00;
        
        //Assumes all particles diameter are the same
        //set box that contains all particles
        //box must start from zero
        double radius = molecule_list[0].colloid_list[0].getRadius();
        
        for (int j=0;j<3;j++){//maybe not three dimensional
            double r = molecule_list[0].colloid_list[0]._center_of_mass[j];
            box.box_lo.push_back( r - offset*radius);
            box.box_hi.push_back(box_length + offset*radius);
            box.box_period.push_back(box_length + 2*offset*radius - r);
        }
        
        
        //Now translate system, so that box_lo will have coordinates 0,0,0
        coord_t x_c = box.box_lo;
        
        untranslate(box.box_lo, x_c);
        untranslate(box.box_hi, x_c);
        
        
        for (unsigned int i=0; i<molecule_list.size(); i++) {
            coord_t x = molecule_list[i].colloid_list[0]._center_of_mass;
            untranslate(x, x_c);
            molecule_list[i].setCenterOfMass(x,box);
        }
        
        _box = box;
        
        
        //now adjust box... so that box origin is in the center
        buildFullColloidListPointer();
        
    }
    
    void MoleculeList::popBackMolecule(){
        molecule_list.pop_back();
        buildFullColloidListPointer(); // do need to worry about pointing to nothing
    }
    
    
#pragma mark GETS
    
    
    double MoleculeList::getBondLength(){
        return _bond_length;
    }
    
    //maybe use unique_ptr
    //This is not what I want, I want to return a reference to the list of center of masses
    //so that they are also update appropriately
    //INEFFFFFFICIENT
    coord_list_t& MoleculeList::getMoleculeListCoord(){
        _center_of_mass_list.clear();
        for (unsigned int i=0; i<molecule_list.size(); i++) {
            _center_of_mass_list.push_back(molecule_list[i].colloid_list[0]._center_of_mass);
        }
        
        return _center_of_mass_list;
        //return cm;
    }
    
    
    coord_list_t& MoleculeList::getFullColloidListCoord(){
        _center_of_mass_colloid_list.clear();
        for (unsigned int i=0; i<molecule_list.size(); i++) {
            for (unsigned int j=0; j<molecule_list[i].colloid_list.size(); j++) {
                _center_of_mass_colloid_list.push_back(molecule_list[i].colloid_list[j]._center_of_mass);
            }

        }
        return _center_of_mass_colloid_list;
        //return cm;
    }
    
    coord_list_t& MoleculeList::getMoleculeListOrientation(){
        return _orientation_list;
    }
    
    coord_list_t& MoleculeList::getMoleculeListOrientation(int index){
        return molecule_list[index].orientation_list;
    }
    
    void MoleculeList::buildFullColloidListPointer(){
        assert(molecule_list.size() > 0);
        full_colloid_list.clear();
        for (unsigned int i=0; i < molecule_list.size(); i++) {
            for (unsigned int j=0; j<molecule_list[i].colloid_list.size(); j++) {
                full_colloid_list.push_back(&molecule_list[i].colloid_list[j]);
            }
            
        }
    }
    
    
#pragma mark SETS
    
    void MoleculeList::setMoleculeListOrientation(coord_list_t& data){
        assert(data.size() == _orientation_list.size());
        _orientation_list.clear();
        for (int i=0; i<molecule_list.size(); i++) {
            _orientation_list.push_back(data[i]);
        }
    }
    
    //needs work
    void MoleculeList::setMoleculeListDiameter(double diameter){
        assert(molecule_list.size() > 0);
        for (unsigned int i=0; i < molecule_list.size(); i++){
            //molecule_list[i].setColloidDiameter(diameter);
        }
    }
    
    void MoleculeList::setMoleculeListOrientation(){
        _orientation_list.clear();
        for (unsigned int i=0; i<molecule_list.size(); i++) {
            _orientation_list.push_back(molecule_list[i].colloid_list[0].orientation); //this needs more work
        }
    }
    
    void MoleculeList::setMoleculeListBondLength(double bond_length){
        _bond_length = bond_length;
    }
    
#pragma mark IO
    std::ostream& operator << (std::ostream& _os, const MoleculeList& system){
        assert(system.molecule_list.size()>0);
        
        _os << (system.nparticle_per_molecule)*(system.molecule_list.size()) << "\n"; //specific to this system
        _os << "Atoms. " << "" << "\n";
        
        for (unsigned int i=0; i<system.molecule_list.size(); i++){
            _os << system.molecule_list[i];
            
        }
        return _os;
    }
    
    
}
