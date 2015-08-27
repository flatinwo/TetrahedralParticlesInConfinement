//
//  simulation.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "simulation.h"
#include "preprocess.h"
#include <cassert>
#include <cmath>

//note deleting destructor so that i can have implicitly defined copy xtor

namespace TetrahedralParticlesInConfinement{
    
    Simulation::Simulation(MoleculeList& molecule_list, Box& box, RandomNumberGenerator& rng):_molecule_list(molecule_list),_box(box),_rng(rng){
        _pair_info = pair_info();
        _steps = 0;
        _max_displacement = 0.;
        buildNeighborList();
    }
    
    
    Simulation::~Simulation(){
        
    }
    
    void Simulation::run(int nstep){
        assert(0);
    }
    
    
#pragma mark SETS

    void Simulation::setMoleculeCoord(TetramerPatchyColloid& molecule, int index){
        _molecule_list.molecule_list[index] = molecule;
    }
    
    void Simulation::setMoleculeList(MoleculeList& molecule_list){
        assert(0);
        //_molecule_list = molecule_list; //define assignment operator
    }
    
    void Simulation::setNMovesPerCycle(double nmovespercycle){
        _nmovespercycle = nmovespercycle;
    }
    
    void Simulation::setEquilibrate(bool flag){
        _equilibrate = flag;
    }
    
    void Simulation::setPairInfo(pair_info& info){
        _pair_info = info;
    }
    
    void Simulation::setNeighborInfo(neighbor_list_info& info){
        _neighbor_list.second = info;
    }
    
    void Simulation::resetSteps(){
        _steps = 0;
    }
#pragma mark GETS
    
    double Simulation::getAcceptanceProbability(){
        assert(0);
        return 0.0;
        
    }
    
    double Simulation::getNMovesPerCycle(){
        return _nmovespercycle;
    }
    
    MoleculeList& Simulation::getMoleculeList(){
        return _molecule_list;
    }
    
    Box& Simulation::getBox(){
        return _box;
    }
    
    //to be tested
    void Simulation::buildNeighborList(){
        _coords_since_last_neighbor_build = _molecule_list.getFullColloidListCoord();
        ::TetrahedralParticlesInConfinement::build_neighbor_list(_coords_since_last_neighbor_build, _box, _neighbor_list);
        assert(_molecule_list.full_colloid_list.size() == _neighbor_list.first.size());
    }
    
    double Simulation::computeEnergy(int index){
        assert(0);
        return 0.;
        
    }
    
    double Simulation::computeEnergy(){
        assert(0);
        return 0.;
    }
    
#pragma mark OTHERS
    
    void Simulation::addToMoleculeList(const coord_t& x_c){
        _molecule_list.addToMoleculeList(x_c);
        buildNeighborList();
    }
    
    
    void Simulation::addMolecule(){
        
        int dim = (int) _box.box_lo.size();
        coord_t xc(dim);
        
        for (int i=0; i<3; i++) {
            xc[i] = _box.box_lo[i] + _rng.randDouble()*_box.box_period[i];
        }
        
        addToMoleculeList(xc);
    }
    
    void Simulation::removeMolecule(int i){
        TetramerPatchyColloid temp;
        int n_elements = (int)_molecule_list.molecule_list.size();
        
        temp = _molecule_list.molecule_list[i];
        _molecule_list.molecule_list[i] = _molecule_list.molecule_list[n_elements-1];
        _molecule_list.molecule_list[n_elements-1] = temp;
        
        _molecule_list.popBackMolecule();
        buildNeighborList();
        
        
    }
    
    void Simulation::removeMolecule(){
        _molecule_list.popBackMolecule();
        buildNeighborList();
    }
    
#pragma mark OPERATORS
    std::ostream& operator << (std::ostream& _os, const Simulation& _simulation){
        //wrap to box first
        coord_list_t x = _simulation._molecule_list.getFullColloidListCoord();
        
        assert(x.size() > 0);
        
        unmapwithfloor(x, _simulation._box.box_lo, _simulation._box);
        //unmap(x, _simulation._box.box_lo, _simulation._box);
        
        _os << x.size() << "\n";
        _os << "Atoms. MC Step: " << _simulation._steps << "\tBox Period\t" << _simulation._box.box_period << std::endl;
        
        //this still needs work
        for (unsigned int i=0; i<x.size(); i++){
            if (i%5==0) {
                _os << "\tN\t" << x[i] << "\t" << _simulation._molecule_list.full_colloid_list[i]->orientation << std::endl;
            }
            else
                _os << "\tO\t" << x[i] << "\t" << _simulation._molecule_list.full_colloid_list[i]->orientation << std::endl;
        }
        return _os;
        
    }
    
    
}