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
        _delta_skin = 0.5*(sqrt(_neighbor_list.second.cut_off_sqd) - sqrt(_pair_info.cut_off_criteria));
        confinement = NULL;
        _cell_list = nullptr;
        buildNeighborList();
    }
    
    
    Simulation::~Simulation(){
        
    }
    
    void Simulation::run(int nstep){
        assert(nstep>0);
        assert(0);
    }
    
    
#pragma mark SETS

    void Simulation::setMoleculeCoord(TetramerPatchyColloid& molecule, int index){
        _molecule_list.molecule_list[index] = molecule;
    }
    
    void Simulation::setMoleculeList(MoleculeList& molecule_list){
        assert(0);
        _molecule_list = molecule_list; //define assignment operator
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
        _delta_skin = 0.5*(sqrt(_neighbor_list.second.cut_off_sqd) - sqrt(_pair_info.cut_off_criteria));
    }
    
    void Simulation::resetSteps(){
        _steps = 0;
    }
    
    void Simulation::setCellList(CellList& info){
        _cell_list = &info;
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
    
    
    pair_info& Simulation::getPairInfo(){
        return _pair_info;
    }
    //to be tested
    void Simulation::buildNeighborList(){
        if (_cell_list != nullptr) {
            buildCellList();
            return;
        }
        _coords_since_last_neighbor_build = _molecule_list.getFullColloidListCoord();
        ::TetrahedralParticlesInConfinement::build_neighbor_list(_molecule_list.getMoleculeListCoord(), _box, _neighbor_list);
        assert(_molecule_list.full_colloid_list.size() == _neighbor_list.first.size());
    }
    
    void Simulation::buildCellList(){
        //assert(_cell_list != nullptr);
        _cell_list->clear();
        _cell_list->setBox(_box);
        _cell_list->setInteractionRange(sqrt(_pair_info.cut_off_criteria));
        _cell_list->setNeighborStyle(CellList::FULL);
        _coords_since_last_neighbor_build = _molecule_list.getFullColloidListCoord();
        for (unsigned int i=0; i<_coords_since_last_neighbor_build.size(); i++) _cell_list->insert(i, _coords_since_last_neighbor_build[i]);
    }
    
    
    double Simulation::computeEnergy(int index){
        assert(index > 0);
        assert(0);
        return 0;
        
    }
    
    double Simulation::computeEnergy(){
        assert(0);
        return 0.;
    }
    
#pragma mark OTHERS
    
    void Simulation::addToMoleculeList(const coord_t& x_c){
        _molecule_list.addToMoleculeList(x_c,_box);
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
    
    
    void Simulation::addMoleculeWithoutOverlaps(){
        
        //allow for max count
        
        addMolecule();
        computeEnergy();
        
        while (_pair_info.overlap) {
            removeMolecule();
            addMolecule();
            computeEnergy();
        }
    }
    
    void Simulation::addPlates(Plates& plate, bool rebuild){
        //first do not allow this to be called if system has been built already
        unsigned int nmolecules = (unsigned int) _molecule_list.molecule_list.size();
        
        if (rebuild) for (unsigned int i=0; i<nmolecules-1; i++) _molecule_list.popBackMolecule();
        
        confinement = &plate; //sets-up plates
        
        //set-up box, i.e. remove periodicity in one of the axis
        assert(confinement->_walls[0].axis == confinement->_walls[1].axis); //make sure we have the walls in the same axis
        assert((unsigned int) confinement->_walls[0].axis < _box.periodic.size());  //make sure axis specification is right
        
        
        int m = confinement->_walls[0].axis;
        _box.periodic[m] = false;
        _box.box_hi[m] = confinement->_walls[1].position;
        _box.box_lo[m] = confinement->_walls[0].position;
        _box.box_period[m] = _box.box_hi[m] - _box.box_lo[m];
        
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
                _os << "\tN\t" << x[i] << "\t" << _simulation._molecule_list.full_colloid_list[i]->quaternion << std::endl;
            }
            else
                _os << "\tO\t" << x[i] << "\t" << _simulation._molecule_list.full_colloid_list[i]->orientation << std::endl;
        }
        return _os;
        
    }
    
    
}
