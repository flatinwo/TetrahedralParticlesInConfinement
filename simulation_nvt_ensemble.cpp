//
//  simulation_nvt_ensemble.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "simulation_nvt_ensemble.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include "pair.h"
#include "moves.h"
#include "operator.h"

//question: it may be better to return on translate or rotate instead
//of passing by reference, because of the acceptance probability
//how often should i build the neighbor list?

//need to assess rotation move

namespace TetrahedralParticlesInConfinement{
    SimulationNVTEnsemble::SimulationNVTEnsemble(MoleculeList& molecule_list, Box& box, RandomNumberGenerator& rng):Simulation(molecule_list, box,rng),_beta(1.0){
        _n = (int) _molecule_list.getFullColloidListCoord().size();
        
        assert(_n>0);
        computeVolume();
        
        _density = (double)_molecule_list.molecule_list.size()/_volume;
        _flag = TRANSLATE;
        _nsubmoves = 7; //6 rotations (4 arms only, 1 core only, 1 molecule whole), plus 1 translation
        _nmovespercycle *= _nsubmoves;
        
        _update_move_frequency_per_cycle = 10;
        _update_neighbors_frequency_per_cycle = 1;//set based on deltamove
        
        _move_info_map[TRANSLATE] = move_info();
        _move_info_map[ROTATE] = move_info();
        
        _move_info_map[TRANSLATE].delta_move = 0.025;
        _move_info_map[TRANSLATE].delta_move_max = 0.3*_molecule_list.molecule_list[0].colloid_list[0].diameter;
        
        //_update_neighbors_frequency_per_cycle = (int) floor(_colloid_list.colloid_list[0].diameter/_move_info_map[TRANSLATE].delta_move);
        
        /*std::cout << _move_info_map[TRANSLATE].delta_move << "\n";
        std::cout << _move_info_map[ROTATE].delta_move << "\n";
        std::cout << _update_move_frequency_per_cycle << "\n";*/
        
        _cos_angle_max = 0.99; //have a set for this
        
        _E = _delE = 0.;
        _core_flag = false;
    }
    
    
    SimulationNVTEnsemble::~SimulationNVTEnsemble(){
        
    }
    
    
#pragma mark SETS
    void SimulationNVTEnsemble::setBeta(double beta){
        _beta = beta;
        
    }
    
    void SimulationNVTEnsemble::setDensity(double density){
        
        assert(density > 0.);
        
        double old_density = getDensity();
        _density = density;
        double s = pow(old_density/_density,1./3.);
        
        rescale(_molecule_list, s);
        rescale(_box, s);
        
        computeVolume();
        buildNeighborList();
    }
    
    void SimulationNVTEnsemble::setUpdateMoveFrequencyPerCycle(int num){
        _update_move_frequency_per_cycle = num*(int)_molecule_list.getFullColloidListCoord().size();
    }
    
    void SimulationNVTEnsemble::setUpdateNeighborFrequencyPerCycle(int num){
        double diameter = _molecule_list.molecule_list[0].colloid_list[0].diameter;
        double avg_distance_moved =
        (double) diameter * _move_info_map[TRANSLATE].delta_move;
        
        assert(avg_distance_moved < diameter);
    }
    
    void SimulationNVTEnsemble::UpdateNeighborFrequencyPerCycle(){
        _update_neighbors_frequency_per_cycle = (int) floor(_molecule_list.molecule_list[0].colloid_list[0].diameter/_move_info_map[TRANSLATE].delta_move);
    }
    void SimulationNVTEnsemble::setDeltaMoveTranslate(double delta){
        _move_info_map[TRANSLATE].delta_move = delta;
    }
    
    void SimulationNVTEnsemble::setDeltaMoveRotate(double delta){
        _move_info_map[ROTATE].delta_move = delta;
    }
    
    void SimulationNVTEnsemble::setCosAngleMax(double cosangle){
        _cos_angle_max = cosangle;
    }
    
#pragma mark GETS
    double SimulationNVTEnsemble::getTemperature(){
        return 1./(_beta);
    }
    
    double SimulationNVTEnsemble::getDensity(){
        _density = (double) _molecule_list.molecule_list.size() / getVolume();
        return _density;
    }
    
    double SimulationNVTEnsemble::getVolume(){
        computeVolume();
        return _volume;
    }
    
    coord_list_t& SimulationNVTEnsemble::getFullColloidListCoord(){
        return _molecule_list.getFullColloidListCoord();
    }
    
    
    std::map<int,move_info> SimulationNVTEnsemble::getMoveInfoMap(){
        return _move_info_map;
    }
    
    
    
#pragma mark OTHERS
    void SimulationNVTEnsemble::run(int nstep){
        
        // UpdateNeighborFrequencyPerCycle();
        assert(_molecule_list.getFullColloidListCoord().size() == _molecule_list.full_colloid_list.size());
        _n = (int) _molecule_list.full_colloid_list.size();
        
        int nmolecules = (int) _molecule_list.molecule_list.size();
        _nmovespercycle = _nsubmoves*nmolecules;
        
        _E = computeEnergy();
        //std::cout << "Initial Total Energy\t" << _E << std::endl;
        for (int i=0; i<nstep; i++) {
            for (int j=0; j< _nmovespercycle; j++) {  //note 1 cycle = _nsubmoves*nmolecules
                int p = _rng.randInt()%_n;
                attemptMove(p);
            }
            
            if (_steps%_update_neighbors_frequency_per_cycle==0) {
                buildNeighborList();
            }
            _steps++;
            if ((!_equilibrate) && (_steps%100==0)) _ofile_energy << _steps << "\t" << computeEnergy()/(double) nmolecules << "\t" << _E / (double) nmolecules << std::endl;
            //std::cout << "Total Energy\t" << _E << std::endl;
        }
        
    }
    
    void SimulationNVTEnsemble::computeVolume(){
        int dim = (int) _box.box_period.size();
        _volume = 1.;
        for (int i=0;i<dim;i++)
            _volume *= _box.box_period[i];
    }
    
    double SimulationNVTEnsemble::computeEnergy(int index){
        _pair_info.overlap = false;
        return compute_pair_energy_full(index, _molecule_list,
                                   _box, _pair_info,
                                   _neighbor_list.first,
                                   _neighbor_list.second);
    }
    
    double SimulationNVTEnsemble::computeEnergy(){
        _pair_info.overlap = false;
        return compute_pair_energy_full(_molecule_list,
                                   _box, _pair_info,
                                   _neighbor_list.first,
                                   _neighbor_list.second);
        
    }
    
    double SimulationNVTEnsemble::computeMoleculeEnergy(int index){
        _pair_info.overlap = false;
        assert(index < _molecule_list.molecule_list.size());
        return compute_pair_molecule_energy_full(index, _molecule_list, _box, _pair_info,_neighbor_list.first, _neighbor_list.second);
    }
    
#pragma mark MOVES
    
    int SimulationNVTEnsemble::attemptMove(int i){
        
        assert(i < _molecule_list.full_colloid_list.size());
        
        if (!(_molecule_list.full_colloid_list[i]->core)) {
            _flag = ROTATE;
            _core_flag = false;
            return attemptSubMove(i);
        }
        else{
            _flag = _rng.randInt()%3; //option to translate or rotate core or rotate molecule
            _core_flag = true;
            return attemptSubMove(i);
        }
        
    }
    
    int SimulationNVTEnsemble::attemptSubMove(int i){

        //_flag = TRANSLATE;
        
        saveConfig(i);
        
        double old_e = 0.0, new_e = 0.0;
        
        //make moves, translate, rotate a colloid, or rotate molecule
        if (_flag == TRANSLATE) {
            int molecule_id = _molecule_list.full_colloid_list[i]->molecule_id;
            old_e = computeMoleculeEnergy(molecule_id);
            Translation(i);
            new_e = computeMoleculeEnergy(molecule_id);
        }
        else if (_flag == ROTATE ){
            
            if (!_core_flag) {
                old_e = computeEnergy(i);
                Rotation(i);
                if (!isRotationGood(i)) {
                    updateMoveInfo(REJECT);
                    revertConfig(i);
                    return 0;
                }
                new_e = computeEnergy(i);
            }
            else{
                Rotation(i);
                if (!isRotationGood(i)){
                    updateMoveInfo(REJECT);
                    revertConfig(i);
                    return 0;
                }
                else{
                    updateMoveInfo(ACCEPT);
                    return 1;
                }
            }
            
        }
        else if(_flag == ROTATEMOLECULE){
            int molecule_id = _molecule_list.full_colloid_list[i]->molecule_id;
            old_e = computeMoleculeEnergy(molecule_id);
            Rotation(i);
            new_e = computeMoleculeEnergy(molecule_id);
        }
        else std::cerr << "ERROR in SimulationNVTEnsemble: Unknown flag for setting _flag\n";
        
        
        //std::cerr << old_e << "\t" << new_e << "\t" << exp(-_beta*(new_e - old_e)) << std::endl;
        
        if (old_e > 749 && !_equilibrate) {
            std::cerr << "Initial Energy during equilibration run is too high" << std::endl;
            exit(0);
        }
        
        if (_pair_info.overlap){
            updateMoveInfo(REJECT);
            _delE = 0.;
            _E += _delE;
            revertConfig(i);
            return 0;
        }
        
        if (new_e < old_e){
            updateMoveInfo(ACCEPT);
            _delE = new_e - old_e;
            _E += _delE;
            return 1;
        }
        
        else if (_rng.randDouble() < exp(-_beta*(new_e - old_e))) {
            updateMoveInfo(ACCEPT);
            _delE = new_e - old_e;
            _E += _delE;
            return 1;
        }
        
        else {
            updateMoveInfo(REJECT);
            _delE = 0.;
            _E += _delE;
            revertConfig(i);
            return 0;
        }
    }
    
    //can make faster maybe by moving updateMoveInfo to attemptMove if statements
    void SimulationNVTEnsemble::updateMoveInfo(int update_flag){
        
        if (update_flag == ACCEPT) {
            _move_info_map[_flag].accepted_moves++;
            _move_info_map[_flag].total_moves++;
        }
        else if (update_flag == REJECT){
            _move_info_map[_flag].rejected_moves++;
            _move_info_map[_flag].total_moves++;            
        }
        else{
            std::cerr << "ERROR in SimulationNVTEnsemble: Unknown flag for update "
            << "moveinfo _flag\n";
        }
        
    }
    
    //this needs to be checked seriously
    void SimulationNVTEnsemble::Translation(int index){
        assert(_molecule_list.full_colloid_list[index]->core);
        
        coord_t x = _molecule_list.full_colloid_list[index]->_center_of_mass;
        translate(x, _box, _move_info_map[TRANSLATE]);
        _molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id].setCenterOfMass(x);
        
    }
    
    void SimulationNVTEnsemble::Rotation(int index){
        if (_flag == ROTATEMOLECULE) {
            rotate(_molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id], _move_info_map[ROTATE]);
        }
        else{
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            if (_core_flag){
                rotate(*_molecule_list.full_colloid_list[index],
                       _molecule_list.molecule_list[molecule_id].orientation_list,
                       _move_info_map[ROTATE]);
            }
            else{
                rotate(*_molecule_list.full_colloid_list[index],
                       _molecule_list.molecule_list[molecule_id].colloid_list[0],_molecule_list.getBondLength(),
                       _move_info_map[ROTATE]);
            }
        }
        
    }
    
    void SimulationNVTEnsemble::TranslationAndRotation(int index){
        Translation(index);
        Rotation(index);
    }
    
    void SimulationNVTEnsemble::saveConfig(int index){
        _old_config.clear();
        
        if (_flag == TRANSLATE)
            _old_config.push_back(_molecule_list.full_colloid_list[index]->_center_of_mass);
        else if (_flag == ROTATE){
            _old_colloid = *(_molecule_list.full_colloid_list[index]);
            if (_core_flag)
                _old_config = _molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id].orientation_list;
        }
        else if (_flag == ROTATEMOLECULE){
            _old_molecule = _molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id];
        }
        else{
            std::cerr << "ERROR in SimulationNVTEnsemble: Unknown flag for saving config\n";
        }
        
    }
    
    void SimulationNVTEnsemble::revertConfig(int index){
        if (_flag == TRANSLATE){
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            _molecule_list.molecule_list[molecule_id].setCenterOfMass(_old_config[0]);
        }
        else if (_flag == ROTATE){
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            int dim = (int) _molecule_list.molecule_list[molecule_id].colloid_list.size();
            int colloid_index = index%dim;
            _molecule_list.molecule_list[molecule_id].colloid_list[colloid_index] = _old_colloid;
            if (_core_flag)
                _molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id].orientation_list =
                _old_config;
        }
        else if (_flag == ROTATEMOLECULE){
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            _molecule_list.molecule_list[molecule_id] = _old_molecule;
        }
        else{
            std::cerr << "ERROR in SimulationNVTEnsemble: Unknown flag for saving config\n";
        }
        
    }
    
    bool SimulationNVTEnsemble::isRotationGood(int index){
        
        int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
        int dim = (int) _molecule_list.molecule_list[molecule_id].colloid_list.size();
        int colloid_index = index%dim;
        
        if (colloid_index != 0){
            coord_t r1 = _molecule_list.molecule_list[molecule_id].colloid_list[colloid_index].orientation; //arm orientation
            coord_t r2 = _molecule_list.molecule_list[molecule_id].orientation_list[colloid_index-1];
            
            
            double cos_angle = cosine_angle(r1,r2);
            
            return (cos_angle > _cos_angle_max);
        }
        else{
            double cos_angle = 1.0;
            for (int i=1;i<dim; i++){
                coord_t r1 = _molecule_list.molecule_list[molecule_id].colloid_list[i].orientation; //arm orientation
                coord_t r2 = _molecule_list.molecule_list[molecule_id].orientation_list[i-1];
                cos_angle = cosine_angle(r1, r2);
                if (cos_angle < _cos_angle_max) return false;
            }
            return (cos_angle > _cos_angle_max);
        }
    }
    
#pragma mark FILE_HANDLING
    
    void SimulationNVTEnsemble::openFile(){
        
        std::ostringstream os;
        std::ostringstream os_energy;
        
        os << std::setprecision(4);
        os << "NVT_ensemble_T="<<getTemperature()<<"_density="<<getDensity()<<".xyz";
        
        os_energy << std::setprecision(4);
        os_energy << "log_energy_NVT_T="<<getTemperature()<<"_density"<<getDensity()<<".txt";
        
        _ofile.open(os.str().c_str(), std::ios::app);
        _ofile_energy.open(os_energy.str().c_str(), std::ios::app);
        
        if (_ofile.fail() || _ofile_energy.fail()) {
            std::cerr << "ERROR: failed to open NVT files\n";
            exit(0);
        }
    }
    
    void SimulationNVTEnsemble::closeFile(){
        _ofile.close();
        _ofile_energy.close();
        
    }
    
    void SimulationNVTEnsemble::writeConfig(){
        _ofile << *this;
    }
    
    
    
}