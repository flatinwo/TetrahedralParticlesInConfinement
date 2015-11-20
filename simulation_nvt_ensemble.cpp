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
#include "spatial.h"

//question: it may be better to return on translate or rotate instead
//of passing by reference, because of the acceptance probability
//how often should i build the neighbor list?

//try low density NVT calculation... are there any failures?
//also fix box to always be from zeros

//need to assess rotation move

namespace TetrahedralParticlesInConfinement{
    SimulationNVTEnsemble::SimulationNVTEnsemble(MoleculeList& molecule_list, Box& box, RandomNumberGenerator& rng):Simulation(molecule_list, box,rng),_beta(1.0){
        _n = (int) _molecule_list.getFullColloidListCoord().size();
        
        assert(_n>0);
        _temperature = 1./_beta;
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
        _move_info_map[TRANSLATE].delta_move_max = 0.5*_molecule_list.molecule_list[0].colloid_list[0].diameter;
        
        
        _cos_angle_max = 0.99; //have a set for this
        
        _E = _delE = 0.;
        _core_flag = false;
        
        pressureCalculation = false;
        
        
        _delta_skin = 0.5*(sqrt(_neighbor_list.second.cut_off_sqd) - sqrt(_pair_info.cut_off_criteria));
        assert(_delta_skin > 0.1);
    }
    
    
    SimulationNVTEnsemble::~SimulationNVTEnsemble(){
        
    }
    
    
#pragma mark SETS
    void SimulationNVTEnsemble::setBeta(double beta){
        _beta = beta;
        _temperature = 1./_beta;
        
    }
    
    void SimulationNVTEnsemble::setDensity(double density){
        
        assert(density > 0.);
        
        double old_density = getDensity();
        _density = density;
        double s = pow(old_density/_density,1./3.);
        
        rescale(_molecule_list, _box, s);
        
        computeVolume();
        if (pressureCalculation) {
            _pressure_log.inverseVolume = 1./_volume;
            _pressure_log.refresh();
        }
        buildNeighborList();
    }
    
    void SimulationNVTEnsemble::setUpdateMoveFrequencyPerCycle(int num){
        _update_move_frequency_per_cycle = num*(int)_molecule_list.getFullColloidListCoord().size();
    }
   /*
    void SimulationNVTEnsemble::setUpdateNeighborFrequencyPerCycle(int num){
        double diameter = _molecule_list.molecule_list[0].colloid_list[0].diameter;
        double avg_distance_moved =
        (double) diameter * _move_info_map[TRANSLATE].delta_move;
        
        assert(avg_distance_moved < diameter);
    }
    */
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
    
    void SimulationNVTEnsemble::setPressureCalculation(bool flag){
        pressureCalculation = flag;
        if (pressureCalculation) {
            _pressure_log.reset();
            _pressure_log.inverseVolume = 1./getVolume();
            _pressure_log.refresh();
        }
    }
    
    void SimulationNVTEnsemble::setPressureCalculationMode(PressureMode mode){
        _pressure_log._pressure_config = mode;
    }
    
    void SimulationNVTEnsemble::setPressureCalculationScaleFactor(double dv){
        assert(dv > 0. && dv < 1.);
        _pressure_log.scale_factor = dv;
        _pressure_log.inverseVolume = 1./getVolume();
        _pressure_log.refresh();
    }
    
    
#pragma mark RESETS
    void SimulationNVTEnsemble::resetPressureTally(){
        _pressure_log.reset();
    }
    
#pragma mark GETS
    double SimulationNVTEnsemble::getTemperature(){
        return _temperature;
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
    
    
    std::map<int,move_info>& SimulationNVTEnsemble::getMoveInfoMap(){
        return _move_info_map;
    }
    
    double SimulationNVTEnsemble::getCosAngleMax(){
        return _cos_angle_max;
    }
    
    double SimulationNVTEnsemble::getPressure(){
        assert(_pressure_log.count > 0);
        return _pressure_log._inverse_scale_factor*_temperature*(log(_pressure_log.pressure_sum /(double) _pressure_log.count));
    }
    
#pragma mark OTHERS
    void SimulationNVTEnsemble::run(int nstep){
        
        // UpdateNeighborFrequencyPerCycle();
        assert(_molecule_list.getFullColloidListCoord().size() == _molecule_list.full_colloid_list.size());
        _n = (int) _molecule_list.full_colloid_list.size();
        
        int nmolecules = (int) _molecule_list.molecule_list.size();
        _nmovespercycle = _nsubmoves*nmolecules;
        
        _E = computeEnergy();
        for (int i=0; i<nstep; i++) {
            for (int j=0; j< _nmovespercycle; j++) {  //note 1 cycle = _nsubmoves*nmolecules
                int p = _rng.randInt()%_n;
                attemptMove(p);
            }
            
            if (_steps%_update_neighbors_frequency_per_cycle==0) {
                //buildNeighborList();
            }
            _steps++;
            if ((!_equilibrate) && (_steps%100==0)) _ofile_energy << _steps << "\t" <<
                "\t" << _E / (double) nmolecules <<  "\t" << getDensity() << std::endl;
            
            if (pressureCalculation && _steps % _pressure_log.frequency == 0) {
                tallyPressure();
                if (_pressure_log.count % _pressure_log.frequency == 0) {
                    _ofile_pressure << _steps << "\t" << getPressure() << std::endl;
                }
            }
            
        }
        
    }
    
    
    void SimulationNVTEnsemble::run(){
        
        int nmolecules = (int) _molecule_list.molecule_list.size();
        
        for (int j=0; j< _nsubmoves; j++) {  //note 1 cycle = _nsubmoves*nmolecules
            int p = _rng.randInt()%_n;
            attemptMove(p);
        }
        if ((!_equilibrate) && (_steps%(nmolecules*100)==0)) _ofile_energy << _steps/nmolecules << "\t" << _E / (double) nmolecules <<  "\t" << getDensity() << std::endl;
        _steps++;

    }
    
    void SimulationNVTEnsemble::computeVolume(){
        int dim = (int) _box.box_period.size();
        _volume = 1.;
        for (int i=0;i<dim;i++)
            _volume *= _box.box_period[i];
    }
    
    double SimulationNVTEnsemble::computeEnergy(int index){
        _pair_info.overlap = false;
        
        double energy = compute_pair_energy_full(index, _molecule_list,
                                                 _box, _pair_info,
                                                 _neighbor_list.first,
                                                 _neighbor_list.second);
        
        if  (confinement == NULL)
            return energy;
        else
            return (energy + compute_pair_energy_plates(index, _molecule_list, _box, *confinement, _pair_info));
        
        /*return compute_pair_energy(index, _molecule_list,
                                        _box, _pair_info);*/
    }
    
    
    //can use function pointer here
    double SimulationNVTEnsemble::computeEnergy(){
        _pair_info.overlap = false;
        double energy = compute_pair_energy_full(_molecule_list,
                                                _box, _pair_info,
                                                _neighbor_list.first,
                                                _neighbor_list.second);
        if  (confinement == NULL)
            return energy;
        else
            return (energy + compute_pair_energy_plates(_molecule_list, _box, *confinement, _pair_info));
        
        /*return compute_pair_energy(_molecule_list,
                                   _box, _pair_info);*/
        
    }
    
    //use function to compute total energy of alternate systems energy using different properties
    double SimulationNVTEnsemble::computeEnergy(MoleculeList& system, Box& box){
        _pair_info.overlap = false;
        
        double energy = compute_pair_energy(system, box,_pair_info);
        
        
        if  (confinement == NULL)
            return energy;
        else
            return (energy + compute_pair_energy_plates(system, box, *confinement, _pair_info));
        
        
    }
    
    //compute Pressure
    //references: Miguel and Jackson, Mol. Phys., vol. 104,3717-3734, 2006
    //            Miguel and Jackson, J. Chem. Phys., 125, 164109, 2006
    double SimulationNVTEnsemble::computePressure(double dv, PressureMode mode){
        
        double old_e = computeEnergy();
        //double old_v = getVolume();
        
        MoleculeList old_list = _molecule_list;
        old_list.buildFullColloidListPointer();
        Box old_box = _box;
        
        double N = (double)_molecule_list.molecule_list.size();
        
        //assert(dv > 0.);
        
        double scale_volume;
        
        if (mode == COMPRESSION) {
            scale_volume = 1.0 - dv;//can move out since this is usually fixed
            double s = pow(scale_volume,1./3.); //same here
            rescale(old_list, old_box, s);
            double new_e = computeEnergy(old_list,old_box);
            
            if (_pair_info.overlap) {
                return 0.;
            }
            else
                return pow(scale_volume,N)*exp(-1.0*_beta*(new_e-old_e));
                
        }
        else if (mode == EXPANSION){
            scale_volume = 1.0 + dv;
            double s = pow(scale_volume,1./3.);
            rescale(old_list, old_box, s);
            double new_e = computeEnergy(old_list,old_box);
            return pow(scale_volume,N)*exp(-1.0*_beta*(new_e-old_e));
        }
        else
            return 0.;
        
        
    }
    
    void SimulationNVTEnsemble::tallyPressure(){
        _pressure_log.count++;
        _pressure_log.pressure_sum += computePressure(_pressure_log.scale_factor, _pressure_log._pressure_config);
    }
    
    
    double SimulationNVTEnsemble::computeMoleculeEnergy(int index){
        _pair_info.overlap = false;
        assert(index < (int)_molecule_list.molecule_list.size());
        double energy = compute_pair_molecule_energy_full(index, _molecule_list, _box, _pair_info,_neighbor_list.first, _neighbor_list.second);
        
        if (confinement == NULL) {
            return energy;
        }
        else
            return (energy + compute_pair_molecule_energy_plates(index, _molecule_list, _box, *confinement, _pair_info));
        
        //return compute_pair_molecule_energy(index, _molecule_list, _box, _pair_info);
    }
    
#pragma mark MOVES
    
    int SimulationNVTEnsemble::attemptMove(int i){
        
        assert(i < (int) _molecule_list.full_colloid_list.size());
        
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
        
        double old_e = 0.0, new_e = 0.0;
        
        saveConfig(i);
        if (checkNeighborList(i)) buildNeighborList();
       

	//std::cerr << "before\t" << i << "\t" << _flag << "\t" << distancesq(_molecule_list.full_colloid_list[236]->_center_of_mass,
	//			      				      _molecule_list.full_colloid_list[561]->_center_of_mass,
	//			     				       _box) << std::endl;
        
 
        //make moves, translate, rotate a colloid, or rotate molecule
        if (_flag == TRANSLATE) {
            int molecule_id = _molecule_list.full_colloid_list[i]->molecule_id;
            old_e = computeMoleculeEnergy(molecule_id);
            if (old_e > 749 && !_equilibrate){
                std::cout << computeEnergy() << std::endl;
                std::cout << *this;
                buildNeighborList();
                old_e = computeMoleculeEnergy(molecule_id);
                std::cout << i << "\t" << old_e << "\t" << computeEnergy() << std::endl;
            }
            Translation(i);
            if (checkNeighborList(i)) buildNeighborList();
            new_e = computeMoleculeEnergy(molecule_id);
        }
        else if (_flag == ROTATE ){
            
            if (!_core_flag) {
                old_e = computeEnergy(i);
                if (old_e > 749  && !_equilibrate){
                    std::cout << computeEnergy() << std::endl;
                    std::cout << *this;
                    buildNeighborList();
                    old_e = computeEnergy(i);
                    std::cout << i << "\t" << old_e << "\t" << computeEnergy() << std::endl;
                }
                Rotation(i);
                if (!isRotationGood(i)) {
                    updateMoveInfo(REJECT);
                    revertConfig(i);
                    return 0;
                }
                if (checkNeighborList(i)) buildNeighborList();
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
        else if (_flag == ROTATEMOLECULE){
            int molecule_id = _molecule_list.full_colloid_list[i]->molecule_id;
            old_e = computeMoleculeEnergy(molecule_id);
            if (old_e > 749  && !_equilibrate){
                std::cout << computeEnergy() << std::endl;
                std::cout << *this;
                buildNeighborList();
                old_e = computeMoleculeEnergy(molecule_id);
                std::cout << i << "\t" << old_e << "\t" << computeEnergy() << std::endl;
            }
            Rotation(i);
            if (checkNeighborList(i)) buildNeighborList();
            new_e = computeMoleculeEnergy(molecule_id);
        }
        else std::cerr << "ERROR in SimulationNVTEnsemble: Unknown flag for setting _flag\n";
        
        //std::cerr << old_e << "\t" << new_e << "\n" ;//<< exp(-_beta*(new_e - old_e)) << std::endl;
        
        
        if (old_e > 749 && !_equilibrate) {
            std::cout << computeEnergy() << "\t" << computeEnergy(i) << "\t" << i << "\t" << _flag <<  std::endl;
            buildNeighborList();
            std::cout << computeEnergy() << "\t" << computeEnergy(i) << "\t" << i << "\t" << _flag <<  std::endl;
            std::cout << *this;
            std::cerr << "Initial Energy during equilibration run is too high" << std::endl;
            std::cerr << "You may need to increase the size of your Neighbor List criterion" << std::endl;
            exit(0);
        }
        
        old_flag.first = _flag;
        old_flag.second = i;
        
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
        
        _new_config = _molecule_list.full_colloid_list[index]->_center_of_mass;
        translate(_new_config, _box, _move_info_map[TRANSLATE]);
        _molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id].setCenterOfMass(_new_config,_box);
        
    }
    
    void SimulationNVTEnsemble::Rotation(int index){
        if (_flag == ROTATEMOLECULE) {
            //rotate(_molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id], _move_info_map[ROTATEMOLECULE]);
            rotate(_molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id], _box, _move_info_map[ROTATEMOLECULE]);
        }
        else{
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            if (_core_flag){
                rotate(*_molecule_list.full_colloid_list[index],
                       _molecule_list.molecule_list[molecule_id].orientation_list,
                       _move_info_map[ROTATE]);
            }
            else{
                //rotate(*_molecule_list.full_colloid_list[index],
                //       _molecule_list.molecule_list[molecule_id].colloid_list[0],_molecule_list.getBondLength(),
                //       _move_info_map[ROTATE]);
                
                rotate(*_molecule_list.full_colloid_list[index],
                       _molecule_list.molecule_list[molecule_id].colloid_list[0],_molecule_list.getBondLength(), _box,
                       _move_info_map[ROTATE]);
            }
        }
        
    }
    
    void SimulationNVTEnsemble::TranslationAndRotation(int index){
        Translation(index);
        Rotation(index);
    }
    
    void SimulationNVTEnsemble::saveConfig(int index){
        
        assert(index < (int) _molecule_list.full_colloid_list.size());
        _old_config.clear();
        
        if (_flag == TRANSLATE){
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            _old_molecule = _molecule_list.molecule_list[molecule_id];
            //_old_config.push_back(_molecule_list.full_colloid_list[index]->_center_of_mass);

        }
            //_old_config.push_back(_molecule_list.full_colloid_list[index]->_center_of_mass);
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
            //assert(_old_config.size()==1);
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            _molecule_list.molecule_list[molecule_id] = _old_molecule;
           // int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
           // _molecule_list.molecule_list[molecule_id].setCenterOfMass(_old_config[0],_box);
        }
        else if (_flag == ROTATE){
            int molecule_id = _molecule_list.full_colloid_list[index]->molecule_id;
            int dim = (int) _molecule_list.molecule_list[molecule_id].colloid_list.size();
            int colloid_index = index%dim;
            _molecule_list.molecule_list[molecule_id].colloid_list[colloid_index] = _old_colloid;
            if (_core_flag){
                _molecule_list.molecule_list[_molecule_list.full_colloid_list[index]->molecule_id].orientation_list =
                _old_config;
            }
            
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
    
    void SimulationNVTEnsemble::computeMaxDisplacement(){
        
        assert(_old_config.size()==1);
        assert(_old_config[0].size() == _new_config.size());
        
        _max_displacement = 0.;
        for (unsigned int i=0; i<_old_config[0].size(); i++){
            _max_displacement = std::max(std::abs(_new_config[i] - _old_config[0][i]),_max_displacement); //should this include, try again with pbc?
        }
    }
    
    
    bool SimulationNVTEnsemble::checkNeighborList(MoleculeList& system, Box& box){
        assert(_coords_since_last_neighbor_build.size() == system.full_colloid_list.size());
        
        for (unsigned int molecule_id=0; molecule_id < system.molecule_list.size(); molecule_id++){
            for (unsigned int k=0; k<5; k++) {
                int j = molecule_id*5 + k;
                _max_displacement = std::max(TetrahedralParticlesInConfinement::distance(system.full_colloid_list[j]->_center_of_mass, _coords_since_last_neighbor_build[j],box), _max_displacement);
                
                if (_max_displacement > _delta_skin) {
                    return true;
                }
            }
        }
        return false;

        
    }
    
    bool SimulationNVTEnsemble::checkNeighborList(int i){
        //conservative test to see if neighbor_list should be rebuilt
        //borrowed from f.19 at www.ccl.net allen and tildsey codes
        
        //or code from frenkel & smit
        assert(_coords_since_last_neighbor_build.size() == _molecule_list.full_colloid_list.size());
        
        _max_displacement = TetrahedralParticlesInConfinement::distance(_molecule_list.full_colloid_list[i]->_center_of_mass, _coords_since_last_neighbor_build[i], _box);
        
        if (_max_displacement > _delta_skin)
            return true;
        
        //fix to allow for different numbers of colloids
        if (_flag == ROTATEMOLECULE || _flag == TRANSLATE) {
            int molecule_id = _molecule_list.full_colloid_list[i]->molecule_id;
            for (unsigned int k=0; k<5; k++) {
                int j = molecule_id*5 + k;
                _max_displacement = std::max(TetrahedralParticlesInConfinement::distance(_molecule_list.full_colloid_list[j]->_center_of_mass, _coords_since_last_neighbor_build[j],_box), _max_displacement);
                
                if (_max_displacement > _delta_skin) {
                    return true;
                }
            }
        }
        
        if (_max_displacement > _delta_skin)
            return true;
        else
            return false;
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
    
    void SimulationNVTEnsemble::openPressureFile(){
        if (!pressureCalculation) {
            std::cerr << "ERROR: cannot open pressure file if pressure calculation mode is off\n";
            exit(0);
        }
        std::ostringstream os;
        os << std::setprecision(4);
        os << "pressure_log_T="<<getTemperature()<<"_density="<<getDensity()<<".txt";
        _ofile_pressure.open(os.str().c_str(), std::ios::app);
        
        if (_ofile_pressure.fail()) {
            std::cerr << "ERROR: failed to open pressure log file\n";
            exit(0);
        }
    }
    
    void SimulationNVTEnsemble::closeFile(){
        _ofile.close();
        _ofile_energy.close();
        
    }
    
    void SimulationNVTEnsemble::closePressureFile(){
        _ofile_pressure.close();
    }
    
    void SimulationNVTEnsemble::writeConfig(){
        _ofile << *this;
    }
    
    
}
