//
//  analyze_simulation_step_size.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/12/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "analyze_simulation_step_size.h"
#include <iomanip>
#include <cassert>

//write npt analysis

namespace TetrahedralParticlesInConfinement{
    
#define SMALL 0.0001
    
    //Constructors
    AnalyzeSimulationStepSize::AnalyzeSimulationStepSize(SimulationNVTEnsemble& NVT):_NVT(&NVT){
        reset();
        _Gibbs = NULL;
        _NPT = NULL;
    }
    
    AnalyzeSimulationStepSize::AnalyzeSimulationStepSize(SimulationNPTEnsemble& NPT):_NPT(&NPT){
        _NVT = &(_NPT->_NVT);
        _Gibbs = NULL;
        reset();
    }
    
    /*AnalyzeSimulationStepSize::AnalyzeSimulationStepSize(SimulationGibbsEnsemble& Gibbs): _Gibbs(&Gibbs){
        _NVT = &(_Gibbs->_NVT1);
        _NPT = NULL;
        reset();
    }*/
    
    void AnalyzeSimulationStepSize::reset(){
        _iteration_count = 0;
        _acceptance_probability = 0.0;
        _upper_bound = 0.40;
        _lower_bound = 0.30;
        _count_max = 30;
    }
    
#pragma mark SETS
    void AnalyzeSimulationStepSize::setLowerBoundProbability(double bound){
        
        assert(bound > 0.01);
        if (bound < 0.19) {
            std::cerr << "WARNING: It appears the lower bound value has been set to"
            << " a rather low value of \t" << bound << std::endl;
        }
        _lower_bound = bound;
    }
    
    void AnalyzeSimulationStepSize::setUpperBoundProbability(double bound){
        
        assert(bound < 0.99);
        if (bound > 0.80) {
            std::cerr << "WARNING: It appears the upper bound value has been set to"
            << " a rather high value of \t" << bound << std::endl;
        }
        _upper_bound = bound;
        
    }
    
#pragma mark OTHERS
    void AnalyzeSimulationStepSize::run(int nsteps){
        
        _NVT->_equilibrate = false;
        openFile();
        
        if (_NPT == NULL && _Gibbs == NULL) {
            nvt_run(nsteps);
        }
        else if (_NPT != NULL){
            npt_run(nsteps);
        }
        else{
            gibbs_run(nsteps);
        }
        
        closeFile();
        
        std::cerr << "Done with step size analysis\n";
        std::cerr << "Here are your results:\n";
        std::cerr << "TRANSLATE";
        std::cerr << _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE];
        
        std::cerr << "Here are your results:\n";
        std::cerr << "ROTATE";
        _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE].compute_move_probability();
        std::cerr << _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE];
        
        std::cerr << "Here are your results:\n";
        std::cerr << "ROTATEMOLECULE";
        _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE].compute_move_probability();
        std::cerr << _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE];
    }
    
    void AnalyzeSimulationStepSize::nvt_run(int nsteps){
        
        const int count_max = _count_max;
        double old_e, new_e;
        std::vector<double> old_delta_move(3,0.);
        int substeps = 100;
        
        int nsub_steps = ceil((double) nsteps/ (double) substeps);
        
        //here update each delta_max angle at the same time
        
        //while (_acceptance_probability > _upper_bound || _acceptance_probability < _lower_bound) {
        while (areAllMoveInfoGood()) {

            old_e = _NVT->computeEnergy(); //compute energy of initial configuration
            
            for (int i=0; i<nsub_steps; i++) {
                _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].reset();
                _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE].reset();
                _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE].reset();
                
                old_delta_move[0] = _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move;
                old_delta_move[1] = _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE].delta_move;
                old_delta_move[2] = _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE].delta_move;
                
                _NVT->run(substeps);
                
                _acceptance_probability = _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].get_move_probability();
                
                //delta_move updates
                double& delta_move = _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move;
                //double& delta_move1 = _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE].delta_move;
                //double& delta_move2 = _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE].delta_move;
                
                
                //can be made a lot cleaner
                updateDeltaMove(_NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE]);
                updateDeltaMove(_NVT->_move_info_map[SimulationNVTEnsemble::ROTATE]);
                updateDeltaMove(_NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE]);
                //updateDeltaMove(delta_move);
                //updateDeltaMove(delta_move1);
                //updateDeltaMove(delta_move2);
                
                //making sure neighbor lists are updated
                _NVT->_update_neighbors_frequency_per_cycle = floor(((double) _NVT->_molecule_list.full_colloid_list[0]->diameter)/delta_move);
                
                if (_NVT->_update_neighbors_frequency_per_cycle > 10) {
                    _NVT->_update_neighbors_frequency_per_cycle = 10;
                }
                else if (_NVT->_update_neighbors_frequency_per_cycle < 2){
                    _NVT->_update_neighbors_frequency_per_cycle = 2;
                }
                
            }
            
            new_e = _NVT->computeEnergy(); //compute energy of final configuration
            _iteration_count++;
            _ofile << _iteration_count << "\t" << old_delta_move << "\t" << old_e << "\t" << new_e << "\t" << _acceptance_probability << std::endl;
            
            if (_iteration_count > count_max) {
                std::cerr << "Maximum number of allowable iterations reached! \n";
                std::cerr << "Statistics for delta move = \t" << old_delta_move <<"\n============\n";
                std::cerr << "Acceptance probability\t" << _acceptance_probability << std::endl;
                std::cerr << "Initial Energy:\t" << old_e << ",\tFinal Energy:\t" << new_e << std::endl;
                
                std::cerr << "Note done with step size analysis\n";
                std::cerr << "Here are your results:\n";
                std::cerr << "TRANSLATE";
                std::cerr << _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE];
                
                std::cerr << "Here are your results:\n";
                std::cerr << "ROTATE";
                _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE].compute_move_probability();
                std::cerr << _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE];
                
                std::cerr << "Here are your results:\n";
                std::cerr << "ROTATEMOLECULE";
                _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE].compute_move_probability();
                std::cerr << _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE];
                
                exit(0);
            }
        }
        
    }
    
    void AnalyzeSimulationStepSize::npt_run(int nsteps){
        const int count_max = 2*_count_max;
        double old_e, new_e;
        double old_density, new_density;
        
        nvt_run(500);
        reset();
        
        while (_acceptance_probability > _upper_bound || _acceptance_probability < _lower_bound) {
            old_e = _NVT->computeEnergy(); //compute energy of initial configuration
            old_density = _NPT->getDensity();
            
            _NPT->run(nsteps);
            
            new_e = _NVT->computeEnergy(); //compute energy of final configuration, maybe return enthalpy
            new_density = _NPT->getDensity(); // new_density
            
            _acceptance_probability = _NPT->_volume_info.get_move_probability();
            
            _iteration_count++;
            double& delta_move = _NPT->_volume_info.delta_move;
            
            _ofile << _iteration_count << "\t" << delta_move << "\t" << old_e << "\t" << new_e << "\t" <<
            old_density << "\t" << new_density << "\t" << _acceptance_probability << std::endl;
            
            if (_iteration_count > count_max) {
                std::cerr << "Maximum number of allowable iterations reached! \n";
                std::cerr << "Statistics for delta move = \t" << delta_move <<"\n============\n";
                std::cerr << "Acceptance probability\t" << _acceptance_probability << std::endl;
                std::cerr << "Initial Energy:\t" << old_e << ",\tFinal Energy:\t" << new_e << std::endl;
                exit(0);
            }
            updateDeltaMove(delta_move);
            _NPT->_volume_info.reset();
            
        }
        
        reset();
        nvt_run(500);
        reset();
        
    }
    
    void AnalyzeSimulationStepSize::gibbs_run(int nsteps){
        //_Gibbs->run(nsteps);
        
    }
    
    void AnalyzeSimulationStepSize::updateDeltaMove(double& delta_move){
        
        if (_acceptance_probability < _lower_bound) {
            if (_acceptance_probability < 0.15) {
                delta_move /= 3.;
            }
            else
                delta_move *= 0.95;
        }
        else if (_acceptance_probability > _upper_bound){
            if (_acceptance_probability > 0.85) {
                delta_move *= 3.;
            }
            else
                delta_move *= 1.05;
        }
        
    }
    
    void AnalyzeSimulationStepSize::updateDeltaMove(move_info& info){
        
        double _acceptance_prob = info.get_move_probability();
        double& delta_move = info.delta_move;
        
        if (_acceptance_prob < _lower_bound) {
            if (_acceptance_prob < 0.15) {
                delta_move /= 3.;
            }
            else
                delta_move *= 0.95;
        }
        else if (_acceptance_prob > _upper_bound){
            if (_acceptance_prob > 0.85) {
                delta_move *= 3.;
            }
            else
                delta_move *= 1.05;
        }
        
        if (delta_move < SMALL) delta_move = SMALL;
        if (delta_move > info.delta_move_max) delta_move = info.delta_move_max;
        
    }
    
    bool AnalyzeSimulationStepSize::areAllMoveInfoGood(){
        
        if (_iteration_count == 0) {
            return true;
        }
        
        double _acceptance_prob;
        _acceptance_prob = _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].get_move_probability();
        
        if (_acceptance_prob > _upper_bound || _acceptance_prob < _lower_bound){
            return true;
        }
        else{
            _acceptance_prob = _NVT->_move_info_map[SimulationNVTEnsemble::ROTATE].get_move_probability();
            if (_acceptance_prob > _upper_bound || _acceptance_prob < _lower_bound){
                return true;
            }
            else{
                _acceptance_prob = _NVT->_move_info_map[SimulationNVTEnsemble::ROTATEMOLECULE].get_move_probability();
                if (_acceptance_prob > _upper_bound || _acceptance_prob < _lower_bound){
                    return true;
                }
                else
                    return false;
            }
        }
        
    }
    
    
    void AnalyzeSimulationStepSize::updateDeltaMove(){
        
        if (_acceptance_probability < _lower_bound) {
            if (_acceptance_probability < 0.15) {
                _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move *= 0.1;
            }
            else
                _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move *= 0.5;
        }
        else if (_acceptance_probability > _upper_bound){
            if (_acceptance_probability > 0.85) {
                _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move *= 10.;
            }
            else
                _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move *= 1.5;
        }
        
        if (_NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move >=
            0.5*_NVT->_molecule_list.full_colloid_list[0]->diameter) {
            _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].delta_move =
            0.5*(_NVT->_rng.randDouble())*
            (_NVT->_molecule_list.full_colloid_list[0]->diameter);
        }
        
        _NVT->_move_info_map[SimulationNVTEnsemble::TRANSLATE].reset();
    }
    
    void AnalyzeSimulationStepSize::openFile(){
        std::ostringstream os;
        
        os << std::setprecision(4);
        os << "results_for_step_size_at_temperature="<<_NVT->getTemperature()<< "_and_density="
        << _NVT->getDensity()<<".txt";
        _ofile.open(os.str().c_str());
        
        if (_ofile.fail()) {
            std::cerr << "ERROR: failed to open " << os.str().c_str() << std::endl;
            exit(0);
        }
    }
    
    void AnalyzeSimulationStepSize::closeFile(){
        _ofile.close();
    }
}