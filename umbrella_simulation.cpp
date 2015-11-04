//
//  umbrella_simulation.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/10/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "umbrella_simulation.h"

namespace TetrahedralParticlesInConfinement {
    UmbrellaSimulation::UmbrellaSimulation(SimulationNPTEnsemble& NPT, RandomNumberGenerator& rng, UmbrellaSpring* umbrella):
    _NPT(&NPT), _rng(&rng), _umbrella(umbrella){
        _nstepsMC = 10;
        _umbrella_type = DENSITY;
        _NVT = &(_NPT->_NVT);
        _E = 1000000.;
        _beta = 1./(_NPT->_NVT.getTemperature());
        //_Gibbs = NULL;
        
    }
    
#pragma mark SETS
    void UmbrellaSimulation::setUmbrellaType(int value){
        _umbrella->umbrella_type = value;
    }
    
    void UmbrellaSimulation::setUmbrellaRestraintValue(double value){
        _umbrella->order_parameter = value;
    }
    
    void UmbrellaSimulation::setUmbrellaRestraintSpringConstant(double value){
        _umbrella->spring_constant = value;
    }
    
    void UmbrellaSimulation::setEquilibrate(bool value){
        _equilibrate = value;
    }
    
    void UmbrellaSimulation::resetRunningAverage(){
        _running_avg = 0.;
        _running_count = 0;
    }
    
#pragma mark GETS
    double UmbrellaSimulation::getRunningAverage(){
        return _running_avg/((double) _running_count);
    }
    
#pragma mark OTHERS
    void UmbrellaSimulation::run(int nstep){
        
        //write getter in simulation class to get any desired parameter for umbrella
        //use overloaded function
        
        std::cout << "===============\n";
        for (int i=0; i<nstep; i++) {
            
            
            //run umbrella_less steps
            _NPT->run(_nstepsMC);
            
            double restrain_param = _NPT->getDensity(); //make this pointer to a function;
            
            //compute umbrella energy
            double E_new = _umbrella->getUmbrellaEnergy(restrain_param);
            
            //          std::cout << E_new << "\t" << _E << "\t*****\n";
            
            //MC accept criteria
            if (E_new <= _E) {
                _restrain_value = restrain_param;
                _E = E_new;
                saveConfig();
                _umbrella_info.accepted_moves++;
            }
            else if (exp(-1.*_beta*(E_new - _E)) < _rng->randDouble()){
                _NVT->_box = _old_configs._old_box;
                _NVT->_molecule_list = _old_configs._old_list;
                _NVT->_molecule_list.buildFullColloidListPointer();
                _umbrella_info.rejected_moves++;
            }
            else{
                _restrain_value = restrain_param;
                _E = E_new;
                saveConfig();
                _umbrella_info.accepted_moves++;
            }
            
            if (i%100) std::cout << _NPT->getDensity() << std::endl;
            
            if (!_equilibrate) {
                _ofile << _steps++ << "\t" << _NPT->getDensity() << std::endl;
                _running_avg += _NPT->getDensity();
                _running_count++;
            }
        }
        
    }
    
#pragma mark SAVES
    void UmbrellaSimulation::saveConfig(double value){
        _old_configs._old_density = value;
        
    }
    
    void UmbrellaSimulation::saveConfig(coord_list_t list){
        _old_configs._old_coord_list = list;
    }
    
    void UmbrellaSimulation::saveConfig(){
        _old_configs._old_list = _NVT->getMoleculeList();
        _old_configs._old_list.buildFullColloidListPointer();
        _old_configs._old_box = _NVT->getBox();
    }
    
    
    void UmbrellaSimulation::openFile(){
        
        std::ostringstream os;
        std::string filename;
        
        if (_umbrella_type == DENSITY) {
            filename = "window_DENSITY_DOP=";
        }
        else if (_umbrella_type == Q6){
            filename = "window_Q6_DOP=";
        }
        else{
            std::cerr << "ERROR: Unknown umbrella_type\n";
            exit(0);
        }
        
        os << filename << _umbrella->order_parameter <<".txt";
        _ofile.open(os.str().c_str(), std::ios::app);
        if (_ofile.fail()) {
            std::cerr << "ERROR: failed to open window " << _umbrella->order_parameter << ".\n";
            exit(0);
        }
    }
    
    void UmbrellaSimulation::closeFile(){
        _ofile.close();
    }
    
    
}