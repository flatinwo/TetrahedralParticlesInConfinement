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
        _nstepsMC = 20;
        _umbrella_type = DENSITY;
        _NVT = &(_NPT->_NVT);
        _E = 1000000.;
        _beta = 1./(_NPT->_NVT.getTemperature());
        _umbrella1 = NULL;
        _nsampleFrequency = 100;
        //_Gibbs = NULL;
        
    }
    
    
    UmbrellaSimulation::UmbrellaSimulation(SimulationNPTEnsemble& NPT, RandomNumberGenerator& rng, UmbrellaSpring* umbrella, UmbrellaSpring* umbrella1):
    _NPT(&NPT), _rng(&rng), _umbrella(umbrella),_umbrella1(umbrella1){
        _nstepsMC = 20;
        _umbrella_type = DENSITYANDQ6;
        _NVT = &(_NPT->_NVT);
        _E = 1000000.;
        _beta = 1./(_NPT->_NVT.getTemperature());
        _q6analysis.reset(new BondStructureAnalysis(_NVT));
        _q6analysis->setMaxNumberOfNearestNeighbors(4);
        _q6analysis->setRcutOff(3.8);
        _nsampleFrequency = 100;
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
    
    void UmbrellaSimulation::setSampleFrequency(int freq){
        assert(freq > 0);
        _nsampleFrequency = freq;
    }
    
    void UmbrellaSimulation::setMCSweepsperUmbrellaSweep(int freq){
        assert(freq>0.);
        _nstepsMC = freq;
    }
    
#pragma mark RESETS
    void UmbrellaSimulation::resetSteps(){
        _steps = 0;
    }
    
#pragma mark GETS
    double UmbrellaSimulation::getRunningAverage(){
        return 0.;
        //return _running_avg/((double) _running_count);
    }
    
#pragma mark OTHERS
    void UmbrellaSimulation::run(int nstep){
        
        //write getter in simulation class to get any desired parameter for umbrella
        //use overloaded function
        _beta = 1./(_NPT->_NVT.getTemperature());
        
        for (int i=0; i<nstep; i++) {
            
            
            //run umbrella_less steps
            _NPT->run(_nstepsMC);
            
            double restrain_param = _NPT->getDensity(); //make this pointer to a function;
            double newQ6 = _q6analysis->getQl(true);
            
            //compute umbrella energy
            double E_new = _umbrella->getUmbrellaEnergy(restrain_param);
            E_new += _umbrella1->getUmbrellaEnergy(newQ6);
            
            
            //MC accept criteria
            if (E_new <= _E) {
                _restrain_value = restrain_param;
                _old_configs._Q6 = newQ6;
                _old_configs._old_density = restrain_param;
                _E = E_new;
                saveConfig();
                _umbrella_info.accepted_moves++;
                _goodQ6 = newQ6;
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
                _old_configs._Q6 = newQ6;
                _old_configs._old_density = restrain_param;
                saveConfig();
                _umbrella_info.accepted_moves++;
                _goodQ6 = newQ6;
            }
            
            if (!_equilibrate && (_steps%_nsampleFrequency==0)) {
                    _ofile << _steps << "\t" << _old_configs._old_density << "\t" << _old_configs._Q6 << std::endl;
            }
            _steps++;
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
            filename = "window_DENSITY=";
        }
        else if (_umbrella_type == Q6){
            filename = "window_Q6=";
        }
        else if (_umbrella_type == DENSITYANDQ6){
            filename = "window_DENSITY=";
        }
        else{
            std::cerr << "ERROR: Unknown umbrella_type\n";
            exit(0);

        }
        
         os << filename << _umbrella->order_parameter;
        if (_umbrella_type==DENSITYANDQ6) os <<"_Q6="<<_umbrella1->order_parameter;
        os << ".txt";
        
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