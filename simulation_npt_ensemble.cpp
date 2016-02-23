//
//  simulation_npt_ensemble.cpp
// :no TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/10/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "simulation_npt_ensemble.h"
#include <cassert>
#include "moves.h"

namespace TetrahedralParticlesInConfinement {
    
    SimulationNPTEnsemble::SimulationNPTEnsemble(SimulationNVTEnsemble& NVT, RandomNumberGenerator& rng, double pressure):
    _NVT(NVT), _rng(rng), _pressure(pressure), _update_volume_move_frequency_per_cycle(10){
        _volume_info = move_info();
        _volume_info.delta_move = 0.005;
        
        vol_move_per_cycle = 2;
        _steps = 0;
        _umbrella = NULL;
        _umbrella_q6 = NULL;
        
    }
    
    //note there is room for understanding delegating constructors
    SimulationNPTEnsemble::SimulationNPTEnsemble(SimulationNVTEnsemble& NVT, RandomNumberGenerator& rng, double pressure, UmbrellaSpring& umbrella):
     _NVT(NVT), _rng(rng), _pressure(pressure), _update_volume_move_frequency_per_cycle(10), _umbrella(&umbrella){
         _volume_info = move_info();
         _volume_info.delta_move = 0.005;
         vol_move_per_cycle = 2;
         _steps = 0;
        
    }
    
    SimulationNPTEnsemble::~SimulationNPTEnsemble(){
        
    }
    
    
#pragma mark SETS
    void SimulationNPTEnsemble::setDensity(double density){
        _NVT.setDensity(density);
    }
    
    void SimulationNPTEnsemble::setPressure(double pressure){
        _pressure = pressure;
    }
    
    void SimulationNPTEnsemble::addUmbrellaSpring(UmbrellaSpring& umbrella){
        if (_umbrella == NULL) _umbrella = &umbrella;
        else {
            _umbrella_q6 = &umbrella;
            oldQ.reset(new BondStructureAnalysis(&old_list,&old_box));
            newQ.reset(new BondStructureAnalysis(&_NVT));
            
            oldQ->setMaxNumberOfNearestNeighbors(4);
            newQ->setMaxNumberOfNearestNeighbors(4);
            oldQ->setRcutOff(3.4);
            newQ->setRcutOff(3.4);
        }
    }
    
#pragma mark GETS
    
    double SimulationNPTEnsemble::getDensity(){
        return _NVT.getDensity();
    }
    
    double SimulationNPTEnsemble::getPressure(){
        return _pressure;
    }
    
    move_info& SimulationNPTEnsemble::getVolumeInfo(){
        return _volume_info;
    }
    
#pragma mark OTHERS
    void SimulationNPTEnsemble::run(int nsteps){
        
        assert(vol_move_per_cycle > 1);
        int nmoleculesplus1 = (int) _NVT.getMoleculeList().molecule_list.size() + 1;
        int nmolecules = nmoleculesplus1 - 1;
        
        assert(nmoleculesplus1 > 2);
        
        _NVT._E = _NVT.computeEnergy();
        
        for (int i=0; i<nsteps; i++) {
            int j = _rng.randInt()%nmoleculesplus1;
            if (j<nmolecules)
                _NVT.run();
            else
                attemptVolumeMoveOptimized();// attemptVolumeMove();
            if (i % nmoleculesplus1 == 0)_steps++;
        }
        
    }
    
    void SimulationNPTEnsemble::sample(){
        std::cout << _steps << "\t" << _NVT.getDensity();
        if (_umbrella_q6 != NULL) std::cout << "\t" << newQ->getQl();
        std::cout << std::endl;
    }
    
    int SimulationNPTEnsemble::attemptVolumeMove(){
        
        double old_e = _NVT.computeEnergy();
        double old_v = _NVT.getVolume();


        old_list = _NVT._molecule_list;
        old_list.buildFullColloidListPointer();
	    old_box = _NVT._box;
 
        if (old_e > 10){
            std::cerr << "Old energy too high in configuration\n";
            std::cout << _NVT;
            std::cout << old_e << std::endl;
            std::cout << _NVT.computeEnergy();
            _NVT.buildNeighborList();
            std::cout << _NVT.computeEnergy();
            exit(2);
        }
        
        double log_new_v = log(old_v) + _volume_info.delta_move*_rng.randNormal();
        double new_v = exp(log_new_v);
        
        double density_ratio = old_v/new_v;
        double old_density = _NVT.getDensity();
        
        double new_density = old_density*density_ratio;
        
        _NVT.setDensity(new_density);
        
        double new_e = _NVT.computeEnergy();
        double Npart = _NVT.getMoleculeList().molecule_list.size();
        
        double arg = (new_e - old_e) + _pressure*(new_v - old_v)
        - ((double) (Npart+ 1))*_NVT.getTemperature()*log(new_v/old_v);
        
        _volume_info.total_moves++;
        
        if (_NVT._pair_info.overlap){
            _NVT._molecule_list = old_list; //fine because reference can only refer to one thing
            _NVT._molecule_list.buildFullColloidListPointer();
	        _NVT._box = old_box;
            _volume_info.rejected_moves++;
            return 0;
        }
        
        if (_rng.randDouble() > exp(-arg/_NVT.getTemperature())) {
            _NVT._molecule_list = old_list;
            _NVT._molecule_list.buildFullColloidListPointer();
	        _NVT._box = old_box;
            _volume_info.rejected_moves++;
            return 0;
        }
        else{
            _NVT._E += (new_e - old_e);
            _volume_info.accepted_moves++;
            return 1;
            
        }
    }
    
    
    
    int SimulationNPTEnsemble::attemptVolumeMoveOptimized(){
        
        double old_e = _NVT.computeEnergy();
        double old_v = _NVT.getVolume();
        
        old_list = _NVT._molecule_list;
        old_list.buildFullColloidListPointer();
        old_box = _NVT._box;
        
        if (old_e > 10){
            std::cerr << "Old energy too high in configuration\n";
            std::cout << _NVT;
            std::cout << old_e << std::endl;
            std::cout << _NVT.computeEnergy();
            _NVT.buildNeighborList();
            std::cout << _NVT.computeEnergy();
            exit(2);
        }
        
        double log_new_v = log(old_v) + _volume_info.delta_move*_rng.randNormal();
        double new_v = exp(log_new_v);
        
        double s = pow(new_v/old_v,1./3.);
        
        //rescale box
        if (!_NVT.confinement) rescale(old_list, old_box, s);
        else rescale2D(old_list, old_box, s);
        
        double new_e = _NVT.computeEnergy(old_list,old_box);
        double Npart = _NVT.getMoleculeList().molecule_list.size();
        
        double arg = (new_e - old_e) + _pressure*(new_v - old_v)
        - ((double) (Npart+ 1))*_NVT.getTemperature()*log(new_v/old_v);
        
        _volume_info.total_moves++;
        
        if (_NVT._pair_info.overlap){
            _volume_info.rejected_moves++;
            return 0;
        }
        
        //Umbrella potential
        if (_umbrella != NULL) {
            //double volume_ratio = new_v/old_v;
            double old_density = _NVT.getDensity();
            double new_density = old_density*old_v/new_v;
            
            arg += (_umbrella->getUmbrellaEnergy(new_density) - _umbrella->getUmbrellaEnergy(old_density)); //it may be more optimal to do this in volume
        }
        
        if (_umbrella_q6 != NULL) {
            oldQ->compute();
            newQ->compute();
            
            std::cout << oldQ->getQl()  << "\t" << newQ->getQl() << std::endl;
            
            arg += (_umbrella_q6->getUmbrellaEnergy(newQ->getQl()) - _umbrella_q6->getUmbrellaEnergy(oldQ->getQl()));
        }
        
        if (_rng.randDouble() > exp(-arg/_NVT.getTemperature())) {
            _volume_info.rejected_moves++;
            return 0;
        }
        else{
            _NVT._E += (new_e - old_e);
            _NVT._molecule_list = old_list;
            _NVT._molecule_list.buildFullColloidListPointer();
            _NVT._box = old_box;
            _NVT._coords_since_last_neighbor_build = _NVT._molecule_list.getFullColloidListCoord();
            _volume_info.accepted_moves++;
            return 1;
            
        }
    }
}
