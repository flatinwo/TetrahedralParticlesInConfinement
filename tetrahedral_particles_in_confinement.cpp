//
//  tetrahedral_particles_in_confinement.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 12/16/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include <stdio.h>
#include "tetrahedral_particles_in_confinement.h"

namespace TetrahedralParticlesInConfinement {
    TPC::TPC(){
        _initialize();
        setup();
        _nvt.reset(new SimulationNVTEnsemble(_system,_box,_rng));
        _nvt->setCosAngleMax(_cosAngleMax);
        _nvt->setBeta(1.0/_temperature);
    }
    
    TPC::TPC(std::vector<std::string>& input){
        
    }
    
    void TPC::setSimulationMode(SimulationMode mode){
        _mode = mode;
    }
    
    void TPC::_initialize(){
        _bondLength = 0.5;
        _temperature = 0.15;
        _density = 0.1;
        _nMolecules = 25;
        _cosAngleMax = 0.9;
        _mode = AUTOSTART;
        _ncyclesEquilibration = 5000;
        _ncyclesProduction = 10000;
        _configOutputFrequency = 100;
        
    }
    
    void TPC::setup(){
        _system.setMoleculeListBondLength(_bondLength);
        
        //Set up desired lattice... maybe make a pointer to lattice
        //so that the user can specify the lattice type
        LatticeFCC lattice;
        lattice.setNumberOfLatticePoints(_nMolecules);
        lattice.setDensity(_density);
        lattice.generateLattice();
        _system.buildMoleculeListAndBox(lattice, _box);
        
    }
    
    void TPC::run(){
        _ofile.open("initial_config.xyz");
        _ofile << *_nvt;
        _ofile.close();
        
        
        std::cout << "Initial Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
        std::cout << "Equilibrating....\n";
        _nvt->setEquilibrate(true);
        _nvt->run(_ncyclesEquilibration);
        std::cout << "Post Equilibration Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
        
        std::cout << "Analyzing (determining optimal step size)....\n";
        AnalyzeSimulationStepSize analysis(*_nvt);
        analysis.run(1000);
        std::cout << "Post Analysis Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
        
        std::cout << "Equilibrating....\n";
        _nvt->setEquilibrate(true);
        _nvt->run(_ncyclesEquilibration);
        std::cout << "Post Equilibration Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
        
        std::cout << "Production....\n";
        _nvt->setEquilibrate(false);
        _nvt->openFile();
        for (unsigned int i=0; i<_ncyclesProduction; i++) {
            if (i%_configOutputFrequency == 0) {
                _nvt->writeConfig();
                _nvt->run(20);
            }
        }
        _nvt->writeConfig();
        _nvt->closeFile();
        
        std::cout << "Post Production Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
        
    }
}