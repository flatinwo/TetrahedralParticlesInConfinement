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
    TPiC::TPiC(){
        _initialize();
        _setup();
        _nullAllPtrs();
        _nvt.reset(new SimulationNVTEnsemble(_system,_box,_rng));
        _nvt->setCosAngleMax(_cosAngleMax);
        _nvt->setBeta(1.0/_temperature);
        
    }
    
    TPiC::TPiC(std::vector<std::string>& input){
        _command_list = input;
    }
    
    TPiC::TPiC(const char* filename){
        std::ifstream readfile(filename);
        std::string str;
        while (!readfile.eof()){
            std::getline(readfile, str,'#');
            _command_list.push_back(str);
        }
        
    }
    void TPiC::setSimulationMode(CalculationMode mode){
        _mode = mode;
    }
    
    void TPiC::_initialize(){
        _bondLength = 0.5;
        _temperature = 0.15;
        _density = 0.1;
        _nMolecules = 25;
        _cosAngleMax = 0.9;
        _mode = AUTOSTART;
        _ncyclesEquilibration = 5000;
        _ncyclesProduction = 10000;
        _ncyclesAnalysis=1000;
        _configOutputFrequency = 100;
        
    }
    
    void TPiC::_setup(){
        _system.setMoleculeListBondLength(_bondLength);
        
        //Set up desired lattice... maybe make a pointer to lattice
        //so that the user can specify the lattice type
        LatticeFCC lattice;
        lattice.setNumberOfLatticePoints(_nMolecules);
        lattice.setDensity(_density);
        lattice.generateLattice();
        _system.buildMoleculeListAndBox(lattice, _box);
        
    }
    
    void TPiC::_nullAllPtrs(){
        _plates = nullptr;
        _density_bias = nullptr;
        _q6_bias = nullptr;
        _nvt = nullptr;
        _npt = nullptr;
        _umbrella = nullptr;
        _analysis = nullptr;
    }
    
    void TPiC::_nvt_equilibrate(){
        std::cout << " Pre Equilibration Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
        std::cout << "Equilibrating....\n";
        _nvt->setEquilibrate(true);
        _nvt->run(_ncyclesEquilibration);
        std::cout << "Post Equilibration Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
    }
    
    void TPiC::_nvt_analysis(){
        std::cout << " Pre Analysis Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
        std::cout << "Analyzing (determining optimal step size)....\n";
        AnalyzeSimulationStepSize analysis(*_nvt);
        analysis.run(_ncyclesAnalysis);
        std::cout << "Post Analysis Energy is\t" << _nvt->computeEnergy()/((double) _nMolecules) << std::endl;
    }
    
    void TPiC::run(){
        _ofile.open("initial_config.xyz");
        _ofile << *_nvt;
        _ofile.close();
        
        _nvt_equilibrate();
        _nvt_analysis();
        _nvt_equilibrate();
        
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