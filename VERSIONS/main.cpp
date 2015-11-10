//
//  main.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include <iostream>
#include <cassert>
#include "tetrahedral_particles_in_confinement.h"

//to do:
//write move type and analyze frame to make sure expected behavior is occurring

using namespace TetrahedralParticlesInConfinement;

int main(int argc, const char * argv[]) {
   
    double density = 0.12;
    double temperature = 0.072191;
    double pressure = 0.0054;
    
    if (argc == 3) {
        temperature = atof(argv[1]);
        pressure = atof(argv[2]);
    }
    
    double bond_length = 0.5;
    MoleculeList System;
    System.setMoleculeListBondLength(bond_length);
    RandomNumberGenerator48 rng, rng2(0,1,3);
    Box box;
    
    std::cout << "reduced density\t" << density << std::endl;
    std::cout << "reduced temperature\t" << temperature << std::endl;
    std::cout << "reduced pressure\t" << pressure << std::endl;
    
    
   
    /* Build system with lattice*/
    LatticeFCC lattice;
    lattice.setNumberOfLatticePoints(150);
    lattice.setDensity(1.1*density);
    lattice.generateLattice();
    System.buildMoleculeListAndBox(lattice,box);
    
    
    
    /* Build NVT system */
    SimulationNVTEnsemble simulation(System, box, rng);
    simulation.setCosAngleMax(0.825);
    simulation.setDensity(density);
    
    //update neighbor list info
    neighbor_list_info n_info;
    n_info.cut_off_sqd = 12.;
    simulation.setNeighborInfo(n_info);
    
    
    SimulationNPTEnsemble npt(simulation, rng2, pressure);
    
    std::cout << box << std::endl;
    
    for (int i=0;i<2;i++){
        
        if (i==0) simulation.setBeta(1./0.1);
        else simulation.setBeta(1./temperature);
        
        double N = (double) System.molecule_list.size();
        
        std::ofstream file0;
        file0.open("initial_config.xyz");
        file0 << simulation;
        file0.close();
        
        std::cout << "Initial Energy is " << simulation.computeEnergy()/N << std::endl;
        
        std::cout << "Equilibrating..." << std::endl;
        simulation.setEquilibrate(true);
        simulation.run(500);
        
        
        AnalyzeSimulationStepSize analysis(npt);
        
        std::cout << "Analyzing the step size..." << std::endl;
        analysis.run(1000);
        
        
        std::cout << "Equilibrating..." << std::endl;
        simulation.setEquilibrate(true);
        simulation.run(500);
        
        
        for (int i=0; i<3; i++){
            std::cout << i << std::endl;
            simulation.getMoveInfoMap()[i].compute_move_probability();
            std::cout << simulation.getMoveInfoMap()[i];
        }
        
        std::cout << "Initial Energy (post equilibration) is " << simulation.computeEnergy()/N << std::endl;
        
        
        std::cout << "Production..." << std::endl;
        
        if (i==0) continue;
        
        simulation.setEquilibrate(false);
        simulation.openFile();
        for (int i=0; i<50000; i++){
            simulation.writeConfig();
            npt.run(10);
            if (i%100==0) std::cout << simulation.getDensity() << std::endl;
        }
        simulation.writeConfig();
        simulation.closeFile();
        
        
        for (int i=0; i<3; i++){
            std::cout << i << std::endl;
            simulation.getMoveInfoMap()[i].compute_move_probability();
            std::cout << simulation.getMoveInfoMap()[i];
        }
        
        std::cout << "VOLUME INFO" << std::endl;
        npt.getVolumeInfo().compute_move_probability();
        std::cout << npt.getVolumeInfo();
        
        
        std::cout << "Final Energy is " << simulation.computeEnergy()/N << std::endl;
        
    }
    
    std::ofstream file1;
    file1.open("orientation_list.txt");
    for (unsigned int i=0; i<System.molecule_list.size(); i++) {
        for (unsigned int j=0; j<System.molecule_list[i].orientation_list.size(); j++) {
            file1 << i << "\t" << j << "\t" <<
            System.molecule_list[i].orientation_list[j] << std::endl;
        }
    }
    file1.close();
    
    return 0;
   

}