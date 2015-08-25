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
    
    MoleculeList system;
    Box box;
    RandomNumberGenerator48 rng;
    double bond_length = 0.5;
    
    if (argc < 2)
        std::cerr "Invalid number of arguments" << std::endl;
    else
        const char * filename = argv[1];
    
    system.setMoleculeListBondLength(bond_length);
    loadxyz(filename, system, box);
    
    std::cout << "Box details\n" << box;
    SimulationNVTEnsemble NVT(system,box,rng);
    
    std::cout << NVT.computeEnergy()/(double (system.molecule_list.size())) << std::endl;
    
    
    simulation.setBeta(1./0.09);
    double N = (double) system.molecule_list.size();
    
    std::ofstream file0;
    file0.open("initial_config.xyz");
    file0 << simulation;
    file0.close();
    
    std::cout << "Initial Energy is " << simulation.computeEnergy()/N << std::endl;
    std::cout << "Equilibrating..." << std::endl;
    simulation.setEquilibrate(true);
    simulation.run(500);
    
    
    AnalyzeSimulationStepSize analysis(simulation);
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
    
    ///Users/Folarin/Library/Developer/Xcode/DerivedData/
    simulation.setEquilibrate(false);
    simulation.openFile();
    for (int i=0; i<500; i++){
        simulation.writeConfig();
        simulation.run(10);
    }
    simulation.writeConfig();
    simulation.closeFile();
    
    
    for (int i=0; i<3; i++){
        std::cout << i << std::endl;
        simulation.getMoveInfoMap()[i].compute_move_probability();
        std::cout << simulation.getMoveInfoMap()[i];
    }
    std::cout << "Final Energy is " << simulation.computeEnergy()/N << std::endl;
    
    return 0;
}

/*
namespace Temp {
    double density;
    
    if (argc < 2)
        density = 0.13;
        else
            density = atof(argv[1]);
            
            double bond_length = 0.5;
            MoleculeList System;
    System.setMoleculeListBondLength(bond_length);
    
    
    RandomNumberGenerator48 rng;
    Box box;
    
    //update density to density*
    std::cout << "rho_star\t" << density << std::endl;
    density /= pow(1.5 + 0.251*0.5, 3.0);
    std::cout << "rho\t" << density << std::endl;
    
    LatticeFCC lattice;
    lattice.setNumberOfLatticePoints(15);
    lattice.setDensity(2.0*density);
    lattice.generateLattice();
    
    System.buildMoleculeListAndBox(lattice,box);
    
    
    SimulationNVTEnsemble simulation(System, box, rng);
    simulation.setCosAngleMax(0.9);
    simulation.setDensity(density);
    
    std::cout << box << std::endl;
    
    for (int i=0;i<2;i++){
        
        if (i==0) simulation.setBeta(1./0.15);
        else simulation.setBeta(1./0.08);
        
        double N = (double) System.molecule_list.size();
        
        std::ofstream file0;
        file0.open("initial_config.xyz");
        file0 << simulation;
        file0.close();
        
        std::cout << "Initial Energy is " << simulation.computeEnergy()/N << std::endl;
        
        std::cout << "Equilibrating..." << std::endl;
        simulation.setEquilibrate(true);
        simulation.run(500);
        
        
        AnalyzeSimulationStepSize analysis(simulation);
        
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
        
        ///Users/Folarin/Library/Developer/Xcode/DerivedData/
        simulation.setEquilibrate(false);
        simulation.openFile();
        for (int i=0; i<500; i++){
            simulation.writeConfig();
            simulation.run(10);
        }
        simulation.writeConfig();
        simulation.closeFile();
        
        
        for (int i=0; i<3; i++){
            std::cout << i << std::endl;
            simulation.getMoveInfoMap()[i].compute_move_probability();
            std::cout << simulation.getMoveInfoMap()[i];
        }
        std::cout << "Final Energy is " << simulation.computeEnergy()/N << std::endl;
        
    }
}
*/