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
    
    double density;
    
    if (argc < 2)
        density = 0.065;
    else
        density = atof(argv[1]);
    
    double bond_length = 0.5;
    MoleculeList System;
    System.setMoleculeListBondLength(bond_length);
    
    std::cout << density << std::endl;
    
    RandomNumberGenerator48 rng;
    Box box;
    
    LatticeFCC lattice;
    lattice.setNumberOfLatticePoints(150);
    lattice.setDensity(0.1);
    lattice.generateLattice();
    
    System.buildMoleculeListAndBox(lattice,box);
    
    
    SimulationNVTEnsemble simulation(System, box, rng);
    simulation.setCosAngleMax(0.9);
    simulation.setDensity(density);
    simulation.setBeta(1./0.08);
    
    double N = (double) System.molecule_list.size();
    
    std::cout << "Initial Energy is " << simulation.computeEnergy()/N << std::endl;
    
    AnalyzeSimulationStepSize analysis(simulation);
    
    analysis.run(10000);
    
    simulation.setEquilibrate(true);
    simulation.run(10);
    
    
    ///Users/Folarin/Library/Developer/Xcode/DerivedData/
    
    simulation.openFile();
    for (int i=0; i<1; i++){
        simulation.writeConfig();
        simulation.run(100);
    }
    simulation.closeFile();
    
    std::cout << "Final Energy is " << simulation.computeEnergy()/N << std::endl;
    
    return 0;
}
