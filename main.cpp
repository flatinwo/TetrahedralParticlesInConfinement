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
//write output in binary to speed up code
//set up calculations so that constructor can be called and set arbitrarily
//make calculations on new config.. saves in efficiency
//why does number of accepted moves become negative.. check memory allocation

using namespace TetrahedralParticlesInConfinement;

int main(int argc, const char * argv[]) {
    TPC Engine;
    Engine.run();
    
    
   /*
    double density = 0.1;
    double temperature = 0.062118;
    double pressure = 0.004602;
    
    if (argc == 3) {
        temperature = atof(argv[1]);
        pressure = atof(argv[2]);
    }
    
    
    
    double bond_length = 0.5;
    MoleculeList System;
    System.setMoleculeListBondLength(bond_length);
    RandomNumberGenerator48 rng, rng2(0,1,3), rng3(0,4,3);
    Box box;
    int nmol = 25;
    
    std::cout << "reduced density\t" << density << std::endl;
    std::cout << "reduced temperature\t" << temperature << std::endl;
    std::cout << "reduced pressure\t" << pressure << std::endl;
    
    
   
    // Build system with lattice
    LatticeFCC lattice;
    lattice.setNumberOfLatticePoints(nmol);
    lattice.setDensity(1.1*density);
    lattice.generateLattice();
    System.buildMoleculeListAndBox(lattice,box);
    
    
    
    // Build NVT system 
    SimulationNVTEnsemble simulation(System, box, rng);
    simulation.setCosAngleMax(0.825);
    simulation.setDensity(density);
    
    //update neighbor list info
    neighbor_list_info n_info;
    n_info.cut_off_sqd = 6.25;
    simulation.setNeighborInfo(n_info);
    
    
    UmbrellaSpring* density_bias = new UmbrellaSpring();
    density_bias->order_parameter = 0.074;
    density_bias->spring_constant = 20000.;
    
    SimulationNPTEnsemble npt(simulation, rng2, pressure);
    
    //UmbrellaSimulation umbrella(npt, rng3, density_bias);
    
    std::cout << box << std::endl;
    
    for (int i=0;i<2;i++){
        
        if (i==0) simulation.setBeta(1./0.15);
        else {
            simulation.setBeta(1./temperature);
            npt.addUmbrellaSpring(*density_bias);
        }
        
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
        for (int i=0; i<500000; i++){
           // if (i%(10*nmol)==0) simulation.writeConfig();
            npt.run(100);
            if (i%(2*nmol)==0) std::cout << simulation.getDensity() << std::endl;
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
    
    //umbrella.run(50000);
    
    delete density_bias;
    */
    return 0;


}