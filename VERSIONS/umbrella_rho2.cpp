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
   
    double density = 0.1;
    double temperature = 0.072191;
    double pressure = 0.0054;

    double spring_constant = 15000.;
    double restraint_value = 0.1;

    int nstep = 500000;
    
    if (argc >= 3) {
        spring_constant = atof(argv[1]);
        restraint_value = atof(argv[2]);
	assert(spring_constant > 1.0);
	assert(restraint_value > 0.);
    }
    
    if (argc >= 4)
	nstep = atoi(argv[3]);
    assert(nstep > 0);   
 
   
   if (argc >= 7){
	temperature = atof(argv[5]);
        pressure = atof(argv[6]);
	assert(temperature > 0.0);
	assert(pressure > 0.); 
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
    std::cout << "spring_constant\t" << spring_constant << std::endl;
    std::cout << "restraint_value\t" << restraint_value << std::endl;
    
   
    /* Build system with lattice*/
    LatticeFCC lattice;
    lattice.setNumberOfLatticePoints(nmol);
    lattice.setDensity(1.1*density);
    lattice.generateLattice();
    System.buildMoleculeListAndBox(lattice,box);
    
    
    
    /* Build NVT system */
    SimulationNVTEnsemble simulation(System, box, rng);
    simulation.setCosAngleMax(0.825);
    simulation.setDensity(density);
    
    //update neighbor list info
    neighbor_list_info n_info;
    n_info.cut_off_sqd = 5.0;
    
    //set cut_off_sqd
    if (argc >= 5)
	n_info.cut_off_sqd = atof(argv[4]);
    assert(n_info.cut_off_sqd > 2.5);

    simulation.setNeighborInfo(n_info);
        
    UmbrellaSpring* density_bias = new UmbrellaSpring();
    density_bias->order_parameter = restraint_value;
    density_bias->spring_constant = spring_constant;
    
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
        
        
        
        if (i==0) continue;
        simulation.run(50000);

        std::cout << "Production..." << std::endl;
	simulation.resetSteps();
        simulation.setEquilibrate(false);
        simulation.openFile();
        for (int i=0; i<nstep; i++){
           // if (i%(10*nmol)==0) simulation.writeConfig();
            npt.run(100);
    //        if (i%(2*nmol)==0) std::cout << simulation.getDensity() << std::endl;
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
    
    return 0;
   

}
