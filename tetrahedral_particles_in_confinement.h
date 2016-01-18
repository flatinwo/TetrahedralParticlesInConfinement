//
//  tetrahedral_particles_in_confinement.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef TetrahedralParticlesInConfinement_tetrahedral_particles_in_confinement_h
#define TetrahedralParticlesInConfinement_tetrahedral_particles_in_confinement_h

#include "struct_def.h"
#include "lattice_fcc.h"
#include "tetramer_patchy_colloid.h"
#include "constants.h"
#include "io.h"
#include "spatial.h"
#include "molecule_list.h"
#include "pair.h"
#include "make_neighbor_list.h"
#include "simulation_nvt_ensemble.h"
#include "moves.h"
#include "operator.h"
#include "preprocess.h"
#include "random_number_generator_48.h"
#include "simulation_npt_ensemble.h"
#include "umbrella_simulation.h"
#include "analyze_simulation_step_size.h"

//how about multistate runs

namespace TetrahedralParticlesInConfinement {
    enum SimulationMode {AUTOSTART, AUTORESTART, MANUAL};
    
    class TPC{
    public:
        TPC();
        TPC(std::vector<std::string>&);
        
        void setSimulationMode(SimulationMode);
        void run();
        void setup();
        void reset();
        
    protected:
        bool _equilibrateFlag;
        bool _analysisFlag;
        bool _pressureFlag;
        
        int _nMolecules;
        unsigned int _ncyclesEquilibration,_ncyclesProduction;
        unsigned int _configOutputFrequency;
        
        double _temperature, _density, _pressure;
        double _bondLength, _cosAngleMax;
        double _plate_separation;
        
        MoleculeList _system;
        RandomNumberGenerator48 _rng,_rng2,_rng3;
        Box _box;
        Plates _plates;
        SimulationMode _mode;
        
        std::unique_ptr<UmbrellaSpring> _density_bias; //maybe make a vector 
        std::unique_ptr<SimulationNVTEnsemble> _nvt;
        std::unique_ptr<SimulationNPTEnsemble> _npt;
        SimulationGibbsEnsemble* _gibbs;
        
        void _printInitialObjectConfig();
        void _printFinalObjectConfig();
        void _initialize();
        
        std::ofstream _ofile;
        
        
    };
}

#endif
