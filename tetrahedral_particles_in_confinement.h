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
#include "radial_distribution_function.hpp"
#include "bond_structure_analysis.hpp"
#include <memory>

//how about multistate runs

namespace TetrahedralParticlesInConfinement {
    enum CalculationMode {AUTOSTART, AUTORESTART, MANUAL};
    
    class TPiC{
    public:
        TPiC();
        TPiC(std::vector<std::string>&);
        TPiC(int argc, const char* argv[]);
        TPiC(const char*);
        
        void setSimulationMode(CalculationMode);
        void run();
        void run_umbrella();
        void reset();
        
    protected:
        bool _equilibrateFlag;
        bool _analysisFlag;
        bool _pressureFlag;
        
        int _nMolecules;
        unsigned int _ncyclesEquilibration,_ncyclesProduction, _ncyclesAnalysis;
        unsigned int _configOutputFrequency;
        
        double _temperature, _density, _pressure;
        double _bondLength, _cosAngleMax;
        double _plate_separation;
        
        MoleculeList _system;
        RandomNumberGenerator48 _rng;
        std::vector<RandomNumberGenerator48*> _rngs;
        void _generateRandomSeeds();
        
        
        Box _box;
        CalculationMode _mode;
        std::vector<move_info> _move_info_list;     //to translate, rotate arms, rotate molecules, restart umbrella change spring constant
        
        
        std::vector<std::string> _command_list;
        
        std::unique_ptr<Plates> _plates;
        
        std::shared_ptr<UmbrellaSpring> _density_bias; //maybe make a vector
        std::shared_ptr<UmbrellaSpring> _q6_bias;
        
        std::shared_ptr<SimulationNVTEnsemble> _nvt;
        void _nvt_equilibrate();
        void _nvt_analysis();
        
        
        std::shared_ptr<SimulationNPTEnsemble> _npt;
        void _npt_equilibrate();
        void _npt_analysis();
        
        
        std::shared_ptr<UmbrellaSimulation> _umbrella;
        
        std::unique_ptr<AnalyzeSimulationStepSize> _analysis;
        
        SimulationGibbsEnsemble* _gibbs;
        
        
        void _setup();
        void _nullAllPtrs();
        
        void _printInitialObjectConfig();
        void _printFinalObjectConfig();
        void _initialize();
        
        std::ofstream _ofile;
        
        
        
    };
}

#endif
