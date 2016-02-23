//
//  simulation_npt_ensemble.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/10/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__simulation_npt_ensemble__
#define __TetrahedralParticlesInConfinement__simulation_npt_ensemble__

#include <stdio.h>
#include "struct_def.h"
#include "simulation_nvt_ensemble.h"
#include "bond_structure_analysis.hpp"

namespace TetrahedralParticlesInConfinement{
    class SimulationNPTEnsemble{
        
        friend class AnalyzeSimulationStepSize;
        friend class UmbrellaSimulation;
        
    public:
        SimulationNPTEnsemble(SimulationNVTEnsemble& NVT, RandomNumberGenerator& rng, double pressure = 1.0);
        SimulationNPTEnsemble(SimulationNVTEnsemble& NVT, RandomNumberGenerator& rng, double pressure, UmbrellaSpring&);
        SimulationNPTEnsemble(SimulationNVTEnsemble& NVT, RandomNumberGenerator& rng, double pressure, UmbrellaSpring&, UmbrellaSpring);
        ~SimulationNPTEnsemble();
        
        int vol_move_per_cycle; //this is chosen -1
        
        void setPressure(double);
        void setDensity(double);
        void sample();

        
        double getPressure();
        double getDensity();
        
        void addUmbrellaSpring(UmbrellaSpring&);
        
        
        move_info& getVolumeInfo();
        SimulationNVTEnsemble& getNVTEnsemble();
        
        void run(int nsteps);
        
    protected:
        SimulationNVTEnsemble& _NVT;
        RandomNumberGenerator& _rng;
        move_info _volume_info;
        MoleculeList old_list;
        Box	old_box;
        double _pressure;
        int _update_volume_move_frequency_per_cycle;
        
        UmbrellaSpring* _umbrella;
        UmbrellaSpring* _umbrella_q6;
        
        int _steps;

        int attemptVolumeMove();
        int attemptVolumeMoveOptimized();
        
        std::shared_ptr<BondStructureAnalysis> oldQ;
        std::shared_ptr<BondStructureAnalysis> newQ;
        
        std::ostringstream _os;
        
    };
}



#endif /* defined(__TetrahedralParticlesInConfinement__simulation_npt_ensemble__) */
