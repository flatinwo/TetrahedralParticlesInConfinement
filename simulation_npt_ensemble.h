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

namespace TetrahedralParticlesInConfinement{
    class SimulationNPTEnsemble{
        
        friend class AnalyzeSimulationStepSize;
        
    public:
        SimulationNPTEnsemble(SimulationNVTEnsemble& NVT, RandomNumberGenerator& rng, double pressure = 1.0);
        ~SimulationNPTEnsemble();
        
        int vol_move_per_cycle; //this is chosen -1
        
        void setPressure(double);
        void setDensity(double);
        
        double getPressure();
        double getDensity();
        move_info getVolumeInfo();
        SimulationNVTEnsemble& getNVTEnsemble();
        
        void run(int nsteps);
        
    protected:
        double _pressure;
        int _update_volume_move_frequency_per_cycle;
        int _steps;
        
        SimulationNVTEnsemble& _NVT;
        RandomNumberGenerator& _rng;
        move_info _volume_info;
        
        int attemptVolumeMove();
    };
}



#endif /* defined(__TetrahedralParticlesInConfinement__simulation_npt_ensemble__) */
