//
//  analyze_simulation_step_size.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/12/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__analyze_simulation_step_size__
#define __TetrahedralParticlesInConfinement__analyze_simulation_step_size__

#include <stdio.h>
// #include "simulation_gibbs_ensemble.h"
#include "simulation_npt_ensemble.h"

namespace TetrahedralParticlesInConfinement {
    class AnalyzeSimulationStepSize{
    public:
        AnalyzeSimulationStepSize(SimulationNVTEnsemble&);
        AnalyzeSimulationStepSize(SimulationNPTEnsemble&);
        //AnalyzeSimulationStepSize(SimulationGibbsEnsemble&);
        
        void setUpperBoundProbability(double=0.40);
        void setLowerBoundProbability(double=0.30);
        
        void run(int=1000);
        
    protected:
        //member variables and objects
        int _iteration_count;
        int _count_max;
        
        double _acceptance_probability;
        double _upper_bound;
        double _lower_bound;
        
        
        SimulationNVTEnsemble* _NVT; //use null to switch the flow
        SimulationNPTEnsemble* _NPT; //use null to switch the flow
        SimulationNPTEnsemble* _Gibbs; //use null to switch the flow (note: this should be a gibbs calculation)
        
        std::ofstream _ofile;
        
        //member functions
        void openFile();
        void closeFile();
        void writeFinalResults();
        void updateDeltaMove();
        void updateDeltaMove(double&);
        void updateDeltaMove(move_info&);
        void reset();
        
        void nvt_run(int);
        void npt_run(int);
        void gibbs_run(int);
        
        bool areAllMoveInfoGood();
        
        
        
    };
}

#endif /* defined(__TetrahedralParticlesInConfinement__analyze_simulation_step_size__) */
