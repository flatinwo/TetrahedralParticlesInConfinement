//
//  umbrella_simulation.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 8/10/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__umbrella_simulation__
#define __TetrahedralParticlesInConfinement__umbrella_simulation__

#include <stdio.h>
#include <stdio.h>
#include "struct_def.h"
#include "simulation_npt_ensemble.h"
#include "bond_structure_analysis.hpp"
#include <memory>

//#include "simulation_gibbs_ensemble.h"

namespace TetrahedralParticlesInConfinement {
    class SimulationGibbsEnsemble;
    
    class UmbrellaSimulation{
    public:
        UmbrellaSimulation(SimulationNPTEnsemble&, RandomNumberGenerator&, UmbrellaSpring*);
        UmbrellaSimulation(SimulationNPTEnsemble&, RandomNumberGenerator&, UmbrellaSpring*, UmbrellaSpring*);
        
        
        void setUmbrellaType(int);
        void setUmbrellaRestraintValue(double);
        void setUmbrellaRestraintSpringConstant(double);
        void setEquilibrate(bool);
        void setSampleFrequency(int);
        void setMCSweepsperUmbrellaSweep(int);
        
        double getRunningAverage();
        
        void openFile();
        void closeFile();
        
        void run(int nstep);
        void run(std::pair<int,int>& step_data);
        void resetRunningAverage();
        void resetSteps();
        
    protected:
        int _nstepsMC;
        unsigned int _nsampleFrequency;
        double _E, _restrain_value, _goodQ6;
        double _beta;
        bool _equilibrate;
        int _umbrella_type;
        
        int _running_count;
        double _running_avg;
        
        move_info _umbrella_info;
        
        
        SimulationNVTEnsemble* _NVT; //use null to switch the flow
        SimulationNPTEnsemble* _NPT; //use null to switch the flow
        SimulationGibbsEnsemble* _Gibbs; //use null to switch the flow
        RandomNumberGenerator* _rng;
        UmbrellaSpring* _umbrella;
        UmbrellaSpring* _umbrella1;
        
        std::unique_ptr<BondStructureAnalysis> _q6analysis;
        
        std::ofstream _ofile;			//< File to write umbrella data to
        
        struct Config_t{
            double _old_density;
            double _Q6;
            double _E;
            coord_list_t _old_coord_list;
            MoleculeList _old_list;
            Box _old_box;
        } _old_configs;
        
        void saveConfig(double); //could use auto or typeid here
        void saveConfig(coord_list_t);
        void saveConfig();
        
        int _steps;
        
        enum {DENSITY=0, Q6=1, DENSITYANDQ6=2};
    };
}


#endif /* defined(__TetrahedralParticlesInConfinement__umbrella_simulation__) */
