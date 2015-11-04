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
//#include "simulation_gibbs_ensemble.h"

namespace TetrahedralParticlesInConfinement {
    class SimulationGibbsEnsemble;
    
    class UmbrellaSimulation{
    public:
        UmbrellaSimulation(SimulationNPTEnsemble&, RandomNumberGenerator&, UmbrellaSpring*);
        
        void setUmbrellaType(int);
        void setUmbrellaRestraintValue(double);
        void setUmbrellaRestraintSpringConstant(double);
        void setEquilibrate(bool);
        
        double getRunningAverage();
        
        void openFile();
        void closeFile();
        
        void run(int nstep);
        void run(std::pair<int,int>& step_data);
        void resetRunningAverage();
        
    protected:
        int _nstepsMC;
        double _E, _restrain_value;
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
        
        std::ofstream _ofile;			//< File to write umbrella data to
        
        struct Config_t{
            double _old_density;
            coord_list_t _old_coord_list;
            MoleculeList _old_list;
            Box _old_box;
        } _old_configs;
        
        void saveConfig(double); //could use auto or typeid here
        void saveConfig(coord_list_t);
        void saveConfig();
        
        int _steps;
        
        enum {DENSITY=0, Q6=1};
    };
}


#endif /* defined(__TetrahedralParticlesInConfinement__umbrella_simulation__) */
