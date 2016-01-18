//
//  structure_analysis.hpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 12/21/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef structure_analysis_hpp
#define structure_analysis_hpp

#include <stdio.h>
#include "struct_def.h"
#include "simulation_nvt_ensemble.h"


namespace  TetrahedralParticlesInConfinement{
    class StructureAnalysis{
    public:
        StructureAnalysis(SimulationNVTEnsemble*);
        ~StructureAnalysis();
        
        void virtual compute();
        
    private:
        
        
    protected:
        MoleculeList* _system;
        Box* _box;
        NeighborList_with_info_t* _neighbors;
        unsigned int _nmolecules;
        
    };
}



#endif /* structure_analysis_hpp */
