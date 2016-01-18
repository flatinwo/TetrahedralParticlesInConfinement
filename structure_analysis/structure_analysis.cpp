//
//  structure_analysis.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 12/21/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "structure_analysis.hpp"

namespace TetrahedralParticlesInConfinement {
    StructureAnalysis::StructureAnalysis(SimulationNVTEnsemble* nvt):_system(&(nvt->_molecule_list)),_box(&(nvt->_box)){
        _nmolecules = (unsigned int) _system->molecule_list.size();
        
    }
    
    StructureAnalysis::~StructureAnalysis(){
        
    }
    
    void StructureAnalysis::compute(){
        
    }
}
