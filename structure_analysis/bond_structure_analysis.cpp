//
//  bond_structure_analysis.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 1/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "bond_structure_analysis.hpp"
#include "spatial.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <cassert>
#include <gsl/gsl_sf_coupling.h>
#include <algorithm>

using namespace boost::math;

namespace TetrahedralParticlesInConfinement {
    
#define MAX_NUMBER_OF_NEIGHBORS 16
    
    BondStructureAnalysis::BondStructureAnalysis(SimulationNVTEnsemble* nvt, int l):
    StructureAnalysis(nvt),
    _l(l),
    _mode(GLOBAL){
        assert(_l>1);
        _rcutoff = 0.74;
        _requireThirdOrderInvaraints = true;
        _useMaxNumberOfNeighbors = false;
        _max_number_of_neighbors = 0;
        _nearest_neighbors.resize(_nmolecules,
                                  double_unsigned_pair1d_t (MAX_NUMBER_OF_NEIGHBORS, std::pair<double, unsigned int>(0.,0))) ;
        
        _number_of_neighbors.resize(_nmolecules,0.);
        _Wl = 0.;
        _Wl_i.resize(_nmolecules, std::complex<double>(0.,0.));
        _qlm_i.resize(_nmolecules);
        
        for (auto& m : _qlm_i) m.resize(2*_l+1, std::complex<double>(0,0));
        _Qlm.resize(2*_l+1,std::complex<double>(0,0));
        
    }
    
    void BondStructureAnalysis::compute(){
        
    }
}

