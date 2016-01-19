//
//  radial_distribution_function.hpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 12/21/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#ifndef radial_distribution_function_hpp
#define radial_distribution_function_hpp

#include <stdio.h>
#include "structure_analysis.hpp"
#include "histogram_dynamic.h"
#include <fstream>
#include <memory>

namespace TetrahedralParticlesInConfinement {
    class RadialDistributionFunction: public StructureAnalysis{
    public:
        RadialDistributionFunction(SimulationNVTEnsemble*, double binsize=0.05);
        
        void setComputeCoreCoreFlag(bool);
        void setComputeCorePatchFlag(bool);
        void setComputePatchPatchFlag(bool);
        void setBinSize(double);
        void setPairInfo(pair_info&);
        
        void compute();
        void print();
        void update();
        
    protected:
        bool_list_t _rdf_flags;
        struct _normalize_info{
            _normalize_info():
            _same_type(false),
            _ntypei(0.),
            _ntypej(0.){
            }
            bool _same_type;
            double _ntypei;
            double _ntypej;
        };
        
        unsigned long _number_of_frames;
        double _bin_size;
        double _rmax;
        double _vol;
        
        std::vector<function1d_t> _GofRs;
        std::vector< utils::HistogramDynamic < double > > _Hists;
        std::vector< _normalize_info > _hist_infos;
        pair_info _info;
        
	std::vector< std::shared_ptr<std::ofstream> > _ofiles;       
 
        void _updateGofRcorecore();
        void _updateGofRcorepatch();
        void _updateGofRpatchpatch();
        void _updateGofRpatchpatchbound();
	void _openFiles();
	void _closeFiles();
        
        void _initialize();
        void _normalize();
        
        function1d_t _normalize(const utils::HistogramDynamic<double>&, _normalize_info&);
        
        enum {CORECORE=0, COREPATCH=1, PATCHPATCH=2, PATCHPATCHBOND=3};
        
    };
}

#endif /* radial_distribution_function_hpp */
