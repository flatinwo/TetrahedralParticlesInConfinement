//
//  radial_distribution_function.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 12/21/15.
//  Copyright © 2015 Folarin Latinwo. All rights reserved.
//

#include "radial_distribution_function.hpp"
#include "spatial.h"
#include "io.h"
#include "operator.h"

namespace TetrahedralParticlesInConfinement {

#define RDFSMALL 0.00001

    RadialDistributionFunction::RadialDistributionFunction(SimulationNVTEnsemble* nvt, double bin_size):
    StructureAnalysis(nvt),
    _rdf_flags(bool_list_t(4,true)),
    _number_of_frames(0),
    _bin_size(bin_size),
    _rmax(-1),
    _vol(100000),
    _GofRs(4),
    _Hists(4),
    _hist_infos(4)
    {
	assert(_bin_size > RDFSMALL);
        _initialize();
    }
    
    
    void RadialDistributionFunction::setPairInfo(pair_info& info){
        _info = info;
        _info.cut_off_criteria = _rmax;
    }
    
    void RadialDistributionFunction::compute(){
        unsigned long natoms = _system->full_colloid_list.size();
        
        for (unsigned int i=0; i<natoms; i++) {
            for (unsigned int j=i+1; j<natoms; j++) {
                if (_system->full_colloid_list[i]->molecule_id == _system->full_colloid_list[j]->molecule_id) continue;
                
                double_coord_t temp = distancesqandvec(_system->full_colloid_list[j]->_center_of_mass,
                                        _system->full_colloid_list[i]->_center_of_mass,
                                        *_box);
                double rsq = temp.first;
                
                
                if (rsq < _rmax) {
                    double r = sqrt(rsq);
                    if (_rdf_flags[CORECORE]) {
                        if (_system->full_colloid_list[i]->core && _system->full_colloid_list[j]->core) _Hists[CORECORE].insert(r, 2.0);
                    }
                    if (_rdf_flags[COREPATCH]) {
                        if (_system->full_colloid_list[i]->core && !_system->full_colloid_list[j]->core) _Hists[COREPATCH].insert(r, 2.0);
                    }
                    if (_rdf_flags[PATCHPATCH]) {
                        if (!_system->full_colloid_list[i]->core && !_system->full_colloid_list[j]->core) _Hists[PATCHPATCH].insert(r,2.0);
                    }
                    if (_rdf_flags[PATCHPATCHBOND]) {
                        //if (_system->full_colloid_list[i]->boundto == j) _Hists[PATCHPATCHBOND].insert(r,2.0);
                        if (test_orientations(_system->full_colloid_list[i]->orientation, _system->full_colloid_list[j]->orientation, temp, _info)) _Hists[PATCHPATCHBOND].insert(r,2.0);
                        //if (dot_product(_system->full_colloid_list[i]->_center_of_mass, _system->full_colloid_list[j]->_center_of_mass) > 0.95) _Hists[PATCHPATCHBOND].insert(r,2.0);
                    }
                }
            }
        }
        _number_of_frames++;
    }
    
    void RadialDistributionFunction::_initialize(){
        //set up dimension to be square of (half of box period in smallest dimension)
        _rmax = *(std::min_element(_box->box_period.begin(),_box->box_period.end()));
        _rmax *= 0.5;
        _rmax *= _rmax;
        
        
        //set bin size
        for (unsigned int i=0; i<4; i++) {
            _Hists[i].setBinSize(_bin_size);
            _Hists[i].insert(0,0);
            if (i == CORECORE || i == PATCHPATCH || i==PATCHPATCHBOND) _hist_infos[i]._same_type = true;
            
            if (i == CORECORE) {
                _hist_infos[i]._ntypei = (double) _system->molecule_list.size();
                _hist_infos[i]._ntypej = (double) (_system->molecule_list.size() - 1);
            }
            else if (i == COREPATCH){
                _hist_infos[i]._ntypei = (double) _system->molecule_list.size();
                _hist_infos[i]._ntypej = (double) (_system->molecule_list.size() * 4);
            }
            else if (i == PATCHPATCH){
                _hist_infos[i]._ntypei = (double) (_system->molecule_list.size() * 4);
                _hist_infos[i]._ntypej = (double) (_system->molecule_list.size() * 4 - 1);
            }
            else if (i == PATCHPATCHBOND) _hist_infos[i] = _hist_infos[PATCHPATCH];
            else{
                std::cerr << "Error in RadialDistributionFunction.hpp:\nHistograms out of range\n";
                exit(-1);
            }
        }
        
    }
    
    void RadialDistributionFunction::_normalize(){
        assert(_system->molecule_list.size() > 0);
        _vol = _box->getVolume();
        
        for (unsigned int i=0; i<4; i++) _GofRs[i] = _normalize(_Hists[i], _hist_infos[i]);
    }
    
    function1d_t RadialDistributionFunction::_normalize(const utils::HistogramDynamic<double>& _hist, _normalize_info& info){
        function1d_t _gofr, normalization, probability;
        utils::HistogramDynamic<double> temp(_hist);
        
        double density = info._ntypei/_vol;
        double delta = 0.5*temp.getBinSize();
        
        for (unsigned int i=0; i<temp.getNumberOfBins(); i++) {
            std::pair<double,double> bin = temp.getBin(i);
            double r = bin.first;
            double count = bin.second;
            
            double vol_shell = 4.0/3.0*M_PI*(pow(r+delta,3) - pow(r-delta,3));
            double np_ideal_gas = density*vol_shell;
            
            if (normalization.find(r) == normalization.end()) {
                normalization[r] = 0.0;
                probability[r] = 0.0;
            }
            
            normalization[r] += np_ideal_gas*info._ntypej;
            probability[r] += count;
        }
        
        for (function1d_t::iterator i=normalization.begin(); i!=normalization.end(); ++i)
            _gofr[i->first] = probability[i->first]/(i->second*_number_of_frames);
        return _gofr;
        
    }
    
    void RadialDistributionFunction::update(){
        compute();
    }
    
    void RadialDistributionFunction::print(){
        _normalize();
        std::cout << "After " << _number_of_frames << " frames the rdfs are in gr files\n";
	_openFiles();
        for (unsigned int i=0; i<4; i++) *(_ofiles[i]) << _GofRs[i];
	_closeFiles();
    }

    void RadialDistributionFunction::_openFiles(){
	_ofiles.resize(4);
	_ofiles[0].reset(new std::ofstream("grcorecore.txt"));
	_ofiles[1].reset(new std::ofstream("grcorepatch.txt"));
	_ofiles[2].reset(new std::ofstream("grpatchpatch.txt"));
	_ofiles[3].reset(new std::ofstream("grpatchpatchbonded.txt"));	
    }

    void RadialDistributionFunction::_closeFiles(){ 
        for (unsigned int i=0; i<4; i++) _ofiles[i]->close();
    }


    
}

