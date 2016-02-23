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
    
#define MAX_NUMBER_OF_NEIGHBORS 36
    
    BondStructureAnalysis::BondStructureAnalysis(SimulationNVTEnsemble* nvt, int l):
    StructureAnalysis(nvt),
    _l(l),
    _mode(GLOBAL){
        _configure();
    }
    
    BondStructureAnalysis::BondStructureAnalysis(MoleculeList* system, Box* box, int l):
    StructureAnalysis(system,box),
    _l(l),
    _mode(GLOBAL){
        _configure();
    }
    
    BondStructureAnalysis::~BondStructureAnalysis(){
        
    }
    
    void BondStructureAnalysis::_configure(){
        assert(_l>1);
        _rcutoff = 0.74;
        _requireThirdOrderInvaraints = false;
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
        _computed = false;
    }
    
    void BondStructureAnalysis::setLvalue(int l){
        assert(l>1);
        _l = l;
    }
    
    void BondStructureAnalysis::setRcutOff(double rcutoff){
        assert(rcutoff > 0.);
        _rcutoff = rcutoff;
    }
    
    void BondStructureAnalysis::setMaxNumberOfNearestNeighbors(unsigned int n_nghbrs){
        assert(n_nghbrs > 0);
        _max_number_of_neighbors = n_nghbrs;
        _useMaxNumberOfNeighbors = true;
        
    }
    
    double BondStructureAnalysis::getQl(bool compute_flag){
        if (compute_flag) compute();
        return _Ql;
    }
    
    double BondStructureAnalysis::getWl(bool compute_flag){
        if (!_requireThirdOrderInvaraints) std::cerr << "ERROR: BondStructureAnalysis::getWl()\n" <<
            "Cannot get Wl without requiring it\n";
        if (compute_flag) compute();
        return _Wl;
        
    }
    
    void BondStructureAnalysis::_computeNearestNeighbors(){
        double rcutsqd = _rcutoff*_rcutoff;
        
        for (unsigned int i=0; i<_com->size(); i++) {
            unsigned int k=0;
            for (unsigned int j=0; j<_com->size(); j++) {
                if (i==j) continue;
                double rsq = distancesq((*_com)[i], (*_com)[j], *_box);
                if (rsq < rcutsqd) {
                    _nearest_neighbors[i][k].first = rsq;
                    _nearest_neighbors[i][k].second = j;
                    k++;
                    assert(k < MAX_NUMBER_OF_NEIGHBORS);
                }
            }
            std::sort(_nearest_neighbors[i].begin(), _nearest_neighbors[i].begin()+k); //implement sort up to
            assert(k >= _max_number_of_neighbors);
            assert(k > 0);  //Number of nearest numbers too small
        }
    }
    
    void BondStructureAnalysis::compute(){
        if (!_useMaxNumberOfNeighbors) _computeWithRcutOff();
        else _computeWithMaxNeighbors();
        _computed = true;
    }
    
    void BondStructureAnalysis::_computeWithMaxNeighbors(){
        _com = &(_system->getMoleculeListCoord());
        _computeNearestNeighbors(); //get information about nearest neighbors
        for (unsigned int i=0; i<_com->size(); i++) {
            for (unsigned int j=0; j < _max_number_of_neighbors; j++) {
                unsigned int k = _nearest_neighbors[i][j].second;
                _computeHarmonics(i, k);
            }
            for (unsigned int m=0; m <_qlm_i[i].size(); m++) {
                assert(_number_of_neighbors[i] > 0);
                assert(_number_of_neighbors[i] == _max_number_of_neighbors);
                _qlm_i[i][m] /= (double) _number_of_neighbors[i];
                _Qlm[m] += _qlm_i[i][m];
            }
            if (_requireThirdOrderInvaraints) _computeWl_i(i);
        }
        _computeQl();
        if (_requireThirdOrderInvaraints) _computeWl();
    }
    
    void BondStructureAnalysis::_computeWithRcutOff(){
        _com = &(_system->getMoleculeListCoord());
        for (unsigned int i=0; i<_com->size(); i++) {
            for (unsigned int j=0; j<_com->size(); j++) {
                if (i==j) continue;
                _computeHarmonics(i, j);
            }
            for (unsigned int m=0; m <_qlm_i[i].size(); m++) {
                assert(_number_of_neighbors[i] > 0);
                _qlm_i[i][m] /= (double) _number_of_neighbors[i];
                _Qlm[m] += _qlm_i[i][m];
            }
            if (_requireThirdOrderInvaraints) _computeWl_i(i);
        }
        _computeQl();
        if (_requireThirdOrderInvaraints) _computeWl();
        
    }
    
    //compute Harmonics for system
    void BondStructureAnalysis::_computeHarmonics(unsigned int i, unsigned int j){
        
        _coord = spherical_orientation((*_com)[i],(*_com)[j],*_box);
        if (_coord[0] > _rcutoff) return;
        else{
            double theta = _coord[1];
            double phi = _coord[2];
            assert(theta <= M_PI);
            assert(phi < 2*M_PI);
            for (int m = -_l; m < _l+1 ; m++) {
                _qlm_i[i][m + _l] += spherical_harmonic(_l, m, theta, phi);
            }
            _number_of_neighbors[i]++;
        }
    }
    
    double BondStructureAnalysis::_ThreeJSymbol(int ja, int jb, int jc, int ma, int mb, int mc){
        return gsl_sf_coupling_3j(2*ja, 2*jb, 2*jc, 2*ma, 2*mb, 2*mc);
    }
    
    //compute Ql of system
    void BondStructureAnalysis::_computeQl(){
        double sum = 0;
        for (unsigned int i=0; i<_Qlm.size();i++) sum += std::norm(_Qlm[i]);
        sum *= (4.*M_PI/(double) (2*_l + 1));
        _Ql = sqrt(sum)/(double) _nmolecules;
        if (!_requireThirdOrderInvaraints) _refresh();
    }
    
    //compute Wl of system
    void BondStructureAnalysis::_computeWl_i(unsigned int i){
        
        for (int m1=-_l; m1 < _l+1; m1++) {
            for (int m2=-_l; m2 < _l+1; m2++) {
                for (int m3=-_l; m3<_l+1; m3++) {
                    if (m1 + m2 + m3 != 0) continue;
                    double w3j = _ThreeJSymbol(_l, _l, _l, m1, m2, m3);
                    _Wl_i[i] += w3j*(_qlm_i[i][m1+_l]*_qlm_i[i][m2+_l]*_qlm_i[i][m3+_l]);
                }
            }
        }
        double normalization = 0.;
        for (auto& m : _qlm_i[i]) normalization += std::norm(m);
        _Wl_i[i] /= std::pow(normalization, 1.5);
        
    }
    
    void BondStructureAnalysis::_computeWl(){
        _Wl = 0.;
        for (auto& i : _Wl_i) _Wl += i.real();
        _Wl /= (double) _Wl_i.size();
        _refresh();
    }

    
    void BondStructureAnalysis::_refresh(){
        //refresh _qlmi
        for (auto& m : _qlm_i)
            for (auto& l : m) l = std::complex<double>(0.,0.);
        
        //refresh _Qlm
        for (auto& i: _Qlm) i = std::complex<double>(0.,0.);
        
        //number of nearest neighbors
        for (auto& n :_number_of_neighbors) n = 0;
        
        //refresh Wl_i
        if (_requireThirdOrderInvaraints)
            for (auto& i : _Wl_i) i = std::complex<double>(0.,0.);
        
    }
}

