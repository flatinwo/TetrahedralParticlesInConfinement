//
//  bond_structure_analysis.hpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 1/18/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef bond_structure_analysis_hpp
#define bond_structure_analysis_hpp

#include <stdio.h>
#include "structure_analysis.hpp"
#include "histogram_dynamic.h"

namespace TetrahedralParticlesInConfinement {
    class BondStructureAnalysis: public StructureAnalysis{
    public:
        BondStructureAnalysis(SimulationNVTEnsemble*, int l=6);
        BondStructureAnalysis(MoleculeList*, Box*, int l=6);
        
        ~BondStructureAnalysis();
        
        enum Calc_t {GLOBAL, LOCAL, GLOBALANDLOCAL};
        
        void setLvalue(int);
        void setRcutOff(double);
        void setMaxNumberOfNearestNeighbors(unsigned int);
        
        double getQl(bool=false);
        double getWl(bool=false);
        
        
        void compute();
        
        void print();
        void update();
        
    protected:
        int _l;                                     //make this a vector also to accommodate several L values
        double _rcutoff;                            //maximum distance for neighbors
        double _rcutoffsqd;                         //maximum distance sqd
        bool _requireThirdOrderInvaraints;          //should I compute third order invariants
        unsigned_list_t _number_of_neighbors;       //count of number of neighbors
        component_list_t _Qlm;                      //make this a vector to accommodate several L values
        component_list_t _Wl_i;                     //ditto
        unsigned int _max_number_of_neighbors;
        bool _useMaxNumberOfNeighbors;
        std::vector<double_unsigned_pair1d_t> _nearest_neighbors;
        bool _computed;
        
        std::vector<component_list_t> _qlm_i;       //make this a vector also to accommodate several LM values
        Calc_t _mode;
        coord_t _coord;
        
        double _Ql,_Wl;
        
        void _computeWithRcutOff();
        void _computeWithMaxNeighbors();            //also considered as maximum number of bonds
        
        void _computeHarmonics(unsigned int, unsigned int);
        double _ThreeJSymbol(int, int, int, int, int, int); //well-defined Wigner-3j symbol
        
        void _computeQl();
        void _computeWl();
        void _computeWl_i(unsigned int);
        void _resize();
        void _computeNearestNeighbors();
        void _refresh();
        void _configure();
        
        coord_list_t* _com;
        
    };
}

#endif /* bond_structure_analysis_hpp */
