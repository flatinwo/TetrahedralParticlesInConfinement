//
//  pair.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__pair__
#define __TetrahedralParticlesInConfinement__pair__

#include <stdio.h>
#include "struct_def.h"
#include "molecule_list.h"
#include "make_neighbor_list.h"

//note: can write such that you store bound or unbound with each pair
//use a bool flag

namespace TetrahedralParticlesInConfinement {
    
    struct pair_info{
        pair_info():
        overlap(false),
        overlap_probability(0.),
        energy_scale(1.0),
        overlap_criteria(1.0),
        cut_off_criteria(1.251*1.251),
        cut_off_orientation(0.95)
        {}
        
        bool overlap;
        double overlap_probability;
        double energy_scale;
        double overlap_criteria;
        double cut_off_criteria;
        double cut_off_orientation;
        
    };
    
    
    struct cell_list_info{
        cell_list_info():
        build(false),
        last_build(0),
        build_frequency(10)
        {}
        
        bool build;
        int last_build;
        int build_frequency;
    };
    
    double compute_pair_energy(Colloid& colloid1, Colloid& colloid2, Box& box, pair_info& info); //single patch colloid calculation
    
    double compute_pair_energy(TetramerPatchyColloid& molecule1, TetramerPatchyColloid& molecule2, Box& box, pair_info& info);
    
    
    
    double compute_pair_energy(MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info);//make template for file, question for Ladapo vector of references or reference of vectors
    
    
    
    
    //can use stl pair to combine NeighborList_t and neighbor_list_info
    double compute_pair_energy(int index, MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info);
    
    
    double compute_pair_energy_full(MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info);//make template for file, question for Ladapo vector of references or reference of vectors
    
    
    
    
    //can use stl pair to combine NeighborList_t and neighbor_list_info
    double compute_pair_energy_full(int index, MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info);
    
    double compute_pair_molecule_energy_full(int index,
                                    MoleculeList& molecule_list,
                                    Box& box,
                                    pair_info& pair_info,
                                    NeighborList_t& neighbor_list,
                                    neighbor_list_info& neighbor_info);
    
    
    // energy computes without neighbor lists
    double compute_pair_energy(MoleculeList& molecule_list,
                               Box& box,
                               pair_info& pair_info);//make template for file, question for Ladapo vector of references or reference of vectors
    
    
    double compute_pair_energy(int index,
                               MoleculeList& molecule_list,
                               Box& box,
                               pair_info& pair_info);
    
    double compute_pair_molecule_energy(int index,
                                        MoleculeList& molecule_list,
                                        Box& box,
                                        pair_info& pair_info);
    
    
    
    double compute_orientation(Colloid& colloid1, Colloid& colloid2, double_coord_t& dx, pair_info& info);
    
    double compute_orientations(Colloid& colloid1, Colloid& colloid2, double_coord_t& dx, pair_info& info);
    
    double evaluate_orientation(double alpha0, double beta0, double cut_off_criteria);
    
    enum {BIG_NUM=750};
}


//double compute_pair_energy(std::vector<TetramerPatchyColloid>& molecule_list, Box& box, pair_info& pair_info, CellList& cell_list, cell_list_info& cell_info);//make template for file, question for Ladapo vector of references or reference of vectors

//double compute_pair_energy(int index, std::vector<TetramerPatchyColloid>& colloid_list, Box& box, pair_info& pair_info, CellList& cell_list, cell_list_info& cell_info);
#endif /* defined(__TetrahedralParticlesInConfinement__pair__) */
