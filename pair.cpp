//
//  pair.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "pair.h"
#include "spatial.h"
#include "io.h"
#include <iostream>
#include <cassert>

namespace TetrahedralParticlesInConfinement {
    
    double compute_pair_energy(Colloid& colloid1, Colloid& colloid2,
                               Box& box, pair_info& info){
        
        if (colloid1.core || colloid2.core) {
            return 0.;
        }
        
        if (colloid1.molecule_id == colloid2.molecule_id)
            return 0.;
        
        
        //std::cout << colloid1.molecule_id << "\t" << colloid2.molecule_id << std::endl;
        
        coord_t x1 = colloid1._center_of_mass;
        coord_t x2 = colloid2._center_of_mass;
        
        double_coord_t temp = distancesqandvec(x2, x1, box);
        double rsq = temp.first;
        coord_t dx = temp.second;
        
        if (rsq < info.overlap_criteria) {
            info.overlap = true;
            return ((double) BIG_NUM);
        }
        else if (info.cut_off_criteria <= rsq){
            info.overlap = false;
            return 0.;
        }
        else{
            info.overlap = false;
            return (-1.*info.energy_scale*compute_orientations(colloid1,colloid2,temp,info));
        }
    }
    
    
    
    double compute_pair_energy(TetramerPatchyColloid& molecule1, TetramerPatchyColloid& molecule2, Box& box, pair_info& info){
        
        double energy = 0.;
        for (unsigned int i=1; i<molecule1.colloid_list.size(); i++) {
            for (unsigned int j=1; j<molecule2.colloid_list.size(); j++) {
                energy += compute_pair_energy(molecule1.colloid_list[i], molecule2.colloid_list[j], box, info);
            }
        }
        return energy;
        
    }
    
    
    double compute_orientation(Colloid& colloid1, Colloid& colloid2, double_coord_t& pair_data, pair_info& info){
        
        double alpha0,beta0;
        
        coord_t r1 = colloid1.orientation;
        coord_t r2 = colloid2.orientation;
        
        double rsq = pair_data.first;
        coord_t dx = pair_data.second;
        
        alpha0 = beta0 = 0.;
        
        for (unsigned int i=0; i<r1.size(); i++) {
            alpha0 += r1[i]*dx[i];
            beta0  -= r2[i]*dx[i];
        }
        
        double cut_off_criteria = sqrt(rsq)*info.cut_off_orientation;
        
        //std::cout << "alpha and beta and cut-off-criteria\t" << alpha0 << "\t" << beta0 << "\t" <<cut_off_criteria << "\n";
        
        return evaluate_orientation(alpha0, beta0, cut_off_criteria);
        
    }
    
    
    
    //check to see if necessary to have more than one patch calculated
    double compute_orientations(Colloid& colloid1, Colloid& colloid2, double_coord_t& pair_data, pair_info& info){
        
        double alpha0,beta0;
        
        coord_t r1 = colloid1.orientation;
        coord_t r2 = colloid2.orientation;
        
        double rsq = pair_data.first;
        coord_t dx = pair_data.second;
        
        double cut_off_criteria = sqrt(rsq)*info.cut_off_orientation;
        
        double sum_orientation = 0.;
        alpha0 = beta0 = 0.;
        for (unsigned int i=0; i<r1.size(); i++) {
            alpha0 += r1[i]*dx[i];
            beta0  -= r2[i]*dx[i];
        }
        
        sum_orientation += evaluate_orientation(alpha0,beta0,cut_off_criteria);
        return sum_orientation;
        
    }
    
    

    
    double evaluate_orientation(double alpha0, double beta0, double cut_off_criteria){
        
        if ((alpha0 >= cut_off_criteria) && (beta0 >= cut_off_criteria) ) {
            return 1.;
        }
        else{
            return 0.;
        }
    }
    
    
    double compute_pair_energy(MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info){
        
        if (!neighbor_info.built || neighbor_info.full_neighbor_list) {//or rebuild frequency
            neighbor_info.full_neighbor_list = false;
            coord_list_t x = molecule_list.getMoleculeListCoord();
            build_neighbor_list(x, box, neighbor_list, neighbor_info);
        }
        
        int n = (int) molecule_list.molecule_list.size();
        double e = 0.;
        
        for (unsigned int i=0; i<neighbor_list.size(); i++) {
            for (unsigned int j=0; j<neighbor_list[i].size(); j++) {
                e += compute_pair_energy(molecule_list.molecule_list[i],molecule_list.molecule_list[neighbor_list[i][j]],box, pair_info);
            }
            if (pair_info.overlap) {
                return ((double) n*BIG_NUM);
            }
        }
        //neighbor_info.built = false;
        return e;
        
    }
    
    double compute_pair_energy(int index, MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info){
        
        if (!neighbor_info.built || !neighbor_info.full_neighbor_list) {//or rebuild frequency
            neighbor_info.full_neighbor_list = true;
            coord_list_t x = molecule_list.getMoleculeListCoord();
            build_neighbor_list(x, box, neighbor_list, neighbor_info);
        }
        
        int n = (int) molecule_list.molecule_list.size();
        double e = 0.;
        
        for (unsigned int j=0; j<neighbor_list[index].size(); j++) {
            int nbr = neighbor_list[index][j];
            e += compute_pair_energy(molecule_list.molecule_list[index],molecule_list.molecule_list[nbr],
                                     box, pair_info);
            if (pair_info.overlap) {
                return ((double) n*BIG_NUM);
            }
        }
        //neighbor_info.built = false;
        return e;
        
    }
    
    
    double compute_pair_energy_full(MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info){
        
        if (!neighbor_info.built || neighbor_info.full_neighbor_list) {//or rebuild frequency
            neighbor_info.full_neighbor_list = false;
            coord_list_t x = molecule_list.getFullColloidListCoord();
            build_neighbor_list(x, box, neighbor_list, neighbor_info);
        }
        
        int n = (int) molecule_list.molecule_list.size();
        double e = 0.;
        
        for (unsigned int i=0; i<neighbor_list.size(); i++) {
            for (unsigned int j=0; j<neighbor_list[i].size(); j++) {
                
                e += compute_pair_energy(*(molecule_list.full_colloid_list[i]),*(molecule_list.full_colloid_list[neighbor_list[i][j]]),box, pair_info);
                
                if (pair_info.overlap) {
                    
                    return ((double) n*BIG_NUM);
                }
            }
        }
        //neighbor_info.built = false;
        
        return e;
        
    }
    
    double compute_pair_energy_full(int index, MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info,
                               NeighborList_t& neighbor_list,
                               neighbor_list_info& neighbor_info){
        
        if (!neighbor_info.built || !neighbor_info.full_neighbor_list) {//or rebuild frequency
            neighbor_info.full_neighbor_list = true;
            coord_list_t x = molecule_list.getFullColloidListCoord();
            build_neighbor_list(x, box, neighbor_list, neighbor_info);
        }
        
        int n = (int) neighbor_list.size();
        assert(index < n);
        double e = 0.;
        
        for (unsigned int j=0; j<neighbor_list[index].size(); j++) {
            int nbr = neighbor_list[index][j];
            e += compute_pair_energy(*(molecule_list.full_colloid_list[index]),*(molecule_list.full_colloid_list[nbr]),
                                     box, pair_info);
            if (pair_info.overlap) {
                return ((double) n*BIG_NUM);
            }
        }
        //neighbor_info.built = false;
        return e;
        
    }
    
    double compute_pair_molecule_energy_full(int index0,
                                    MoleculeList& molecule_list,
                                    Box& box,
                                    pair_info& pair_info,
                                    NeighborList_t& neighbor_list,
                                    neighbor_list_info& neighbor_info){
        
        if (!neighbor_info.built || !neighbor_info.full_neighbor_list) {//or rebuild frequency
            neighbor_info.full_neighbor_list = true;
            coord_list_t x = molecule_list.getFullColloidListCoord();
            build_neighbor_list(x, box, neighbor_list, neighbor_info);
        }
        
        int n = (int) neighbor_list.size();
        double e = 0.;
        int dim = (int) molecule_list.molecule_list[index0].colloid_list.size();

        for (unsigned int i=1; i < dim; i++){
            int index = dim*index0 + i;
            for (unsigned int j=0; j<neighbor_list[index].size(); j++) {
                int nbr = neighbor_list[index][j];
                e += compute_pair_energy(*(molecule_list.full_colloid_list[index]),*(molecule_list.full_colloid_list[nbr]),box, pair_info);
                if (pair_info.overlap) {
                    return ((double) n*BIG_NUM);
                }
            }
        }
        

        //neighbor_info.built = false;
        return e;
        
    }
    
    
#pragma COMPUTE ENERGY WITHOUT NEIGHBOR LISTS
    
    //energy between particle i and system
    double compute_pair_energy(int index, MoleculeList& molecule_list,
                               Box& box, pair_info& pair_info){
        
        double e = 0.;
        
        for (unsigned int j=0; j<molecule_list.full_colloid_list.size(); j++) {
            if (j!=index)
                e += compute_pair_energy(*(molecule_list.full_colloid_list[index]),*(molecule_list.full_colloid_list[j]),
                                         box, pair_info);
            if (pair_info.overlap) {
                return ((double) BIG_NUM);
            }
        }
        return e;
        
    }
    
    //energy between molecule i and system
    double compute_pair_molecule_energy(int index0,
                                        MoleculeList& molecule_list,
                                        Box& box,
                                        pair_info& pair_info){
        
        double e = 0.;
        int dim = (int) molecule_list.molecule_list[index0].colloid_list.size();
        
        for (unsigned int i=1; i < dim; i++){
            int index = dim*index0 + i;
            for (unsigned int j=0; j<molecule_list.full_colloid_list.size(); j++) {
                if (j!=index)
                    e += compute_pair_energy(*(molecule_list.full_colloid_list[index]),*(molecule_list.full_colloid_list[j]),
                                             box, pair_info);
                if (pair_info.overlap) {
                    return ((double) BIG_NUM);
                }
            }
        }
        return e;
    }
    
    //total system energy
    double compute_pair_energy(MoleculeList& molecule_list,
                                Box& box, pair_info& pair_info){
        
        double e = 0.;
        
        for (unsigned int i=0; i<molecule_list.full_colloid_list.size(); i++) {
            for (unsigned int j=i+1; j<molecule_list.full_colloid_list.size(); j++) {
                
                e += compute_pair_energy(*(molecule_list.full_colloid_list[i]),*(molecule_list.full_colloid_list[j]),box, pair_info);
                
                if (pair_info.overlap) {
                    
                    return ((double) BIG_NUM);
                }
            }
        }
        
        return e;
        
    }
    
    
}

namespace TetrahedralParticlesInConfinement {
    
    //should one do a move rescale
    /*
     double compute_pair_energy(std::vector<TetramerPatchyColloid>& molecule_list, Box& box, pair_info& pair_info,
     CellList& cell_list, cell_list_info& cell_info){
     
     if (!cell_info.build || !cell_list.good()) {//or rebuild frequency
     cell_list.clear();
     cell_list.setBox(box.box_lo,box.box_hi);
     double range = sqrt(pair_info.cut_off_criteria);
     cell_list.setInteractionRange(range);
     cell_list.setNeighborStyle(CellList::FULL);
     cell_info.build = true;
     //std::cout << "built\n";
     }
     
     int n = (int) colloid_list.size();
     double e = 0.;
     
     //std::cout << cell_info.build << "\t" << cell_list.good() << "\n";
     
     if (cell_list.good()) {
     for (int i=0; i<n; i++) {
     cell_list.insert(i, colloid_list[i]._center_of_mass);
     }
     cell_list.resetIterator();
     for (; ; ) {
     std::pair<int, int> p = cell_list.nextPair();
     if (cell_list.end()) {
     break;
     }
     e += compute_pair_energy(colloid_list[p.first],colloid_list[p.second],
     box, pair_info);
     if (pair_info.overlap) {
     return ((double) n*BIG_NUM);
     }
     }
     }
     else{
     std::cerr << "ERROR: no cell_list\n";
     return ((double) n*BIG_NUM);
     }
     
     return e;
     }
     
     */
    

    
    /**
     \brief Computes interaction energy between particle with index index and colloid_list
     Note that this function assumes that the cell_list has been built already
     */
    

    
    /*
     double compute_pair_energy(int index, std::vector<PatchyColloid>& colloid_list, Box& box, pair_info& pair_info, CellList& cell_list, cell_list_info& cell_info){
     
     if (!cell_info.build || !cell_list.good()) {//or rebuild frequency
     cell_list.clear();
     cell_list.setBox(box.box_lo,box.box_hi);
     double range = sqrt(pair_info.cut_off_criteria);
     cell_list.setInteractionRange(range);
     cell_list.setNeighborStyle(CellList::FULL);
     cell_info.build = true;
     }
     
     //std::cout << cell_info.build << "\t" << cell_list.good() << "\n";
     
     std::vector<int> nbr;
     cell_list.getNeighborsOf(index, nbr);
     double e = 0;
     if (cell_list.good()) {
     for (unsigned int j=0; j<nbr.size(); j++) {
     e += compute_pair_energy(colloid_list[index], colloid_list[nbr[j]], box, pair_info);
     
     if (pair_info.overlap) {
     return ((double) BIG_NUM);
     }
     }
     }
     else{
     std::cerr << "ERROR: no cell_list\n";
     return ((double) BIG_NUM);
     }
     
     return e;
     }
     */

}

