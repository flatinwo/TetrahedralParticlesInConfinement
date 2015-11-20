//
//  make_neighbor_list.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "make_neighbor_list.h"
#include "spatial.h"
#include <cassert>
#include <iostream>

//kd-trees at some point

namespace TetrahedralParticlesInConfinement {
    
    void build_neighbor_list(coord_list_t& x, Box& box, int_2Dlist_t& neighbor_list,
                             neighbor_list_info& info){
        
        neighbor_list.clear();
        info.number_of_pairs = 0;
        
        assert(x.size()  != 0);
        assert(box.box_lo.size() != 0);
        
        int dim = (int) x.size();
        
        //first assumes a HALF neighbor list
        
        double rmagsqd = 0;
        
        for (int i=0; i<dim; i++) {
            int_list_t neighbor_i;
            for (int j=i+1; j<dim; j++) {
                rmagsqd = distancesq(x[i], x[j], box);
                if (rmagsqd < info.cut_off_sqd) {
                    neighbor_i.push_back(j);
                    info.number_of_pairs++;
                }
            }
            neighbor_list.push_back(neighbor_i);
        }
        
        //build full neighbor list
        //maybe use insert here
        //info.full_neighbor_list = false;
        
        
        if (info.full_neighbor_list) {
            int_list_t index;
            for (unsigned int i=0; i<neighbor_list.size(); i++) {
                index.push_back((int) neighbor_list[i].size());
            }
            
            for (unsigned int i=0; i<neighbor_list.size(); i++) {
                for (int j=0; j<index[i]; j++) {
                    int index_j = neighbor_list[i][j];
                    neighbor_list[index_j].push_back(i);
                }
            }
        }
        
        info.built = true;
        
        
    }
    
    void build_neighbor_list(coord_list_t& x, Box& box, NeighborList_with_info_t& n_info){
        build_neighbor_list(x, box, n_info.first, n_info.second);
        
    }
    
    void build_sorted_neighbor_list(coord_list_t& x, Box& box, int_2Dlist_t& neighbor_list,
                                    neighbor_list_info& info);
    
    
    
}