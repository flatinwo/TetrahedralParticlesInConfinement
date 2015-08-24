//
//  make_neighbor_list.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__make_neighbor_list__
#define __TetrahedralParticlesInConfinement__make_neighbor_list__

#include <stdio.h>
#include "struct_def.h"

namespace TetrahedralParticlesInConfinement {
    /**
     \brief quick and dirty Neighbor list to be deprecated
     This neighbor list is written in procedural
     programming style
     */
    
    struct neighbor_list_info{
        neighbor_list_info():
        cut_off_sqd(9.00),
        number_of_pairs(0),
        built(false),
        full_neighbor_list(true)
        {}
        
        double cut_off_sqd;
        bool built, sorted;
        int number_of_pairs;
        bool full_neighbor_list;
    };
    
    typedef std::vector<int> int_list_t;
    typedef std::vector< int_list_t > int_2Dlist_t;
    typedef int_2Dlist_t NeighborList_t;
    
    typedef std::pair<NeighborList_t,neighbor_list_info> NeighborList_with_info_t;
    
    void build_neighbor_list(coord_list_t& x, Box& box, int_2Dlist_t& neighbor_list,
                             neighbor_list_info& info);
    void build_neighbor_list(coord_list_t& x, Box&, NeighborList_with_info_t&);
    
    void build_sorted_neighbor_list(coord_list_t& x, Box& box, int_2Dlist_t& neighbor_list,
                                    neighbor_list_info& info);
    void build_sorted_neighbor_list(coord_list_t& x, Box&, NeighborList_with_info_t&);
}

#endif /* defined(__TetrahedralParticlesInConfinement__make_neighbor_list__) */
