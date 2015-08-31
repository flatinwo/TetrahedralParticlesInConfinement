//
//  moves.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__moves__
#define __TetrahedralParticlesInConfinement__moves__

#include <stdio.h>
#include "struct_def.h"
#include "molecule_list.h"

//update moves in simulation

namespace TetrahedralParticlesInConfinement {
    
    typedef std::pair<coord_t, coord_list_t> coord_pair;
    
    struct move_info{
        
        move_info():
        delta_move(0.1),
        delta_move_max(1.0),
        move_probability(1.0),
        accepted_moves(0),
        rejected_moves(0),
        total_moves(0){}
        
        double delta_move; // step size of move
        double delta_move_max;
        double move_probability;
        int accepted_moves;
        int total_moves;
        int rejected_moves;
        
        void compute_move_probability();
        void update_delta_move();//
        void reset();
        double get_move_probability();//
        
        friend std::ostream& operator<< (std::ostream&, const move_info&);
    };
    
    
    void translate(coord_t& x, Box& box, move_info& translate_info); //done... not tested
    void translate(coord_list_t& x, Box& box, move_info& translate_info); //done... not tested
    
    void rotate(coord_t& orientation, move_info& rotate_info); //done... not tested
    void rotate_sites_per_molecule(coord_list_t& orientation_list, move_info& rotate_info);
    coord_t rotate_sites(coord_list_t& orientation_list, move_info& rotate_info);
    
    
    void rotate(coord_list_t& orientation_list, move_info& rotate_info); //done... not tested
    
    void rotate(TetramerPatchyColloid&, move_info& rotate_info);
    void rotate(Colloid&, move_info&);
    void rotate(Colloid&, coord_pair&);
    void rotate(Colloid& colloid1, Colloid& colloid_ref, move_info&);
    void rotate(Colloid& colloid1, Colloid& colloid_ref, double bond_length, move_info&);
    
    void rotate(int, TetramerPatchyColloid&, move_info&);
    void rotate(Colloid&, coord_list_t&, move_info&);
    
    void rescale(coord_t& x, double s);
    void rescale(coord_list_t& x, double s);
    void rescale(Box& box, double s);
    void rescale(MoleculeList& System, double s);
    
    void anisotropic_box_move(coord_list_t& x, Box& box, move_info& box_move_info);
    void isotropic_box_move(coord_list_t& x, Box& box, move_info& box_move_info);
    
    double gasdev(); //old school gaussian number generator... perhaps create a new class for it
    
    coord_pair generateQuaternionPair(move_info& rotate_info); //done... not tested
}

#endif /* defined(__TetrahedralParticlesInConfinement__moves__) */
