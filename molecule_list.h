//
//  molecule_list.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__molecule_list__
#define __TetrahedralParticlesInConfinement__molecule_list__

#include <stdio.h>
#include "struct_def.h"
#include "tetramer_patchy_colloid.h"

namespace TetrahedralParticlesInConfinement {
    class MoleculeList{
    public:
        MoleculeList();
        ~MoleculeList(); //why does destructor code stop working
        
        const int nparticle_per_molecule;
        
        std::vector<TetramerPatchyColloid> molecule_list;
        std::vector<Colloid*> full_colloid_list;
        
        //make a list of pointers so you can have whatever you want
        //std::vector< PatchyColloid* > pointer_colloid_list;
        //std::vector< coord_t& > _cm_list;
        
        //to use this pointer I need to write my own custom copy constructor
        //and assignment operator
        
        
        void buildMoleculeList(Lattice& lattice_list);
        void buildMoleculeListAndBox(Lattice& lattice_list, Box& box);
        void popBackMolecule();
        
        void addToMoleculeList(TetramerPatchyColloid&);
        void addToMoleculeList(const coord_t& center_of_mass);
        void addToMoleculeList(const coord_list_t& center_of_mass, const coord_list_t& orientation, const Box& box);
        void addToMoleculeList(const coord_t& center_of_mass, const Box& box);
        
        coord_list_t& getMoleculeListCoord();
        coord_list_t& getMoleculeListOrientation();
        coord_list_t& getMoleculeListOrientation(int index);
        double getBondLength();
        
        coord_list_t& getFullColloidListCoord(); 
        
        void setMoleculeListOrientation(coord_list_t&);
        void setMoleculeListOrientation();
        void setMoleculeListDiameter(double);
        void setMoleculeListBondLength(double);
        
        void buildFullColloidListPointer();
        
        bool are_there_overlaps();
        
        friend std::ostream& operator << (std::ostream&, const MoleculeList&);
        
    protected:
        double _bond_length;
        bool _is_good;
        bool _is_bad;
        Box _box;
        
        coord_list_t _center_of_mass_list;
        coord_list_t _center_of_mass_colloid_list;
        coord_list_t _orientation_list;
        
    };
}
#endif /* defined(__TetrahedralParticlesInConfinement__molecule_list__) */
