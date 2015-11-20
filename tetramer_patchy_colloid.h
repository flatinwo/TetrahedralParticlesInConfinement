//
//  tetramer_patchy_colloid.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__tetramer_patchy_colloid__
#define __TetrahedralParticlesInConfinement__tetramer_patchy_colloid__

#include <stdio.h>
#include "struct_def.h"
#include "constants.h"
#include <iostream>

namespace TetrahedralParticlesInConfinement {
    class TetramerPatchyColloid{
    public:
        //TetramerPatchyColloid();
        TetramerPatchyColloid(double bond_length=1.0);
        
        ~TetramerPatchyColloid();
        
        std::vector<Colloid> colloid_list; // this is the list of colloids
        coord_list_t orientation_list;      // orientation list .. this should be perhaps be a reference
                                            // a reference for the patches on the core
        
        void setBondLength(double);
        void setPatchAngle(double);
        void setColloidDiameter(double);
        void setVariableAngle(double);
        
        void setCenterOfMass(const coord_t&);  //need to use pbc here
        void setCenterOfMass(const coord_t&, const Box&);
        void setCenterOfMasses(const coord_list_t&, const coord_list_t&);
        
        void setOrientationRandomly();
        void setMoleculeID(int index);
        
        friend std::ostream& operator << (std::ostream&, const TetramerPatchyColloid&);
        
    protected:
        int _number_of_particles; //alternative use const here
        tetrahedron_t Template; //alternative use const here
        
        //note: if i want use const then i need to define my own assignment operator
        tetrahedron_t variableTemplate;
        
        //coord_t* _center_of_mass;
        double _bond_length;
        double _theta;
        double _phi;
        double _sigma;
        
        void buildMolecule();
        void rebuildMolecule();
        
        
    };
}

#endif /* defined(__TetrahedralParticlesInConfinement__tetramer_patchy_colloid__) */
