//
//  simulation.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__simulation__
#define __TetrahedralParticlesInConfinement__simulation__

#include <stdio.h>
#include "struct_def.h"
#include "molecule_list.h"
#include "make_neighbor_list.h"
#include "random_number_generator.h"
#include "pair.h"
#include "io.h"
#include "cell_list.hpp"
#include <fstream>

namespace TetrahedralParticlesInConfinement {
    
    class AnalyzeSimulationStepSize; //forward declaration for friendships
    
    class Simulation{
    public:
        Simulation(MoleculeList&, Box&, RandomNumberGenerator&);
        ~Simulation();
        
        virtual void run(int);
        
        void setMoleculeCoord(TetramerPatchyColloid&, int index);
        void setMoleculeList(MoleculeList&);
        void setEquilibrate(bool);
        
        void setPairInfo(pair_info& info);
        void setNeighborInfo(neighbor_list_info& info);
        void setCellList(CellList& info);
        void setNMovesPerCycle(double);
        void resetSteps();
        
        
        void addToMoleculeList(const coord_t&);
        void addMolecule();
        void addPlates(Plates& plates, bool rebuild=true);
        void addMoleculeWithoutOverlaps();
        
        void removeMolecule(int);
        void removeMolecule();
        
        double getAcceptanceProbability();
        double getNMovesPerCycle();
        
        MoleculeList& getMoleculeList();
        Box& getBox();
        pair_info& getPairInfo();
        
        void buildNeighborList();
        void buildCellList();
        
        virtual double computeEnergy(int);
        virtual double computeEnergy();
        
        friend std::ostream& operator << (std::ostream&, const Simulation&);
        
    protected:
        int _nmovespercycle;
        int _steps;
        double _max_displacement;
        double _delta_skin;
        
        bool _equilibrate;
        int _update_neighbors_frequency_per_cycle;
        
        MoleculeList& _molecule_list;
        Box&         _box;
        TetramerPatchyColloid _old_colloid_coords;
        
        NeighborList_with_info_t _neighbor_list;
        CellList*    _cell_list;
        
        pair_info _pair_info;
        RandomNumberGenerator& _rng;
        coord_list_t _coords_since_last_neighbor_build;
        Plates* confinement;
        
        //build two neighborlists, use in umbrella
        std::ofstream _ofile;
        std::ofstream _ofile_energy;
        
    };
}

#endif /* defined(__TetrahedralParticlesInConfinement__simulation__) */
