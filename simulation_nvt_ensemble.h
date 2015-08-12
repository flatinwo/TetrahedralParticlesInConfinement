//
//  simulation_nvt_ensemble.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/30/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef __TetrahedralParticlesInConfinement__simulation_nvt_ensemble__
#define __TetrahedralParticlesInConfinement__simulation_nvt_ensemble__

#include <stdio.h>
#include "struct_def.h"
#include "simulation.h"
#include <iomanip>
#include "moves.h"
#include <string>

//question save config with static keyword
//using function pointers or templates, i need to think.. we can compute

namespace TetrahedralParticlesInConfinement{
    class SimulationNVTEnsemble : public Simulation {
        
        friend class AnalyzeSimulationStepSize;
        
    public:
        SimulationNVTEnsemble(MoleculeList&, Box&, RandomNumberGenerator&);
        ~SimulationNVTEnsemble();
        
        void setBeta(double);
        void setDensity(double);
        void setUpdateMoveFrequencyPerCycle(int);
        void setDeltaMoveTranslate(double);
        void setDeltaMoveRotate(double);
        void setCosAngleMax(double);
        void setUpdateNeighborFrequencyPerCycle(int);
        void UpdateNeighborFrequencyPerCycle();
        
        double getTemperature();
        double getDensity();
        double getVolume();
        coord_list_t& getFullColloidListCoord();
        std::map<int,move_info> getMoveInfoMap();
        
        
        void run(int);
        
        //file handling
        void openFile();
        void closeFile();
        void writeConfig();
        
        double computeEnergy(int);
        double computeEnergy();
        double computeMoleculeEnergy(int);
        
        void   computeVolume(); //make return
        double computePressure();
        
        
    protected:
        int _n;
        double _beta;
        double _density;
        double _volume;
        int _nsubmoves;
        int _update_move_frequency_per_cycle;
        double _E, _delE;
        bool _core_flag;
        double _cos_angle_max;
        
        int attemptMove(int);
        int attemptSubMove(int);
        int _flag;
        
        coord_list_t _old_config;
        TetramerPatchyColloid _old_molecule;
        Colloid _old_colloid;
        
        void Translation(int);
        void Rotation(int);
        void TranslationAndRotation(int);
        bool isRotationGood(int);
        
        void    saveConfig(int index);
        void    revertConfig(int index);
        void    updateMoveInfo(int);
        
        std::map<int,move_info> _move_info_map;
        
        enum {TRANSLATE = 0, ROTATE = 1, ROTATEMOLECULE=2};
        enum {REJECT = 0, ACCEPT = 1};
        
        double computePairEnergy(int);
    };
}

#endif /* defined(__TetrahedralParticlesInConfinement__simulation_nvt_ensemble__) */
