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
#include <cassert>

//question save config with static keyword
//using function pointers or templates, i need to think.. we can compute
//write out production or equilibrating depending on setEquilibrate
//assert to test orientation_list that orientation good

namespace TetrahedralParticlesInConfinement{
    class SimulationNVTEnsemble : public Simulation {
        
        friend class AnalyzeSimulationStepSize;
        friend class TestSimulationNVTEnsemble;
        friend class SimulationNPTEnsemble;
        friend class UmbrellaSimulation;
        
    public:
        SimulationNVTEnsemble(MoleculeList&, Box&, RandomNumberGenerator&);
        ~SimulationNVTEnsemble();
        
        enum PressureMode {COMPRESSION=0, EXPANSION=1, COMPRESSIONANDEXPANSION=2};
        
        void setBeta(double);
        void setDensity(double);
        void setUpdateMoveFrequencyPerCycle(int);
        void setDeltaMoveTranslate(double);
        void setDeltaMoveRotate(double);
        void setCosAngleMax(double);
        void setUpdateNeighborFrequencyPerCycle(int);
        void setPressureCalculation(bool);
        void setPressureCalculationScaleFactor(double);
        void setPressureCalculationMode(PressureMode);
        void UpdateNeighborFrequencyPerCycle();
        
        double getTemperature();
        double getDensity();
        double getVolume();
        double getCosAngleMax();
        double getPressure();
        coord_list_t& getFullColloidListCoord();
        std::map<int,move_info>& getMoveInfoMap();
        
        void resetPressureTally();
        
        void run(int);
        void run();
        
        //file handling
        void openFile();
        void closeFile();
        void writeConfig();
        void openPressureFile();
        void closePressureFile();
        
        double computeEnergy(int);
        double computeEnergy();
        double computeMoleculeEnergy(int);
        double computeEnergy(MoleculeList& system, Box& box);
        

        
        void   computeVolume(); //make return
        double computePressure(double dv=0.01, PressureMode mode=COMPRESSION);
        void   tallyPressure();
        
        
    protected:
        int _n;
        double _beta,_temperature;
        double _density;
        double _volume;
        int _nsubmoves;
        int _update_move_frequency_per_cycle;
        double _E, _delE;
        bool _core_flag;
        double _cos_angle_max;
        double computePairEnergy(int);
        
        int attemptMove(int);
        int attemptSubMove(int);
        int _flag;
        
        coord_list_t _old_config;
        coord_t _new_config;
        TetramerPatchyColloid _old_molecule;
        Colloid _old_colloid;
        
        void Translation(int);
        void Rotation(int);
        void TranslationAndRotation(int);
        bool isRotationGood(int);
        
        void    saveConfig(int index);
        void    revertConfig(int index);
        void    updateMoveInfo(int);
        
        bool checkNeighborList(int);
        bool checkNeighborList(MoleculeList&, Box&);
        void computeMaxDisplacement();
        bool pressureCalculation;
        
        
        std::map<int,move_info> _move_info_map;
        
        std::pair<int,int> old_flag;
        
        struct PressureLog{
            PressureLog():
            count(0),pressure_sum(0),scale_factor(0.01),
            _pressure_config(COMPRESSION),inverseVolume(1.),
            frequency(100)
            {
                _inverse_scale_factor = -1.0*inverseVolume/(scale_factor);
            };
            
            unsigned long count;
            unsigned int frequency;
            double pressure_sum;
            double scale_factor, _inverse_scale_factor;
            PressureMode _pressure_config;
            double inverseVolume;
            
            
            void reset(){
                count = pressure_sum = 0;
            }
            
            void refresh(){
                assert(scale_factor > 0.);
                _inverse_scale_factor = inverseVolume/scale_factor;
                
                if (_pressure_config == COMPRESSION)
                    _inverse_scale_factor *= -1.0;
            }
            
        } _pressure_log;
        
        enum {TRANSLATE = 0, ROTATE = 1, ROTATEMOLECULE=2};
        enum {REJECT = 0, ACCEPT = 1};
        
        std::ofstream _ofile_pressure;
        
    };
}

#endif /* defined(__TetrahedralParticlesInConfinement__simulation_nvt_ensemble__) */
