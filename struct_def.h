//
//  struct_def.h
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/28/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#ifndef TetrahedralParticlesInConfinement_struct_def_h
#define TetrahedralParticlesInConfinement_struct_def_h


#include <vector>
#include <string>
#include <complex>
#include <map>

namespace TetrahedralParticlesInConfinement {
    
    //perhaps quaternion and orientation should be a pair
    
    typedef std::vector<double> coord_t;
    typedef std::vector<bool> bool_list_t;
    typedef std::vector<coord_t> coord_list_t;
    typedef std::map<double,double> function1d_t;
    
    typedef std::map<int,coord_list_t> vector1d_t;
    
    typedef std::pair<double,coord_t> double_coord_t;
    
    typedef std::complex<double> component_t;
    typedef std::vector<component_t> shpdesc_t;
    
    typedef std::vector< std::complex<double> > component_list_t;

    
    
    typedef void* arg_t; //what is a void pointer?
    
    //note quaternion and orientation should be pair
    
    struct Colloid{
        
        Colloid():
        _center_of_mass(3,0.0),
        orientation(3),
        quaternion(4),
        diameter(1.),
        molecule_id(0),
        core(false),
        bound(false)
        {}
        
        coord_t _center_of_mass;
        coord_t orientation;
        coord_t quaternion;
        double diameter;
        int molecule_id;
        bool core;
        bool bound;
        
        void setCenterOfMass(const coord_t& com){
            _center_of_mass = com;
        }
        
        void setOrientationVector(const coord_t& _orientation){
            orientation = _orientation;
        }
        
        double getRadius(){ return 0.5*diameter;};
        
        double getDiameterSq(){return diameter*diameter;}
        
    };
    
    struct Box{
        
        Box():
        box_lo(3,0.0),
        box_hi(3,10.0),
        box_period(3,10.0),
        periodic(3,true)
        {}
        
        coord_t box_lo;
        coord_t box_hi;
        coord_t box_period;
        bool_list_t periodic;
        
        
        void clear(){
            box_lo.clear();
            box_hi.clear();
            box_period.clear();
        }
        
        void updatePeriod(){
            for (unsigned int i=0; i< box_lo.size(); i++)
                box_period[i] = box_hi[i] - box_lo[i];
        }
        
        
    };
    
    
    struct Lattice{
    public:
        Lattice():
        density(1.0),
        number_of_lattice_points(216)
        {}
        
        double box_length;
        
        virtual void generateLattice(){};
        coord_list_t getLattice() {return points;};
        std::string getLatticeType(){return lattice_type;};
        
        int getNumberOfLatticePoints(){
            return number_of_lattice_points;};
        
        void setNumberOfLatticePoints(int n){number_of_lattice_points = n;};
        void setDensity(double dens){density = dens;};
        
        
    protected:
        coord_list_t points;
        coord_list_t base_vectors;
        std::string lattice_type;
        double _lattice_spacing;
        double density;
        int number_of_lattice_points;
        
        
    };
    
    
    struct UmbrellaSpring{
        
        UmbrellaSpring():
        order_parameter(0.5),
        spring_constant(10.),
        umbrella_type("density"){
        }
        
        double  order_parameter; //in particular, this corresponds to average extension, assumes harmonic
        double  spring_constant;
        std::string umbrella_type;
        
        //methods
        double getUmbrellaEnergy(double val){
            return 0.5*spring_constant*pow(val - order_parameter, 2.0);
        }
        
    };
    
    struct Wall{
        
        Wall():
        type("WCA"),
        position(0.),
        location(BOTTOM),
        axis(Z){
            
        }
        
        enum Axis {X, Y, Z};
        enum Location {TOP, BOTTOM};

        std::string type;
        double position;
        Location location;
        Axis axis;
        
    };
    
    struct Plates{
        friend class Simulation;
    public:
        
        Plates():
        _separation(6.0){
            for (unsigned int i=0; i<2; i++) _walls.push_back(Wall());
            
            _walls[1].position = _separation;
            _walls[1].location = Wall::TOP;
        }
        
        void setPlateSeparation(double sep){ _separation = sep; _walls[1].position = _walls[0].position + _separation;}
        void setPlateAxis(Wall::Axis axis){_walls[0].axis = _walls[1].axis = axis;}
        std::vector<Wall>& getWalls(){return _walls;};
        
        
    protected:
        double _separation;
        std::vector<Wall> _walls;
        
    };
    
}

#endif
