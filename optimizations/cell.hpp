//
//  Cell.hpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 4/1/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef Cell_hpp
#define Cell_hpp

#include <stdio.h>
#include "struct_def.h"

namespace TetrahedralParticlesInConfinement{
    
    /**
     @brief Contains a list of integers
     @author Aaron Keys (FL borrowed from ASK)
     @ingroup opt
     
     Cell represents one cell in class CellList, which contains and manages many Cells
     */
    class Cell{
    public:
        Cell();                                 //Constructor
        virtual ~Cell();                        //Destructor
        
        virtual void clear();                   //clear the cell memory space (contents, neighbors)
        virtual void clearContents();           //clear the contents of the cell but not neighbors
        
        virtual void insert(int);               //insert an object at the back of the cell
        virtual void remove(int);               //remove an object from the cell
        
        void addNeighbor(int);                  //add a neighboring cell to the list of neighbors
        
        int getNumberOfObjects();               //return the number of objects in the cell
        int getNumberOfNeighbors() const;       //return the number of neighboring cells
        std::vector<int> getNeighbors();        //return a list of neighbors
        std::vector<int> getContents();         //return a list of contents
        
    protected:
        std::vector<int> _neighbor;             //list containing the neighboring cells of a cell
        std::vector<int> _tag;                  //a container to hold the tags of the objectes of the cell
        std::map<int,int> _index_of;            //index of a given tag within the tag array
        
        
    };
}

#endif /* Cell_hpp */
