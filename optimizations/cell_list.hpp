//
//  cell_list.hpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 4/4/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#ifndef cell_list_hpp
#define cell_list_hpp

#include <stdio.h>
#include "cell.hpp"
#include "struct_def.h"

namespace TetrahedralParticlesInConfinement {
    class CellList{
    public:
        CellList();                                 ///<Constructor
        virtual ~CellList();                        ///<Destructor
        
        int getNumbeOfCells() const;                ///<Returns number of Cells in the CellList
        int getCellContaining(coord_t&) const;      ///<Returns the cell containing a position in the box
        int getCellContaining(int);                 ///<Returns the cell containing an object
        int getNeighborStyle();                     ///<Returns the neighbor style
        
        
        coord_t getPositionOfCell(int);             ///<Returns cell at index i
        Cell getCell(int);                          ///
        
        void setBox(coord_t&);                      ///<sets the box dimensions
        void setBox(coord_t&, coord_t&);            ///<sets the box min,max
        void setBox(coord_t&, coord_t&,
                    bool_list_t&);                  ///<sets the box min, max, periodicity
        void setBox(Box&);                          ///<sets full box info
        
        void setMinimumCellSize(double);            ///<sets the smallest possible cell dimension
        void setInteractionRange(double);           ///<sets the range of the interaction potential
        void setNeighborStyle(int);                 ///<sets the neighbor style, or the way the cell neighbor list is built
        
        void clear();                               ///<Clears the contents of all the cells, their memory space, and the list of cells
        virtual void clearContents();               ///<Clears the contents of all the cells
        void placeInCorrectCell(int, coord_t&);     ///<If object has changed position, removes it from old cell, places it in new cell
        void insert(int, coord_t&);                 ///<Inserts and object into the cell list
        void remove(int);                           ///<Removes an object from the cell list
        
        bool good() const;                          ///<Indicates whether the cell list is valid
        bool fail() const;                          ///<Indicates whether the cell list is invalid
        
        bool end() const;                           ///<
        void resetIterator();
        std::pair<int,int> nextPair();
        
        /**
         @brief Special function to get the neighbors of an object
         */
        void getNeighborsOf(int, std::vector<int>&);
        
        /**
         @brief List of modes for setNeighborStyle
         
         For an MD-style cell list where interactions are pairwise, HALF is
         typically ised, since we only need half of the cell neighbors to
         compute the total force. For an MC-style cell list, FULL is typically
         used, since we often need to look at an individual particle and compute
         the force over particles in all neighboring cells
         */
        enum {FULL=0, HALF=1};
        
    protected:
        virtual void rebuild();                     ///<Builds the cell list
        virtual void divideBoxIntoCells();          ///<Divides the box into cells
        double _interaction_range;                  ///<The range of the interaction potential
        double _min_cell_size;                      ///<The smallest possible cell dimension (usually _interaction_range).
        int _neighbor_style;                        ///<An integer indicating how to build the cell list
        bool _force_rebuild;                        ///<If true, then rebuild the list no matter what next time Rebuild is called.
        bool _is_good;                              ///<Indicates whether the cell list is valid
        
        std::vector<Cell> _cell;                    ///<List of all cells in the cell list
        std::map<int, int> _cell_containing;        ///<The cell containing an object
        
        int _ncellyz;                               ///<number of cells in yz plane
        int _ncellx;                                ///<number of cells in x dim
        int _ncelly;                                ///<number of cells in y dim
        int _ncellz;                                ///<number of cells in z dim
        
        double _cell_size[3];                       ///<the size of a cell
        double _inv_cell_size[3];                   ///<the inverse of the size of a cell
        coord_t _box_min, _box_max;                 ///<The box boundary
        bool_list_t _periodic;                      ///<periodicity
        
    private:
        int _i, _j, _k, _n, _c;
        std::pair<int,int> _current_pair;
        bool _same_cell, _end;
        std::vector<int> _tag;
        std::map<int,int> _index_of_tag;            ///<The index of a given tag in the _tag array
        
        
        
        
        
    };
}

#endif /* cell_list_hpp */
