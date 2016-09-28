//
//  cell_list.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 4/4/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "cell_list.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

namespace TetrahedralParticlesInConfinement {
#pragma mark CONSTRUCT/DESTRUCT
    /**
     */
    CellList::CellList():
    _interaction_range(10.),
    _min_cell_size(-1),
    _neighbor_style(FULL),
    _is_good(false),
    _box_min(coord_t(3,0.)),
    _box_max(coord_t(3,0.)),
    _periodic(bool_list_t(3,true)){
        
    }
    
    /**
     */
    CellList::~CellList(){
        
    }
    
#pragma mark SET
    /**
     @param style is an int representing a given neighbor style
    */
    void CellList::setNeighborStyle(int style){
        assert(style == 0 || style == 1);
        _neighbor_style = style;
        _force_rebuild = true;
        rebuild();
    }
    
    /**
     @param box is the size of the boundary in x,y,z
    */
    void CellList::setBox(coord_t& box){
        for (unsigned int i=0; i<box.size(); i++) {
            _box_max[i] = box[i]/2.0;
            _box_min[i] = -box[i]/2.0;
        }
        rebuild();
    }
    
    /**
     @param box is the sie of the boundary in x,y,z
    */
    void CellList::setBox(coord_t& box_min, coord_t& box_max){
        for (unsigned int i=0; i<box_min.size(); i++) {
            _box_max[i] = box_max[i]/2.0;
            _box_min[i] = box_min[i]/2.0;
        }
        rebuild();
    }
    
    /**
     @param box is the sie of the boundary in x,y,z
    */
    void CellList::setBox(coord_t& box_min, coord_t& box_max, bool_list_t& periodic){
        for (unsigned int i=0; i<box_min.size(); i++) {
            _box_max[i] = box_max[i];
            _box_min[i] = box_min[i];
            _periodic[i] = periodic[i];
        }
        rebuild();
    }
    
    /**
     @param box is the sie of the boundary in x,y,z
     */
    void CellList::setBox(Box& box){
        _box_max = box.box_hi;
        _box_min = box.box_lo;
        _periodic = box.periodic;
        rebuild();
    }
    
    
    /**
     The interaction range determines which cells are neighbors.
     Normally, the cell size is slightly larger than the interaction
     range such that each cell has 26 neighbors. However, we can make
     the cells smaller and loop over more cell neighbors to more 
     closely approximate a spherical interaction range and avoid
     calculating interactions between non-interacting particles.
     
     @param range is the interaction range
     */
    void CellList::setInteractionRange(double range){
        _interaction_range = range;
        rebuild();
    }
    
    /**
     If the minimum cell size is not set, the default value is the interaction
     range.
     
     @param size is the minimum cell size
     */
    
    void CellList::setMinimumCellSize(double size){
        _min_cell_size = size;
        rebuild();
    }
    
#pragma mark GET
    /**
     @return number of cells in the list
     */
    int CellList::getNumbeOfCells() const{
        return (int)_cell.size();
    }
    
    /**
     */
    int CellList::getCellContaining(int tag){
        return _cell_containing[tag];
    }
    
    /**
     */
    int CellList::getCellContaining(coord_t& x) const {
        return ( _ncellyz * (int)((x[0] - _box_min[0])*_inv_cell_size[0]) +
                 _ncellz  * (int)((x[1] - _box_min[1])*_inv_cell_size[1]) +
                            (int)((x[2] - _box_min[2])*_inv_cell_size[2]));
    }
    
    /**
     Returns lower bounds for cell
     @param icell is cell index
     */
    coord_t CellList::getPositionOfCell(int icell){
        coord_t x_cell(3,0.);
        
        x_cell[0] = icell/(_ncelly*_ncellz)%_ncellx*_cell_size[0] + _box_min[0];
        x_cell[1] = icell/(_ncellz)%_ncelly*_cell_size[1] + _box_min[1];
        x_cell[2] = icell%_ncellz*_cell_size[2] + _box_min[2];
        
        return x_cell;
    }
    
    /**
     Returns a reference to a cell in the list. (FL. Is this really a reference?
     It seems like copy by value)
     @param i is the number of the referenced cell
     */
    Cell CellList::getCell(int i){
        return _cell[i];
    }
    
    /**
     Returns the neighbor style for the cell list
     */
    int CellList::getNeighborStyle(){
        return _neighbor_style;
    }
    
#pragma mark DO
    /**
     */
    void CellList::clear(){
        _cell.clear();
        _cell_containing.clear();
        _tag.clear();
    }
    
    /**
     */
    void CellList::clearContents(){
        for (auto& i : _cell) i.clearContents();
        _cell_containing.clear();
        _tag.clear();
    }
    
    /**
     @param object to be inserted
     */
    void CellList::insert(int tag, coord_t& x){
        int cell_number = getCellContaining(x);
        assert(cell_number >= 0 && cell_number < (int) _cell.size());
        _cell_containing[tag] = cell_number;
        _cell[cell_number].insert(tag);
        
        
        //keep track of a list of tags (used for next... functions)
        //book keeping
        _index_of_tag[tag] = (int)_tag.size();
        _tag.push_back(tag);
    }
    
    /**
     @param object is object to be removed
     */
    void CellList::remove(int tag){
        int c = _cell_containing[tag];
        _cell[c].remove(tag);
        _cell_containing.erase(tag);
        
        //keep track of a list of tags (used for next ... functions)
        //book keeping stuff
        int replace_index = _index_of_tag[tag];
        int last_tag = _tag[_tag.size()-1];
        _tag[replace_index] = last_tag;
        _index_of_tag[last_tag] = replace_index;  //likely use tag infor to cycle through
        _tag.pop_back();
    }
    
    /**
     @param object is an object in the cell list
     */
    void CellList::placeInCorrectCell(int tag, coord_t& x){
        int c = getCellContaining(x);
        if ( c != _cell_containing[tag]) { // no need to insert if object already exists in cell
            remove(tag);
            insert(tag, x);
        }
    }
    
#pragma mark CHECK
    /**
     */
    bool CellList::good() const {
        return _is_good;
    }
    
    /**
     */
    bool CellList::fail() const {
        return !_is_good;
    }
    
#pragma mark SPECIAL FUNCTIONS
    /**
     This function is called whenever the physical boundaries, cutoffs, etc.
     are changed.
     */
    void CellList::rebuild(){
        _is_good = true;
        coord_t period(3,0.); //maybe use box here
        for (unsigned int i=0; i<3; i++) period[i] = _box_max[i] - _box_min[i];
        
        double irx = period[0]/_interaction_range;
        double iry = period[1]/_interaction_range;
        double irz = period[2]/_interaction_range;
        
        if (!(irx >= 3.0 || !_periodic[0]) ||
            !(iry >= 3.0 || !_periodic[1]) ||
            !(irz >= 3.0 || !_periodic[2]) ) {
            _is_good = false;
        }
        
        if (_min_cell_size > 0.0) {
            _ncellx = (int)(period[0]/_min_cell_size);
            _ncelly = (int)(period[1]/_min_cell_size);
            _ncellz = (int)(period[2]/_min_cell_size);
        }
        else{
            _ncellx = (int)(period[0]/_interaction_range);
            _ncelly = (int)(period[1]/_interaction_range);
            _ncellz = (int)(period[2]/_interaction_range);
        }
        
        _ncellyz = _ncelly*_ncellz;
        
        _cell_size[0] = period[0]/_ncellx;
        _cell_size[1] = period[1]/_ncelly;
        _cell_size[2] = period[2]/_ncellz;
        
        for (unsigned int i=0; i<3; i++) _inv_cell_size[i] = 1.0/_cell_size[i];
        
        if (_ncellx*_ncelly*_ncellz != (int)_cell.size()) _force_rebuild = true;
        
        if (_force_rebuild) {
            clear();
            _cell.resize(_ncellx*_ncelly*_ncellz);
            divideBoxIntoCells();
            _force_rebuild = false;
        }
    }
    
    /**
     This routine physically divides the box into cells and determines the
     neighbor of each cell
     */
    void CellList::divideBoxIntoCells(){
        
        //Determine the number of neighbors surrounding a cell in each
        //dimension. Handle the case when (interaction range)/(cell dimension)
        //is an integer value using rounding (i.e., the "temp" variable
        //operations).
        //shell_thickness = number of layers of neighboring cells surrounding
        //each cell within the interaction range (int). This defines a sub box
        //within the region. Some cells on the corners of the sub box may not
        //be within the interaction range. We get rid of these later.
        
        int shell_thickness[3];
        // =  number of layers of neighboring cells surrounding each cell
        // within the interaction range
        for (int i=0; i<3; i++) shell_thickness[i] = (int)(ceil(_interaction_range*_inv_cell_size[i]));
        
        //loop over cells
        for (int i=0; i<_ncellx; i++) {
            for (int j=0; j<_ncelly; j++) {
                for (int k=0; k<_ncellz; k++) {
                    //identify current cell
                    //cells start filling up in z dimension, then y, then x
                    //
                    //- first term: k has no factor, just counts position in z
                    //              dimension
                    //- second term: j counts the number of times we have filled
                    //                  up the z dimension.
                    //- third term: i counts the number of times we have filled the
                    //              yz plane.
                    
                    int c = k + j*_ncellz + i*_ncellyz;
                    
                    //collect near neighbors of cell c
                    for (int x=-shell_thickness[0]; x<=shell_thickness[0]; x++) {
                        for (int y=-shell_thickness[1]; y<=shell_thickness[1]; y++) {
                            for (int z=-shell_thickness[2]; z<=shell_thickness[2]; z++) {
                                //don't include the cell itself as a neighbor
                                if (x==0 && y==0 && z==0) continue;
                                if (_neighbor_style == HALF) {
                                    if (x >0) continue;
                                    if (x == 0) {
                                        if (y > z) continue;
                                        if (z < 0 && y==z) continue;
                                    }
                                }
                                
                                
                                //handle non-interacting cells here
                                //(looped over region is a rectangular prism but interaction
                                //region is a sphere). If the center-to-center distance
                                //between cells is greater than interaction energy minus
                                //1 cell length in each dimension the cells do not interact
                                
                                //distance away from cell cell in x direction
                                double term1 = x*_cell_size[0];
                                if (x<0) term1 += _cell_size[0];
                                else if (x>0) term1 -= _cell_size[0];
                                
                                //distance away from cell cell in y direction
                                double term2 = y*_cell_size[1];
                                if (y<0) term2 += _cell_size[1];
                                else if (y>0) term2 -= _cell_size[1];
                                
                                //distance away from cell cell in z direction
                                double term3 = z*_cell_size[2];
                                if (z<0) term3 += _cell_size[2];
                                else if (z>0) term3 -= _cell_size[2];
                                
                                if (term1*term1 + term2*term2 + term3*term3 >
                                    _interaction_range*_interaction_range) continue;
                                
                                //check for boundary conditions
                                int tmp;
                                
                                //x:
                                tmp = i+x;
                                if (!_periodic[0]){
                                    if (tmp < 0 || tmp >= _ncellx) continue;
                                }
                                
                                //y:
                                tmp = j+y;
                                if (!_periodic[1]){
                                    if (tmp < 0 || tmp >= _ncelly) continue;
                                }
                                
                                //z:
                                tmp = k+z;
                                if (!_periodic[2]){
                                    if (tmp < 0 || tmp >= _ncellz) continue;
                                }
                                
                                //Integer division is for pbcs.
                                int neighbor =  (i+x+_ncellx)%_ncellx*_ncellyz +
                                                (j+y+_ncelly)%_ncelly*_ncellz +
                                                (k+z+_ncellz)%_ncellz;
                                _cell[c].addNeighbor(neighbor);
                                
                            }
                        }
                    }
                }
            }
        }
    }
    
    void CellList::getNeighborsOf(int i, std::vector<int>& nbrs){
        nbrs.clear();
        int c = _cell_containing[i];
        std::vector<int> nbrs_c = _cell[c].getContents();
        for (unsigned int j=0; j<nbrs_c.size(); j++) {
            if (nbrs_c[j] != i) nbrs.push_back(nbrs_c[j]);
        }
        
        std::vector<int> surr_cell = _cell[c].getNeighbors();
        for (unsigned int n=0; n<surr_cell.size(); n++) {
            std::vector<int> nbrs_n = _cell[surr_cell[n]].getContents();
            nbrs.insert(nbrs.end(), nbrs_n.begin(),nbrs_n.end());
        }
    }
    
    bool CellList::end() const{
        return _end;
    }
    
    std::pair<int, int> CellList::nextPair(){
        //if we are looking in the same cell for neighbors
        if (_same_cell) {
            //if j exceeds the number of objects in cell c
            if (!(_j < _cell[_c].getNumberOfObjects())) {
                //reset j to zero for next time it enters
                _j = 0;
                //reset neighbor cell number and neighbor number to zero
                _n = 0;
                _k = 0;
                _same_cell = false;
                nextPair();
                return _current_pair;
            }
            else{
                //otherwise skip objects w/ >= sequencer
                if (_cell[_c].getContents()[_j] >= _current_pair.first) {
                    _j++;
                    nextPair();
                    return _current_pair;
                }
                else{
                    _current_pair.second = _cell[_c].getContents()[_j];
                    _j++;
                    return _current_pair;
                }
            }
        }
        else{
            //if k exceeds the number of objects in cell c neighbor n
            if (!(_k < _cell[_cell[_c].getNeighbors()[_n]].getNumberOfObjects())) {
                _n++;           //move on to the next neighbor
                _k = 0;         //reset k to zero
                
                //if n exceeds the number of neighboring cells
                if (!(_n < _cell[_c].getNumberOfNeighbors())) {
                    _i++;       //move on to the next object
                    
                    //if the next object does not exist
                    if (!(_i < (int)_tag.size())) {
                        _end = true;
                        return _current_pair;
                    }
                    
                    //otherwise move on to the next cell
                    _current_pair.first = _tag[_i];
                    _same_cell = true;
                    _c = _cell_containing[_current_pair.first];
                }
                nextPair();
                return _current_pair;
            }
            else{
                //get the next object in the cell
                _current_pair.second = _cell[_cell[_c].getNeighbors()[_n]].getContents()[_k++];
                return _current_pair;
            }
        }
    }
    
    void CellList::resetIterator(){
        _end = false;
        _same_cell = true;
        _i = _j = _k = _n = 0;
        int t = _tag[_i];
        _c = _cell_containing[t];
        _current_pair.first = t;
    }
    
}