//
//  Cell.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 4/1/16.
//  Copyright © 2016 Folarin Latinwo. All rights reserved.
//

#include "cell.hpp"

namespace TetrahedralParticlesInConfinement {
    Cell::Cell(){
        
    }
    
    Cell::~Cell(){
        
    }
    
    void Cell::clear(){
        _neighbor.clear();
        _tag.clear();
        _index_of.clear();
    }
    
    void Cell::clearContents(){
        _tag.clear();
        _index_of.clear();
    }
    
#pragma mark GETS
    /**
     @return the number of neighboring cell around the cell
     */
    int Cell::getNumberOfNeighbors() const{
        return (int)_neighbor.size();
    }
    
    /**
     @return the number of tags in the cell
     */
    int Cell::getNumberOfObjects(){
        return (int)_tag.size();
    }
    
    std::vector<int> Cell::getNeighbors(){
        return _neighbor;
    }
    
    std::vector<int> Cell::getContents(){
        return _tag;
    }
    
#pragma mark SETS
#pragma mark OTHER
    /**
     @param tag is the tag to be inserted
     */
    void Cell::insert(int tag){
        //save the index for fast removal
        //question would this be faster than find that works with stl
        _index_of[tag] = (int)_tag.size();
        _tag.push_back(tag);
    }
    
    /**
     @param tag is the tag to be removed
     */
    void Cell::remove(int tag){
        int index_of_removed_tag = _index_of[tag]; //an assert here would be useful
        int tmp_tag = _tag[_tag.size()-1];
        _tag[index_of_removed_tag] = tmp_tag;
        _tag.pop_back();
        
        _index_of[tmp_tag] = index_of_removed_tag;
        _index_of.erase(tag);
        
    }
    
#pragma mark NEIGHBOR OPERATIONS
    /**
     @param nbr is the index of a neighboring cell
     */
    void Cell::addNeighbor(int nbr){
        _neighbor.push_back(nbr);
    }
    
    
}
