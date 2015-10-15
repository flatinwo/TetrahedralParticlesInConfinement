//
//  io.cpp
//  TetrahedralParticlesInConfinement
//
//  Created by Folarin Latinwo on 7/29/15.
//  Copyright (c) 2015 Folarin Latinwo. All rights reserved.
//

#include "io.h"
#include "spatial.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <map>
#include <assert.h>
#include <stack>

namespace TetrahedralParticlesInConfinement {
    
    //Brief: Load raw coordinates
    void load(const char* filename, coord_list_t& x, arg_t arg){
        
        std::ifstream file(filename);
        if (file.fail()) {
            std::cerr << "Error: load: can't openfile "
            << filename << ".\n";
            exit(1);
        }
        
        std::string str;
        while (std::getline(file,str)){
            std::istringstream iss(str.c_str());
            double temp;
            coord_t xi;
            while (iss >> temp) {
                xi.push_back(temp);
            }
            if (xi.size() != 0) {
                x.push_back(xi);
            }
        }
    }
    
    
    // delete type from coordinate list
    void deltype(coord_list_t& x, std::vector<std::string>& types, std::string t){
        assert(x.size() == types.size());
        
        coord_list_t xtmp = x;
        std::vector<std::string> typetmp = types;
        
        types.clear();
        x.clear();
        
        for (unsigned int i=0; i<xtmp.size(); i++) {
            if (typetmp[i] != t) {
                x.push_back(xtmp[i]);
                types.push_back(typetmp[i]);
            }
        }
    }
    
    
    //converts file to ifstream and loads!
    
    void loadxyz(const char* filename, xyzfile& xyz){
        std::ifstream file(filename);
        assert(file.good());
        loadxyz(file,xyz);
        file.close();
    }
    
    
    void loadxyz(std::istream& is, xyzfile& xyz){
        int N=0;
        is >> N;
        xyz.n = N;
        std::getline(is,xyz.commentstr); // will get rest of first line
        std::getline(is,xyz.commentstr);
        
        xyz.x.resize(N,std::vector<double>(3,0.0));
        xyz.x.resize(N);
        
        for (int i=0; i<N; i++) {
            double x,y,z;
            std::string typestr;
            is >> typestr >> x >> y >> z;
            xyz.type[i] = typestr;
            xyz.x[i][0] = x;
            xyz.x[i][1] = y;
            xyz.x[i][2] = z;
        }
        
    }
    
    
    void loadxyz(const char* filename, coord_list_t& x, arg_t arg){
        xyz_info* info = (xyz_info*) arg;
        xyz_info temp_info;
        
        if (arg == 0x0){
            info = &temp_info;
        }
        
        std::istream* is;
        std::ifstream ifs;
        
        if (info->instream != 0x0){
            is = info->instream;
        }
        else{
            ifs.open(filename);
            is = &ifs;
        }
        
        assert(is->good());
        
        int n = 0;
        std::string str;
        std::getline(*is, str);
        std::istringstream iss(str.c_str());
        iss >> n;
        std::getline(*is,str);
        x.clear();
        info->type.clear();
        
        for (int i=0; i<n; i++) {
            std::string typei;
            double xi, yi, zi;
            *is >> typei >> xi >> yi >> zi;
            info->type.push_back(typei);
            coord_t xyz(3);
            xyz[0] = xi; xyz[1] = yi; xyz[2] = zi;
            x.push_back(xyz);
        }
        std::getline(*is,str);
        
        if (info->instream == 0x0) {
            ifs.close();
        }
        
    }
    
    
    void loadxyz(const char* filename, MoleculeList& system, Box& box, arg_t arg){
        xyz_info* info = (xyz_info*) arg;
        xyz_info temp_info;
        
        if (arg == 0x0){
            info = &temp_info;
        }
        
        std::istream* is;
        std::ifstream ifs;
        
        if (info->instream != 0x0){
            is = info->instream;
        }
        else{
            ifs.open(filename);
            is = &ifs;
        }
        
        assert(is->good());
        
        //get total number of particles
        int n = 0;
        std::string str, word;
        std::getline(*is, str);
        std::istringstream iss(str.c_str());
        iss >> n;
        assert(n>0);
        
        //now attempt to get box size
        //a quick long winded approach and assumes that box_lo = zeroes
        std::getline(*is,str);
        std::stack<std::string> text;
        std::istringstream iss1(str);
        while (iss1.good()) {
            iss1 >> word;
            text.push(word);
        }
        
        //it is iterator
        for (auto it=box.box_hi.rbegin(); it != box.box_hi.rend(); ++it) {
            *it = std::stod(text.top());
            text.pop();
        }
        box.updatePeriod();
        
        
        coord_list_t x,r;
        info->type.clear();
        
        for (int i=0; i<n; i++) {
            std::string typei;
            double xi, yi, zi, rxi, ryi, rzi, wxi = 0.0;
            
            if (i%5 != 0) *is >> typei >> xi >> yi >> zi >> rxi >> ryi >> rzi;
            else *is >> typei >> xi >> yi >> zi >> wxi >> rxi >> ryi >> rzi;
            
            info->type.push_back(typei);
            coord_t xyz(3), rxyz(3);
            xyz[0] = xi; xyz[1] = yi; xyz[2] = zi;
            
            if (i%5 != 0){
               rxyz[0] = rxi; rxyz[1] = ryi; rxyz[2] = rzi;
            }
            else{
                rxyz.resize(4);
                rxyz[0] = wxi; rxyz[1] = rxi; rxyz[2] = ryi; rxyz[3] = rzi;
            }

            x.push_back(xyz);
            
            r.push_back(rxyz);
            if (x.size() == 5) {
                system.addToMoleculeList(x, r, box);
                x.clear();
                r.clear();
            }
        }
        std::getline(*is,str);
        
        if (info->instream == 0x0) {
            ifs.close();
        }
        
    }
    
    
    //saves data with variate types
    
    void savevarxyz(const char* filename, coord_list_t& x, xyz_info& info){
        assert(info.type.size()>0);
        assert(info.reservoir.size()>0);
        std::ofstream xyzfile(filename);
        
        std::ostream *os;
        
        if (info.outstream == NULL) {
            assert(xyzfile.good());
            os = &xyzfile;
        }
        else{
            assert(info.outstream->good());
            os = info.outstream;
        }
        
        int n=0;
        typedef std::map<std::string, int> reservoir_t;
        reservoir_t::iterator r;
        
        for (r = info.reservoir.begin(); r!=info.reservoir.end(); ++r){
            n += r->second;
        }
        
        *os << n << "\n\n";
        
        for (r=info.reservoir.begin(); r!=info.reservoir.end();++r){
            std::string type = r->first;
            int npad = r->second;
            for (unsigned int i=0; i<x.size(); i++) {
                if (info.type[i] == type) {
                    *os << type << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";
                    npad--;
                }
            }
            for (int i=0; i<npad; i++) {
                *os << type << "\t" << info.xres[0] <<"\t"<< info.xres[1] <<"\t"<< info.xres[2] <<"\n";
            }
        }
        xyzfile.close();
    }
    
    
    //have a different version with xyz type
    void savexyz(const char* filename, coord_list_t& x, xyz_info& info){
        unsigned long n = x.size();
        
        std::ofstream xyzfile(filename);
        
        std::ostream *os;
        
        if (info.outstream == NULL){
            assert(xyzfile.good());
            os = &xyzfile;
        }
        else{
            assert(info.outstream->good());
            os = info.outstream;
        }
        
        *os << n << "\n\n";
        if (info.type.size() > 0) {
            for (int i=0; i<n; i++){
                *os << info.type[i] << "\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";
            }
        }
        else{
            for (int i=0; i<n; i++){
                *os << "H\t" << x[i][0] << "\t" << x[i][1] << "\t" << x[i][2] << "\n";
            }
        }
        xyzfile.close();
    }
    
    std::ostream& operator << (std::ostream& os, Box& box){
        os << box.box_lo << "\n" << box.box_hi << "\n" << box.box_period << "\n";
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, coord_list_t& x){
        for (unsigned int i=0; i<x.size(); i++){
            os << x[i] << "\n";
        }
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, coord_t& x){
        unsigned long sm1 = x.size()-1;
        for (unsigned long k=0; k<x.size();k++){
            os << x[k];
            if (k!=sm1) {
                os << "\t";
            }
        }
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, shpdesc_t& sd){
        unsigned long sm1 = sd.size() - 1;
        for (unsigned long k=0; k<sd.size(); k++) {
            os << sd[k];
            if (k!=sm1) {
                os << "\t";
            }
        }
        return os;
    }
    
    std::ostream& operator << (std::ostream& os, function1d_t f){
        for (function1d_t::iterator i=f.begin(); i!=f.end(); ++i){
            os << i->first << "\t" << i->second << "\n";
        }
        return os;
    }
    
    std::istream& operator >> (std::istream& is, coord_list_t& x){
        return is;
    }
    
    std::istream& operator >> (std::istream& is, coord_t& x){
        return is;
    }
    
    std::istream& operator >> (std::istream& is, shpdesc_t& sd){
        return is;
    }
}