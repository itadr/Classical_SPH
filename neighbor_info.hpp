#ifndef _NEIGHBOR_INFO_HPP
#define _NEIGHBOR_INFO_HPP
#include"vector2.hpp"

struct neighbor_info{
    int index;
    double dist_sq;
    vector2 v;
    neighbor_info(){}
    void set(const int i,const double r_sq,const vector2& vec){
        index=i;
        dist_sq=r_sq;
        v=vec;
    }
};

#endif
