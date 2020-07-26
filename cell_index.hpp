#ifndef _CELL_INDEX_HPP
#define _CELL_INDEX_HPP

// index of cell linked list for neighbor particle search
struct cell_index{
    int x,y;
    cell_index();
    cell_index(const int ix,const int iy);
    void set(const int ix,const int iy);
    bool is_same(const int ix,const int iy);
};

#endif
