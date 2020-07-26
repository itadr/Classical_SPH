#include "cell_index.hpp"

cell_index::cell_index() {}
cell_index::cell_index(const int ix, const int iy) {
  x = ix;
  y = iy;
}
void cell_index::set(const int ix, const int iy) {
  x = ix;
  y = iy;
}
bool cell_index::is_same(const int ix, const int iy) {
  return (x == ix) && (y == iy);
}
