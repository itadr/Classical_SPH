#include "vector2.hpp"

#include <cmath>
vector2::vector2() {}
vector2::vector2(const double vx, const double vy) {
  x = vx;
  y = vy;
}
void vector2::set(const double vx, const double vy) {
  x = vx;
  y = vy;
}
double vector2::square() const { return x * x + y * y; }
double vector2::abs() const { return std::sqrt(this->square()); }
double vector2::inner_product(const vector2& v) const {
  return x * v.x + y * v.y;
}
vector2 vector2::operator+(const vector2& v) const {
  return vector2(x + v.x, y + v.y);
}
vector2 vector2::operator-(const vector2& v) const {
  return vector2(x - v.x, y - v.y);
}
vector2 vector2::operator*(const double a) const {
  return vector2(a * x, a * y);
}
vector2 vector2::operator/(const double a) const {
  const double a_inv = 1 / a;
  return vector2(a_inv * x, a_inv * y);
}
vector2& vector2::operator+=(const vector2& v) {
  x += v.x;
  y += v.y;
  return *this;
}
vector2& vector2::operator-=(const vector2& v) {
  x -= v.x;
  y -= v.y;
  return *this;
}
vector2& vector2::operator*=(const double a) {
  x *= a;
  y *= a;
  return *this;
}
vector2& vector2::operator/=(const double a) {
  const double a_inv = 1 / a;
  x *= a_inv;
  y *= a_inv;
  return *this;
}
