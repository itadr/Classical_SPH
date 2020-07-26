#ifndef _VECTOR2_HPP
#define _VECTOR2_HPP

struct vector2{
    double x,y;
    vector2();
    vector2(const double vx,const double vy);
    void set(const double vx,const double vy);
    double square()const;
    double abs()const;
    double inner_product(const vector2& v)const;
    vector2 operator+(const vector2& v)const;
    vector2 operator-(const vector2& v)const;
    vector2 operator*(const double a)const;
    vector2 operator/(const double a)const;
    vector2& operator+=(const vector2& v);
    vector2& operator-=(const vector2& v);
    vector2& operator*=(const double a);
    vector2& operator/=(const double a);
};

#endif
