//
// Created by saman on 10/29/20.
//

#ifndef SIMPLEWIRE_VEC_H
#define SIMPLEWIRE_VEC_H

#include <iostream>
#include <vector>

#include <cmath>

using namespace std;

/* Simple 2D cartesian vector object for simple mass/spring calculations */
class Vec {
public:
    Vec() {
        this->x = 0.0;
        this->y = 0.0;
    };

    Vec(double x1, double y1) {
        this->x = x1;
        this->y = y1;
    };

    // Because implicitly-declaration is deprecated
    Vec(const Vec &obj) {
        this->x = obj.x;
        this->y = obj.y;
    };

    inline Vec& operator=(Vec const& v) {
        if (this != &v){
            x = v.x;
            y = v.y;
        }
        return *this;
    };

    inline Vec& operator=(double const &val) {
        this->x = val;
        this->y = val;

        return *this;
    };

    inline Vec operator+(const Vec& v) const {
        return Vec(x + v.x, y + v.y);
    };

    inline Vec& operator*(const double& val) {
        this->x = val*x;
        this->y = val*y;
        return *this;
    };

    inline Vec setToScale(const double& scale) const {
        return Vec(x*scale, y*scale);

    };

    inline Vec& operator+=(Vec &v) {
        this->x = x + v.x;
        this->y = y + v.y;
        return *this;
    };

    inline Vec& operator*=(double const& scale) {
        this->x = x*scale;
        this->y = y*scale;
        return *this;
    };

    inline Vec operator-(const Vec& v) const {
        return Vec(x - v.x, y - v.y);
    };

    double Magnitude() const {
        return sqrt(x*x + y*y);
    };

    Vec UnitVec() const {
      return Vec(x/Magnitude(), y/Magnitude());
    };

    void PrintVec() const{
        std::cout <<"("<< x << ", " << y << ")" << std::endl;
    };

    double Dot(const Vec& v) const {
        return x*v.x + y*v.y;
    };

    ~Vec() = default;

    double y;
    double x;
};




#endif //SIMPLEWIRE_VEC_H
