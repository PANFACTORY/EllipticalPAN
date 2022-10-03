//*****************************************************************************
// Title		:src/cpp/Point.h
// Author	:Tanabe Yuta
// Date		:2019/10/25
// Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#pragma once
#include <array>
#include <cmath>

namespace EllipticalPAN {
template <class T>
class Point {
   public:
    Point() {
        this->x[0] = T();
        this->x[1] = T();
    };
    ~Point(){};
    Point(T _x, T _y) {
        this->x[0] = _x;
        this->x[1] = _y;
    }

    const Point<T> operator+(const Point<T>& point) const {
        return Point<T>(this->x[0] + point.x[0], this->x[1] + point.x[1]);
    }
    const Point<T> operator-(const Point<T>& point) const {
        return Point<T>(this->x[0] - point.x[0], this->x[1] - point.x[1]);
    }
    const Point<T> operator*(T _a) const {
        return Point<T>(this->x[0] * _a, this->x[1] * _a);
    }
    const Point<T> operator/(T _a) const {
        return Point<T>(this->x[0] / _a, this->x[1] / _a);
    }

    T& operator[](unsigned int index) { return this->x[index]; }

    T Dot(const Point<T>& _point) const {
        return this->x[0] * _point.x[0] + this->x[1] * _point.x[1];
    }
    T Norm() const { return sqrt(this->Dot(*this)); }

   private:
    T x[2];
};

template <class T>
inline Point<T> operator*(T _a, const Point<T>& _point) {
    return _point * _a;
}
}  // namespace EllipticalPAN
