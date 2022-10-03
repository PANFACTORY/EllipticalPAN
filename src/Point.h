//*****************************************************************************
// Title		:src/cpp/Point.h
// Author	:Tanabe Yuta
// Date		:2019/10/25
// Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>

namespace EllipticalPAN {
template <class T>
class Point {
   public:
    Point() {
        this->x[0] = T();
        this->x[1] = T();
    }
    Point(T _x, T _y) {
        this->x[0] = _x;
        this->x[1] = _y;
    }
    Point(const Point<T>& point) {
        this->x[0] = point.x[0];
        this->x[1] = point.x[1];
    }
    ~Point(){};

    bool operator==(const Point<T>& point) const {
        return this->x[0] == point.x[0] && this->x[1] == point.x[1];
    }
    bool operator!=(const Point<T>& point) const { return !(*this == point); }

    Point<T>& operator=(const Point<T>& point) {
        if (this != &point) {
            this->x[0] = point.x[0];
            this->x[1] = point.x[1];
        }
        return *this;
    }
    Point<T>& operator+=(const Point<T>& point) {
        this->x[0] += point.x[0];
        this->x[1] += point.x[1];
        return *this;
    }
    Point<T>& operator-=(const Point<T>& point) {
        this->x[0] -= point.x[0];
        this->x[1] -= point.x[1];
        return *this;
    }
    Point<T>& operator*=(T a) {
        this->x[0] *= a;
        this->x[1] *= a;
        return *this;
    }
    Point<T>& operator/=(T a) {
        this->x[0] /= a;
        this->x[1] /= a;
        return *this;
    }

    const Point<T> operator+(const Point<T>& point) const {
        Point<T> ret = *this;
        return ret += point;
    }
    const Point<T> operator-(const Point<T>& point) const {
        Point<T> ret = *this;
        return ret -= point;
    }
    const Point<T> operator-() const { return Point<T>() - *this; }
    const Point<T> operator*(T a) const {
        Point<T> ret = *this;
        return ret *= a;
    }
    const Point<T> operator/(T a) const {
        Point<T> ret = *this;
        return ret *= a;
    }

    T& operator[](unsigned int index) { return this->x[index]; }

    T Dot(const Point<T>& point) const {
        return this->x[0] * point.x[0] + this->x[1] * point.x[1];
    }
    T Norm() const { return sqrt(this->Dot(*this)); }

   private:
    T x[2];
};

template <class U>
inline std::ostream& operator<<(std::ostream& out, const Point<U>& point) {
    out << point[0] << std::endl << point[1] << std::endl;
    return out;
}

template <class U>
inline Point<U> operator*(U a, const Point<U>& point) {
    return point * a;
}
}  // namespace EllipticalPAN
