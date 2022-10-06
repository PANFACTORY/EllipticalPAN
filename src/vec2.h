/**
 * @file vec2.h
 * @author PANFACTORY (github/PANFACTORY)
 * @brief Define and implement the 2D vector class
 * @version 0.1
 * @date 2022-10-06
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include <cmath>
#include <iostream>

namespace EllipticalPAN {
/**
 * @brief Vec2 class
 *
 * @tparam T Type of element for example double
 */
template <class T>
class Vec2 {
    template <class U>
    friend std::ostream& operator<<(std::ostream& out, const Vec2<U>& vec2);

   public:
    /**
     * @brief Default constructor whose elements are (0, 0)
     *
     */
    Vec2() {
        this->x[0] = T();
        this->x[1] = T();
    }

    /**
     * @brief Construct a new Vec2 object whose elements are (_x, _y)
     *
     * @param _x Value of x coordinate
     * @param _y Value of y coordinate
     */
    Vec2(T _x, T _y) {
        this->x[0] = _x;
        this->x[1] = _y;
    }

    /**
     * @brief Copy constructor
     *
     * @param vec2 Copy source
     */
    Vec2(const Vec2<T>& vec2) {
        this->x[0] = vec2.x[0];
        this->x[1] = vec2.x[1];
    }

    ~Vec2(){};

    /**
     * @brief Check the same
     *
     * @param vec2     Comparison
     * @return true     Same with vec2
     * @return false    Not the same
     */
    bool operator==(const Vec2<T>& vec2) const {
        return this->x[0] == vec2.x[0] && this->x[1] == vec2.x[1];
    }

    /**
     * @brief Check not the same
     *
     * @param vec2     Comparison
     * @return true     Not the same
     * @return false    Same with vec2
     */
    bool operator!=(const Vec2<T>& vec2) const { return !(*this == vec2); }

    /**
     * @brief Assign vec2 to this object
     *
     * @param vec2         Assignor object
     * @return Vec2<T>&    Reference of this object
     */
    Vec2<T>& operator=(const Vec2<T>& vec2) {
        if (this != &vec2) {
            this->x[0] = vec2.x[0];
            this->x[1] = vec2.x[1];
        }
        return *this;
    }

    /**
     * @brief Compound assignment operator for addition
     *
     * @param vec2         Added object
     * @return Vec2<T>&    Reference of this object
     */
    Vec2<T>& operator+=(const Vec2<T>& vec2) {
        this->x[0] += vec2.x[0];
        this->x[1] += vec2.x[1];
        return *this;
    }

    /**
     * @brief Compound assignment operator for subtraction
     *
     * @param vec2         Subtracted object
     * @return Vec2<T>&    Reference of this object
     */
    Vec2<T>& operator-=(const Vec2<T>& vec2) {
        this->x[0] -= vec2.x[0];
        this->x[1] -= vec2.x[1];
        return *this;
    }

    /**
     * @brief Compound assignment operator for multiplication
     *
     * @param a             Coefficient of multiplication
     * @return Vec2<T>&    Reference of this object
     */
    Vec2<T>& operator*=(T a) {
        this->x[0] *= a;
        this->x[1] *= a;
        return *this;
    }

    /**
     * @brief Compound assignment operator for division
     *
     * @param a             Coefficient of division
     * @return Vec2<T>&    Reference of this object
     */
    Vec2<T>& operator/=(T a) {
        this->x[0] /= a;
        this->x[1] /= a;
        return *this;
    }

    /**
     * @brief Addition operator
     *
     * @param vec2             Adding vector
     * @return const Vec2<T>   Added vector
     */
    const Vec2<T> operator+(const Vec2<T>& vec2) const {
        Vec2<T> ret = *this;
        return ret += vec2;
    }

    /**
     * @brief Subtraction operator
     *
     * @param vec2             Subtracting vector
     * @return const Vec2<T>   Subtracted vector
     */
    const Vec2<T> operator-(const Vec2<T>& vec2) const {
        Vec2<T> ret = *this;
        return ret -= vec2;
    }

    /**
     * @brief Get reverse vector
     *
     * @return const Vec2<T>   Reverse vector
     */
    const Vec2<T> operator-() const { return Vec2<T>() - *this; }

    /**
     * @brief Get multipled vector
     *
     * @param a                 Multiple coefficient
     * @return const Vec2<T>   Multipled vector
     */
    const Vec2<T> operator*(T a) const {
        Vec2<T> ret = *this;
        return ret *= a;
    }

    /**
     * @brief Get dividev vector
     *
     * @param a                 Divide coefficient
     * @return const Vec2<T>   Devided vector
     */
    const Vec2<T> operator/(T a) const {
        Vec2<T> ret = *this;
        return ret *= a;
    }

    /**
     * @brief Get i th element value
     *
     * @param index Index of element
     * @return T&   Reference of i th element
     */
    T& operator[](unsigned int index) { return this->x[index]; }

    /**
     * @brief Get i th element value
     *
     * @param index Index of element
     * @return T&   Reference of i th element
     */
    T& operator()(unsigned int index) { return this->x[index]; }

    /**
     * @brief Get scalar product
     *
     * @param vec2 Vector used scalar product
     * @return T    Scalar product
     */
    T dot(const Vec2<T>& vec2) const {
        return this->x[0] * vec2.x[0] + this->x[1] * vec2.x[1];
    }

    /**
     * @brief Get norm of this vector
     *
     * @return T    Norm of this vector
     */
    T norm() const { return sqrt(this->dot(*this)); }

   private:
    T x[2];
};

/**
 * @brief Add vector to stream
 *
 * @tparam U                Type of coefficient and vector
 * @param out               Stream
 * @param vec2              Source vector
 * @return std::ostream&    Stream updated
 */
template <class U>
inline std::ostream& operator<<(std::ostream& out, const Vec2<U>& vec2) {
    out << vec2.x[0] << " " << vec2.x[1] << std::endl;
    return out;
}

/**
 * @brief Get multipled vector
 *
 * @tparam U        Type of coefficient and vector
 * @param a         Multiple coefficient
 * @param vec2      Applied vector
 * @return Vec2<U>  Multipled vector
 */
template <class U>
inline Vec2<U> operator*(U a, const Vec2<U>& vec2) {
    return vec2 * a;
}
}  // namespace EllipticalPAN
