//*****************************************************************************
//Title		:src/cpp/Point.h
//Author	:Tanabe Yuta
//Date		:2019/10/25
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <array>
#include <cmath>


namespace EllipticalPAN{
	template<class T>
	class Point
	{
public:
		Point();
		~Point();
		Point(T _x, T _y);


		std::array<T, 2> x;


		Point<T> operator+(const Point<T>& _point);		
		Point<T> operator-(const Point<T>& _point);		
		Point<T> operator*(T _a);	
		Point<T> operator/(T _a);	


		T Norm();						
	};


	template<class T>
	Point<T>::Point(){}


	template<class T>
	Point<T>::~Point(){}


	template<class T>
	Point<T>::Point(T _x, T _y){
		this->x[0] = _x;
		this->x[1] = _y;
	}


	template<class T>
	Point<T> Point<T>::operator+(const Point<T> &_point){
		return Point<T>(this->x[0] + _point.x[0], this->x[1] + _point.x[1]);
	}


	template<class T>
	Point<T> Point<T>::operator-(const Point<T> &_point) {
		return Point<T>(this->x[0] - _point.x[0], this->x[1] - _point.x[1]);
	}


	template<class T>
	Point<T> Point<T>::operator*(T _a){
		return Point<T>(this->x[0] * _a, this->x[1] * _a);
	}


	template<class T>
	Point<T> Point<T>::operator/(T _a){
		return Point<T>(this->x[0] / _a, this->x[1] / _a);
	}


	template<class T>
	T Point<T>::Norm(){
		return sqrt(pow(this->x[0], 2.0) + pow(this->x[1], 2.0));
	}


	template<class T>
	Point<T> operator*(T _a, const Point<T>& _point) {
		return Point<T>(_point.x[0] * _a, _point.x[1] * _a);
	}
}