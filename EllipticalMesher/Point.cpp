//*****************************************************************************
//Title		:EllipticalMesher/Point.cpp
//Author	:Tanabe Yuta
//Date		:2019/04/16
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#include "Point.h"


Point::Point(){}


Point::~Point(){}


Point::Point(double _x, double _y){
	this->x[0] = _x;
	this->x[1] = _y;
}


Point Point::operator+(const Point &_point){
	return Point(this->x[0] + _point.x[0], this->x[1] + _point.x[1]);
}


Point Point::operator-(const Point &_point) {
	return Point(this->x[0] - _point.x[0], this->x[1] - _point.x[1]);
}


Point Point::operator*(const double &a){
	return Point(this->x[0] * a, this->x[1] * a);
}


Point operator*(const double &a, const Point& _point) {
	return Point(_point.x[0] * a, _point.x[1] * a);
}


Point Point::operator/(const double &b){
	return Point(this->x[0] / b, this->x[1] / b);
}


Point operator/(const double &b, const Point& _point) {
	return Point(_point.x[0] / b, _point.x[1] / b);
}


double Point::Norm(){
	return sqrt(pow(this->x[0], 2.0) + pow(this->x[1], 2.0));
}