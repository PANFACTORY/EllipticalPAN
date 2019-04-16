//*****************************************************************************
//Title		:EllipticalMesher/Point.h
//Author	:Tanabe Yuta
//Date		:2019/04/16
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#pragma once
#include <array>
#include <cmath>

class Point
{
public:
	Point();
	~Point();
	Point(double, double);

	std::array<double, 2> x;

	Point operator+(const Point &);		//ベクトルの和
	Point operator-(const Point &);		//ベクトルの差
	Point operator*(const double &);	//ベクトルの実数倍
	Point operator/(const double &);	//ベクトルの実数商

	double Norm();						//ベクトルのノルムを計算
};

Point operator*(const double &, const Point &);

Point operator/(const double &, const Point &);