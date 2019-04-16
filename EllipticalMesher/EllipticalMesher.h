//*****************************************************************************
//Title		:EllipticalMesher/EllipticalMesher.h
//Author	:Tanabe Yuta
//Date		:2019/04/16
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************


#pragma once
#include <vector>
#include <iostream>
#include "Point.h"


class EllipticalMesher
{
public:
	EllipticalMesher();
	~EllipticalMesher();
	EllipticalMesher(int, int);
	   
	
	void SetBoundaryPoint(int, Point);		//境界上の点を追加
	void MakeMesh();						//メッシュ生成
	void ExportMesh(std::string);			//メッシュ表示


private:
	const int inum;							//x軸方向の座標点数
	const int jnum;							//y軸方向の座標点数


	std::vector<Point> pb;					//境界上の点
	std::vector<std::vector<Point>> pin;	//内部の点


	static const int ITRMAX = 100000;		//連立方程式ソルバの最大反復回数
	const double DEPS = 1.0e-8;				//座標のノルムの許容残差
};

