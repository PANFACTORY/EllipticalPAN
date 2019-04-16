//*****************************************************************************
//Title		:EllipticalMesher/main.cpp
//Author	:Tanabe Yuta
//Date		:2019/04/16
//Copyright	:(C)2019 TanabeYuta
//*****************************************************************************

#include <iostream>
#include "EllipticalMesher.h"

int main() {
	EllipticalMesher Mesh(17, 5);		//メッシャーオブジェクト

	


	Mesh.MakeMesh();					//メッシュ生成
	Mesh.ExportMesh("mesh");			//VTKファイルに出力
	return 0;
}