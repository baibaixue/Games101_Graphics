//
// Created by LEI XU on 4/9/19.
//

#ifndef RASTERIZER_GLOBAL_H
#define RASTERIZER_GLOBAL_H
typedef unsigned char u08;
#define MY_PI 3.1415926
#define TWO_PI (2.0* MY_PI)
#define WIDTH 700
#define HIGHT 700
enum TextureType
{
	NormalShader = 0,
	BlinnPhong = 1,
	TextureShader = 2,
	BumpMapping = 3,
	DisplacementMapping = 4,
};
#endif //RASTERIZER_GLOBAL_H
