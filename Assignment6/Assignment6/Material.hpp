//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_MATERIAL_H
#define RAYTRACING_MATERIAL_H

#include "Vector.hpp"
// 材质类型
enum MaterialType { DIFFUSE_AND_GLOSSY, REFLECTION_AND_REFRACTION, REFLECTION };
// 材质
class Material{
public:
    // 材质类型
    MaterialType m_type;
    // 颜色
    Vector3f m_color;
    // 发出光线
    Vector3f m_emission;
    // 折射率
    float ior;
    // 漫反射系数和镜面反射系数
    float Kd, Ks;
    // 镜面反射指数
    float specularExponent;
    //Texture tex;

    inline Material(MaterialType t=DIFFUSE_AND_GLOSSY, Vector3f c=Vector3f(1,1,1), Vector3f e=Vector3f(0,0,0));
    inline MaterialType getType();
    inline Vector3f getColor();
    inline Vector3f getColorAt(double u, double v);
    inline Vector3f getEmission();


};

Material::Material(MaterialType t, Vector3f c, Vector3f e){
    m_type = t;
    m_color = c;
    m_emission = e;
}

MaterialType Material::getType(){return m_type;}
Vector3f Material::getColor(){return m_color;}
Vector3f Material::getEmission() {return m_emission;}

Vector3f Material::getColorAt(double u, double v) {
    return Vector3f();
}
#endif //RAYTRACING_MATERIAL_H
