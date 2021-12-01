//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_INTERSECTION_H
#define RAYTRACING_INTERSECTION_H
#include "Vector.hpp"
#include "Material.hpp"
class Object;
class Sphere;
// 光线和物体表面接触点
struct Intersection
{
    Intersection(){
        happened=false;
        coords=Vector3f();
        normal=Vector3f();
        distance= std::numeric_limits<double>::max();
        obj =nullptr;
        m=nullptr;
    }
    // 是否接触
    bool happened;
    // 接触点坐标
    Vector3f coords;
    // 接触点法向量
    Vector3f normal;
    // 接触点到光源位置
    double distance;
    // 接触物体
    Object* obj;
    // 接触材质
    Material* m;
};
#endif //RAYTRACING_INTERSECTION_H
