//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_INTERSECTION_H
#define RAYTRACING_INTERSECTION_H
#include "Vector.hpp"
#include "Material.hpp"
class Object;
class Sphere;
// ���ߺ��������Ӵ���
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
    // �Ƿ�Ӵ�
    bool happened;
    // �Ӵ�������
    Vector3f coords;
    // �Ӵ��㷨����
    Vector3f normal;
    // �Ӵ��㵽��Դλ��
    double distance;
    // �Ӵ�����
    Object* obj;
    // �Ӵ�����
    Material* m;
};
#endif //RAYTRACING_INTERSECTION_H
