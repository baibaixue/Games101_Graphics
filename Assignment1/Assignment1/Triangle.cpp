//
// Created by LEI XU on 4/11/19.
//

#include "Triangle.h"
#include <algorithm>
#include <array>
#include <stdexcept>
using namespace std;
// ���캯��
Triangle::Triangle()
{
    v[0] << 0, 0, 0;
    v[1] << 0, 0, 0;
    v[2] << 0, 0, 0;

    color[0] << 0.0, 0.0, 0.0;
    color[1] << 0.0, 0.0, 0.0;
    color[2] << 0.0, 0.0, 0.0;

    tex_coords[0] << 0.0, 0.0;
    tex_coords[1] << 0.0, 0.0;
    tex_coords[2] << 0.0, 0.0;
}

void Triangle::setVertex(int ind, Eigen::Vector3f ver) { v[ind] = ver; }

void Triangle::setNormal(int ind, Vector3f n) { normal[ind] = n; }

void Triangle::setColor(int ind, float r, float g, float b)
{
    if ((r < 0.0) || (r > 255.) || (g < 0.0) || (g > 255.) || (b < 0.0) ||
        (b > 255.)) {
        throw std::runtime_error("Invalid color values");
    }

    color[ind] = Vector3f((float)r / 255., (float)g / 255., (float)b / 255.);
    return;
}
void Triangle::setTexCoord(int ind, float s, float t)
{
    tex_coords[ind] = Vector2f(s, t);
}
//�������εĶ�������ת��Ϊ��ά��������ʽ������η��̣�
array<Vector4f, 3> Triangle::toVector4() const
{
    array<Vector4f, 3> res;
    transform(begin(v), end(v), res.begin(), [](auto& vec) {
        return Vector4f(vec.x(), vec.y(), vec.z(), 1.f);
        });
    return res;
}
// transform ����
/*
begin(v)��ԭ��������ʼ��ַ
end(v):ԭ��������ֹ��ַ
res.begin(): resΪ��Ž���ĵ�������begin() Ŀ����������ʼ��ַ
[](auto& vec){...} ����������ָ��
auto �� �Զ��ƶϱ������ͣ������ڱ�����
*/