//
// Created by LEI XU on 4/11/19.
//

#ifndef RASTERIZER_TRIANGLE_H
#define RASTERIZER_TRIANGLE_H

#include <eigen3/Eigen/Eigen>

using namespace Eigen;
using namespace std;
class Triangle
{
public:
    Vector3f v[3]; /*the original coordinates of the triangle, v0, v1, v2 in
                      counter clockwise order*/ //�����ε�ԭʼ���ꡣv0 ,v1, v2 ������ʱ�뷽������
                      /*Per vertex values*/ // ÿ��������Ϣ
    Vector3f color[3];      // color at each vertex; // ������ÿ���������ɫ
    Vector2f tex_coords[3]; // texture u,v    //����
    Vector3f normal[3];     // normal vector for each vertex  // ������ÿ������ķ�����

    // Texture *tex;
    Triangle(); // ���캯��

    Eigen::Vector3f a() const { return v[0]; }
    Eigen::Vector3f b() const { return v[1]; }
    Eigen::Vector3f c() const { return v[2]; }

    void setVertex(int ind, Vector3f ver); /*set i-th vertex coordinates ���õ�i����������*/
    void setNormal(int ind, Vector3f n);   /*set i-th vertex normal vector ���õ�i������ķ�����*/
    void setColor(int ind, float r, float g, float b); /*set i-th vertex color ���õ�i���������ɫ*/
    void setTexCoord(int ind, float s,
        float t); /*set i-th vertex texture coordinate���õ�i���������������*/
    array<Vector4f, 3> toVector4() const;
};

#endif // RASTERIZER_TRIANGLE_H
