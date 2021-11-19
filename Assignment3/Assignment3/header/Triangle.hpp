//
// Created by LEI XU on 4/11/19.
//

#ifndef RASTERIZER_TRIANGLE_H
#define RASTERIZER_TRIANGLE_H

#include <eigen3/Eigen/Eigen>
#include "Texture.hpp"

using namespace Eigen;
class Triangle{

public:
    Vector4f v[3]; 
    /*the original coordinates of the triangle, v0, v1, v2 in counter clockwise orderԭʼ�����ε����꣬����˳ʱ�뷽������*/
    /*Per vertex valuesÿ�����ص�ֵ*/
    //color at each vertex;ÿ�����ص���ɫ
    Vector3f color[3]; 
    //texture u,v �����uv����
    Vector2f tex_coords[3]; 
    //normal vector for each vertex ÿ������ķ�����
    Vector3f normal[3]; 
    // ����
    Texture *tex= nullptr;  
    Triangle();

    Eigen::Vector4f a() const { return v[0]; }
    Eigen::Vector4f b() const { return v[1]; }
    Eigen::Vector4f c() const { return v[2]; }
    /*set i-th vertex coordinates ���õ�i�����������*/
    void setVertex(int ind, Vector4f ver); 
    /*set i-th vertex normal vector ���õ�i������ķ�����*/
    void setNormal(int ind, Vector3f n); 
    /*set i-th vertex color ���õ�i�������rgbֵ*/
    void setColor(int ind, float r, float g, float b); 
    // ������������ķ�����
    void setNormals(const std::array<Vector3f, 3>& normals);   
    // ���������������ɫ
    void setColors(const std::array<Vector3f, 3>& colors);  
    /*set i-th vertex texture coordinate ���õ�i�������uv����*/
    void setTexCoord(int ind,Vector2f uv ); 
    std::array<Vector4f, 3> toVector4() const;
};






#endif //RASTERIZER_TRIANGLE_H
