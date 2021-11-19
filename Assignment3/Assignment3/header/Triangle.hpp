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
    /*the original coordinates of the triangle, v0, v1, v2 in counter clockwise order原始三角形的坐标，按照顺时针方向排列*/
    /*Per vertex values每个像素的值*/
    //color at each vertex;每个像素的颜色
    Vector3f color[3]; 
    //texture u,v 纹理的uv坐标
    Vector2f tex_coords[3]; 
    //normal vector for each vertex 每个顶点的法向量
    Vector3f normal[3]; 
    // 纹理
    Texture *tex= nullptr;  
    Triangle();

    Eigen::Vector4f a() const { return v[0]; }
    Eigen::Vector4f b() const { return v[1]; }
    Eigen::Vector4f c() const { return v[2]; }
    /*set i-th vertex coordinates 设置第i个顶点的坐标*/
    void setVertex(int ind, Vector4f ver); 
    /*set i-th vertex normal vector 设置第i个顶点的法向量*/
    void setNormal(int ind, Vector3f n); 
    /*set i-th vertex color 设置第i个顶点的rgb值*/
    void setColor(int ind, float r, float g, float b); 
    // 设置三个顶点的法向量
    void setNormals(const std::array<Vector3f, 3>& normals);   
    // 设置三个顶点的颜色
    void setColors(const std::array<Vector3f, 3>& colors);  
    /*set i-th vertex texture coordinate 设置第i个顶点的uv坐标*/
    void setTexCoord(int ind,Vector2f uv ); 
    std::array<Vector4f, 3> toVector4() const;
};






#endif //RASTERIZER_TRIANGLE_H
