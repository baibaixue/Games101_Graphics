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
                      counter clockwise order*/ //三角形的原始坐标。v0 ,v1, v2 按照逆时针方向排列
                      /*Per vertex values*/ // 每个顶点信息
    Vector3f color[3];      // color at each vertex; // 三角形每个定点的颜色
    Vector2f tex_coords[3]; // texture u,v    //纹理
    Vector3f normal[3];     // normal vector for each vertex  // 三角形每个顶点的法向量

    // Texture *tex;
    Triangle(); // 构造函数

    Eigen::Vector3f a() const { return v[0]; }
    Eigen::Vector3f b() const { return v[1]; }
    Eigen::Vector3f c() const { return v[2]; }

    void setVertex(int ind, Vector3f ver); /*set i-th vertex coordinates 设置第i个顶点坐标*/
    void setNormal(int ind, Vector3f n);   /*set i-th vertex normal vector 设置第i个顶点的法向量*/
    void setColor(int ind, float r, float g, float b); /*set i-th vertex color 设置第i个顶点的颜色*/
    void setTexCoord(int ind, float s,
        float t); /*set i-th vertex texture coordinate设置第i个顶点的纹理坐标*/
    array<Vector4f, 3> toVector4() const;
};

#endif // RASTERIZER_TRIANGLE_H
