//
// Created by goksu on 4/6/19.
//
// rasterizer : 光栅化
#pragma once //避免一个文件被引用多次

#include "Triangle.h"
#include <algorithm>
#include <eigen3/Eigen/Eigen>
using namespace Eigen;

namespace rst {
    enum class Buffers
    {
        Color = 1,
        Depth = 2
    };

    inline Buffers operator|(Buffers a, Buffers b)
    {
        return Buffers((int)a | (int)b);
    }

    inline Buffers operator&(Buffers a, Buffers b)
    {
        return Buffers((int)a & (int)b);
    }
    //绘制的图元类型，线或三角形
    enum class Primitive
    {
        Line,
        Triangle
    };

    /*
     * For the curious : The draw function takes two buffer id's as its arguments.
     * These two structs make sure that if you mix up with their orders, the
     * compiler won't compile it. Aka : Type safety
     * 绘制函数使用两个缓冲区的id作为参数，这两个结构体保证了如果混淆了两者的顺序，编译器会报错。
     * 即为：类型保护。
     * */
    struct pos_buf_id
    {
        int pos_id = 0;
    };

    struct ind_buf_id
    {
        int ind_id = 0;
    };

    class rasterizer
    {
    public:
        rasterizer(int w, int h);
        // 位置缓冲区
        pos_buf_id load_positions(const std::vector<Eigen::Vector3f>& positions);
        // 索引缓冲区(渲染顺序)
        ind_buf_id load_indices(const std::vector<Eigen::Vector3i>& indices);
        // 将内部的模型矩阵作为参数传递给光栅化器
        void set_model(const Eigen::Matrix4f& m);
        // 将视图变换矩阵设置为内部视图矩阵
        void set_view(const Eigen::Matrix4f& v);
        // 将内部的投影矩阵设为给定矩阵p,并传递给光栅化器
        void set_projection(const Eigen::Matrix4f& p);

        //将屏幕像素点（x,y）设为（r,g,b）的颜色， 并写入相应的帧缓冲区位置
        void set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color);

        // 清空缓冲区
        void clear(Buffers buff);
        //绘制缓冲区内容
        void draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, Primitive type);

        std::vector<Eigen::Vector3f>& frame_buffer() { return frame_buf; }

    private:
        // 划线
        void draw_line(Eigen::Vector3f begin, Eigen::Vector3f end);
        // 光栅化一个三角形的栅框，即三次划线
        void rasterize_wireframe(const Triangle& t);

    private:
        Eigen::Matrix4f model;  // 模型矩阵
        Eigen::Matrix4f view;   // 视图矩阵
        Eigen::Matrix4f projection; // 投影矩阵
        // 位置缓冲区
        std::map<int, std::vector<Eigen::Vector3f>> pos_buf;
        std::map<int, std::vector<Eigen::Vector3i>> ind_buf;

        // 绘制窗口中的每一个像素都需要颜色数据和深度数据
        // 帧缓冲对象，用于存储需要在屏幕上绘制的颜色数据
        std::vector<Eigen::Vector3f> frame_buf ;
        // 深度缓冲区
        std::vector<float> depth_buf ;
        int get_index(int x, int y);

        // 需要绘制的窗口的宽和高
        int width, height;

        int next_id = 0;
        int get_next_id() { return next_id++; }
    };
} // namespace rst
