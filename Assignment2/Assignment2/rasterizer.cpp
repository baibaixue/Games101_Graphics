// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.h"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f>& positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return { id };
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i>& indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return { id };
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f>& cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return { id };
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

//判断点（x,y）是否在三角形中
static bool insideTriangle(float x, float y, const Vector3f* _v)
{
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    Eigen::Vector3f point = Eigen::Vector3f(x, y, 1.f);

    float res[3];

    for (int i = 0; i < 3; i++)
    {
        Eigen::Vector3f VP = point - _v[i];
        Eigen::Vector3f VV = _v[(i + 1) % 3] - _v[i];
        res[i] = VP.x() * VV.y()* 1.0 - VP.y() * VV.x() * 1.0;
    }
    for (int i = 0; i < 3; i++)
    {
        int j = (i + 1) % 3;
        if (res[i] * res[j] <= 0) return false;
    }
    return true;
}
/* 计算重心，c1,c2,c3分别表示点（x,y）在三角形v0，v1,v2三点的差值权重，即：
 * x = c1 * v[0].x() + c2 * v[1].x() + c3 * v[2].x();
 * y = c1 * v[0].y() + c2 * v[1].y() + c3 * v[2].y();
 * 1 = c1 + c2 + c3;
 * eg: tuple元组应用 : c++11中的新特性，让函数返回多个值
 * tuple是一个可以容纳不同元素的容器
 */
static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * y + v[1].x() * v[2].y() - v[2].x() * v[1].y()) / (v[0].x() * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * v[0].y() + v[1].x() * v[2].y() - v[2].x() * v[1].y());
    float c2 = (x * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * y + v[2].x() * v[0].y() - v[0].x() * v[2].y()) / (v[1].x() * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * v[1].y() + v[2].x() * v[0].y() - v[0].x() * v[2].y());
    float c3 = (x * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * y + v[0].x() * v[1].y() - v[1].x() * v[0].y()) / (v[2].x() * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * v[2].y() + v[0].x() * v[1].y() - v[1].x() * v[0].y());
    return { c1,c2,c3 };
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        //std::cout <<"buf-(z,w): "<< buf[i[0]].z()<<",1" << std::endl;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //std::cout << "afterMVP-(z,w): " << v[0].z() << "," << v[0].w() << std::endl;
        //Homogeneous division
        for (auto& vec : v) {
            vec.x() /= vec.w();
            vec.y() /= vec.w();
            vec.z() /= vec.w();
        }
        //std::cout << "after/w-(z,w): " << v[0].z() << "," << v[0].w() << std::endl;
        //Viewport transformation
        for (auto& vert : v)
        {
            vert.x() = 0.5 * width * (vert.x() + 1.0);
            vert.y() = 0.5 * height * (vert.y() + 1.0);
            vert.z() = vert.z() * f1 + f2;
        }
        //std::cout << "after-viewport-(z,w): " << v[0].z() << "," << v[0].w() << std::endl;
        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i]);
            //t.setVertex(i, v[i].head<3>());
        }
        //std::cout << v[0].z() << std::endl;
        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.v;
    Eigen::Vector3f v3f[3];
    for (int i = 0; i < 3; i++)
    {
        v3f[i] = t.v[i].head<3>();
    }
    // Bounding Box
    float minx = width, miny = height, maxx = 0, maxy = 0;
    for (auto i = 0; i < 3; i++)
    {

        minx = std::min(minx, v[i].x());
        miny = std::min(miny, v[i].y());
        maxx = std::max(maxx, v[i].x());
        maxy = std::max(maxy, v[i].y());
    }
    minx = std::max(0.f, minx);
    miny = std::max(0.f, miny);
    maxx = std::max(0.f, maxx);
    maxy = std::max(0.f, maxy);
    //std::cout << "rasterize_triangle" << minx << "," << miny << "," << maxx<< "," << maxy << std::endl;
    // TODO : Find out the bounding box of current triangle.找到当前三角形的Bounding Box
    // iterate through the pixel and find if the current pixel is inside the triangle 遍历BoundingBox中的每个像素点，判断该点是否在三角形内部

    // If so, use the following code to get the interpolated z value.   // 计算插入值的z深度值
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.设置三角形颜色

    for (int w = minx; w < std::min((float)width, maxx); w++)
    {
        for (int h = miny; h < std::min((float)height, maxy); h++)
        {
            //std::cout <<"minx:"<<minx<<" maxx:"<<maxx<<" miny:"<<miny<<" maxy:"<<maxy<<" (x,y):"<< w << "," << h << std::endl;
            Eigen::Vector3f point = Eigen::Vector3f(w, h, 1.f);
            if (insideTriangle(w , h , v3f))
            {
                float _alpha, _beta, _gamma;
                std::tie(_alpha, _beta, _gamma) = computeBarycentric2D(w , h , v3f);
                float w_reciprocal = 1.0 / (_alpha / v[0].w() + _beta / v[1].w() + _gamma / v[2].w());
                float z_interpolated = _alpha * v[0].z() / v[0].w() + _beta * v[1].z() / v[1].w() + _gamma * v[2].z() / v[2].w();
                //float z_interpolated = _alpha * v[0].z() + _beta * v[1].z()  + _gamma * v[2].z();
                z_interpolated *= w_reciprocal;
                //std::cout << w_reciprocal << ", "<< z_interpolated << std::endl;
                MSAA(z_interpolated, point, t, v3f);
            }
        }
    }
}

void rst::rasterizer::MSAA(float z_interpolated,const Eigen::Vector3f& point, const Triangle& t, const Eigen::Vector3f* v3f)
{
    //Eigen::Vector3f v3f[3];
    //for (int i = 0; i < 3; i++)
    //{
    //    v3f[i] = t.v[i].head<3>();
    //}
    auto ind = get_index(point.x(), point.y());
    if (not anti_alising || sampling_value == 0) 
    {
        if (z_interpolated > depth_buf[ind])
        {
            set_depth(point, z_interpolated);
            set_pixel(point, t.getColor());
        }
        return;
    }
    int insideTriCount = 0;
    float samplePointCounts = sampling_value * sampling_value * 1.0;
    float minZ_interpolated = std::numeric_limits<float>::infinity();
    for (int i = 0; i < sampling_value; i++)
    {
        for (int j = 0; j < sampling_value; j++)
        {
            ind = get_index(point.x(), point.y(), i, j);
            auto v = t.toVector4();
            float deltax = (i * 2.f + 1.f) / (sampling_value * 2.f);                 
            float deltay = (j * 2.f + 1.f) / (sampling_value * 2.f);
            Eigen::Vector3f sample_point = Eigen::Vector3f(point.x() - 0.5 + deltax, point.y() - 0.5 + deltay, 1.f);
            float _alpha, _beta, _gamma;
            std::tie(_alpha, _beta, _gamma) = computeBarycentric2D(sample_point.x(), sample_point.y(), v3f);
            float w_reciprocal = 1.0 / (_alpha / v[0].w() + _beta / v[1].w() + _gamma / v[2].w());
            float sample_z_interpolated = _alpha * v[0].z() / v[0].w() + _beta * v[1].z() / v[1].w() + _gamma * v[2].z() / v[2].w();
            sample_z_interpolated *= w_reciprocal;
            //std::cout << point.x() <<" ,"<< point.y() <<","<<i<<","<<j<<","<<ind<<","<<depth_buf.size()<< std::endl;
            if (insideTriangle(sample_point.x(), sample_point.y(), v3f) && sample_z_interpolated > depth_buf[ind])
            {
                insideTriCount++;
                set_depth(point, sample_z_interpolated, i, j);
            }
        }
    }
    //std::cout << insideTriCount << std::endl;
    //if (insideTriCount != 4 && insideTriCount != 0) std::cout << insideTriCount << std::endl;
    //Eigen::Vector3f resColor = t.getColor() * (insideTriCount * 1.0 / (sampling_value * sampling_value));
    //if (minZ_interpolated <= depth_buf[ind])
    {
        float n = (insideTriCount * 1.0 / (samplePointCounts));
        
        Eigen::Vector3f resColor = t.getColor() * n + frame_buf[get_index(point.x(), point.y())] ;
        //std::cout << n << " color:"<< resColor << std::endl;
        //set_depth(point, minZ_interpolated);
        set_pixel(point, resColor);
    }
    return;
}
void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{ 0, 0, 0 });
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), -1 * std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    if (not anti_alising)
    {
        depth_buf.resize(w * h);
    }
    else {
        depth_buf.resize(w * h * sampling_value * sampling_value);
    }
    
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height - 1 - y) * width + x;
}

int rst::rasterizer::get_index(int x, int y, int sx, int sy)
{
    int _x = x + sx * width;
    int _y = y + sy * height;
    int _width = width * sampling_value;
    int _height = height * sampling_value;
    return (_height - 1 - _y) * _width + _x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] = color;

}

void rst::rasterizer::set_depth(const Eigen::Vector3f& point, const float depth, const int sample_index_x, const int sample_index_y)
{
    if (not anti_alising)
    {
        auto ind = get_index(point.x(), point.y());
        depth_buf[ind] = depth;
    }
    else 
    {
        auto ind = get_index(point.x(),point.y(),sample_index_x,sample_index_y);
        depth_buf[ind] = depth;
    }
}

// clang-format on