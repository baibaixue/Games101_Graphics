//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include "../header/rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_normals(const std::vector<Eigen::Vector3f>& normals)
{
    auto id = get_next_id();
    nor_buf.emplace(id, normals);

    normal_id = id;

    return {id};
}


// Bresenham's line drawing algorithm
// 布鲁斯汉姆划线方法
void rst::rasterizer::draw_line(Eigen::Vector3f begin, Eigen::Vector3f end)
{
    auto x1 = begin.x();
    auto y1 = begin.y();
    auto x2 = end.x();
    auto y2 = end.y();

    Eigen::Vector3f line_color = {255, 255, 255};

    int x,y,dx,dy,dx1,dy1,px,py,xe,ye,i;

    dx=x2-x1;
    dy=y2-y1;
    dx1=fabs(dx);
    dy1=fabs(dy);
    px=2*dy1-dx1;
    py=2*dx1-dy1;

    if(dy1<=dx1)
    {
        if(dx>=0)
        {
            x=x1;
            y=y1;
            xe=x2;
        }
        else
        {
            x=x2;
            y=y2;
            xe=x1;
        }
        Eigen::Vector2i point = Eigen::Vector2i(x, y);
        set_pixel(point,line_color);
        for(i=0;x<xe;i++)
        {
            x=x+1;
            if(px<0)
            {
                px=px+2*dy1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    y=y+1;
                }
                else
                {
                    y=y-1;
                }
                px=px+2*(dy1-dx1);
            }
//            delay(0);
            Eigen::Vector2i point = Eigen::Vector2i(x, y);
            set_pixel(point,line_color);
        }
    }
    else
    {
        if(dy>=0)
        {
            x=x1;
            y=y1;
            ye=y2;
        }
        else
        {
            x=x2;
            y=y2;
            ye=y1;
        }
        Eigen::Vector2i point = Eigen::Vector2i(x, y);
        set_pixel(point,line_color);
        for(i=0;y<ye;i++)
        {
            y=y+1;
            if(py<=0)
            {
                py=py+2*dx1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    x=x+1;
                }
                else
                {
                    x=x-1;
                }
                py=py+2*(dx1-dy1);
            }
//            delay(0);
            Eigen::Vector2i point = Eigen::Vector2i(x, y);
            set_pixel(point,line_color);
        }
    }
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}
// 判断平面内点(x,y)是否在空间内三角形在xy平面所在的投影内
static bool insideTriangle(float x, float y, const Vector4f* _v){
    /*
    Vector3f v[3];
    for(int i=0;i<3;i++)
        v[i] = {_v[i].x(),_v[i].y(), 1.0};
    Vector3f f0,f1,f2;
    f0 = v[1].cross(v[0]);
    f1 = v[2].cross(v[1]);
    f2 = v[0].cross(v[2]);
    Vector3f p(x,y,1.);
    if((p.dot(f0)*f0.dot(v[2])> 0) && (p.dot(f1)*f1.dot(v[0])> 0) && (p.dot(f2)*f2.dot(v[1])> 0))
        return true;
    return false;*/
    Eigen::Vector2f point = Eigen::Vector2f(x, y);
    float res[3];
    Eigen::Vector2f V[3];
    for (int i = 0; i < 3; i++)
    {
        V[i] = Eigen::Vector2f(_v[i].x(), _v[i].y());
    }
    for (int i = 0; i < 3; i++)
    {
        Eigen::Vector2f VP = point - V[i];
        Eigen::Vector2f VV = V[(i + 1) % 3] - V[i];
        res[i] = VP.x() * VV.y() * 1.0 - VP.y() * VV.x() * 1.0;
    }
    for (int i = 0; i < 3; i++)
    {
        int j = (i + 1) % 3;
        if (res[i] * res[j] < 0) return false;
    }
    return true;
}
// 重心坐标插值
static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector4f* v){
    float xa = v[0].x(), ya = v[0].y(), xb = v[1].x(), yb = v[1].y(), xc = v[2].x(), yc = v[2].y();
    float c1 = (-1 * (x - xb) * (yc - yb) + (y - yb) * (xc - xb)) / (-1 * (xa - xb) * (yc - yb) + (ya - yb) * (xc - xb));
    float c2 = (-1 * (x - xc) * (ya - yc) + (y - yc) * (xa - xc)) / (-1 * (xb - xc) * (ya - yc) + (yb - yc) * (xa - xc));
    float c3 = (1 - c1 - c2) < 0 ? 0 : (1 - c1 - c2);
    return {c1,c2,c3};
}
// 绘制三角形（先进行三角形的mvp变换和视口变化，后渲染三角形）
void rst::rasterizer::draw(std::vector<Triangle *> &TriangleList) {

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (const auto& t:TriangleList)
    {
        Triangle newtri = *t;
        // mm(vector4)，viewspace_pos(vector3)视图空间的三角形顶点坐标（只做模型变换和视图变化，不做投影变换）
        std::array<Eigen::Vector4f, 3> mm {
                (view * model * t->v[0]),
                (view * model * t->v[1]),
                (view * model * t->v[2])
        };

        std::array<Eigen::Vector3f, 3> viewspace_pos;

        std::transform(mm.begin(), mm.end(), viewspace_pos.begin(), [](auto& v) {
            return v.template head<3>();
        });
        // 三角形的mvp变换
        Eigen::Vector4f v[] = {
                mvp * t->v[0],
                mvp * t->v[1],
                mvp * t->v[2]
        };
        // 归一化
        //Homogeneous division
        for (auto& vec : v) {
            vec.x() /= vec.w();
            vec.y() /= vec.w();
            vec.z() /= vec.w();
        }
        // 这里是求经过mv变换后的空间三角形的法向量
        Eigen::Matrix4f inv_trans = (view * model).inverse().transpose();
        Eigen::Vector4f n[] = {
                inv_trans * to_vec4(t->normal[0], 0.0f),
                inv_trans * to_vec4(t->normal[1], 0.0f),
                inv_trans * to_vec4(t->normal[2], 0.0f)
        };
        //if (n[0].x() == 0 && n[0].y() == -0.0132841)
        // 视口变换
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            // 屏幕空间的顶点
            //screen space coordinates
            newtri.setVertex(i, v[i]);
        }
        // 视图空间的法线
        for (int i = 0; i < 3; ++i)
        {
            //view space normal
            newtri.setNormal(i, n[i].head<3>());
        }
        // 设置每个三角形顶点的颜色
        newtri.setColor(0, 148,121.0,92.0);
        newtri.setColor(1, 148,121.0,92.0);
        newtri.setColor(2, 148,121.0,92.0);

        // Also pass view space vertice position
        // 不仅需要传递需要光栅化的三角形信息，也需要传递三角形在的视图空间的顶点信息
        rasterize_triangle(newtri, viewspace_pos);
    }
}
// 3维向量插值
static Eigen::Vector3f interpolate(float alpha, float beta, float gamma, const Eigen::Vector3f& vert1, const Eigen::Vector3f& vert2, const Eigen::Vector3f& vert3, float weight)
{
    return (alpha * vert1 + beta * vert2 + gamma * vert3) * weight;
}
// 2维向量插值
static Eigen::Vector2f interpolate(float alpha, float beta, float gamma, const Eigen::Vector2f& vert1, const Eigen::Vector2f& vert2, const Eigen::Vector2f& vert3, float weight)
{
    auto u = (alpha * vert1[0] + beta * vert2[0] + gamma * vert3[0]);
    auto v = (alpha * vert1[1] + beta * vert2[1] + gamma * vert3[1]);

    u *= weight;
    v *= weight;

    return Eigen::Vector2f(u, v);
}

//Screen space rasterization
// 渲染空间三角形
void rst::rasterizer::rasterize_triangle(const Triangle& t, const std::array<Eigen::Vector3f, 3>& view_pos) 
{
    // TODO: From your HW3, get the triangle rasterization code.
    // TODO: Inside your rasterization loop:
    //    * v[i].w() is the vertex view space depth value z.
    //    * Z is interpolated view space depth for the current pixel
    //    * zp is depth between zNear and zFar, used for z-buffer
    auto v = t.v;
    // Bounding Box
    float minx = width, miny = height, maxx = 0, maxy = 0;
    for (auto i = 0; i < 3; i++)
    {

        minx = std::min(minx, v[i].x());
        miny = std::min(miny, v[i].y());
        maxx = std::max(maxx, v[i].x());
        maxy = std::max(maxy, v[i].y());
    }
    minx = std::floor(std::max(0.f, minx));
    miny = std::ceil(std::max(0.f, miny));
    maxx = std::floor(std::max(0.f, maxx));
    maxy = std::ceil(std::max(0.f, maxy));

    Eigen::MatrixXf TB(2, 3);
    Eigen::Matrix2f UV;
    Eigen::MatrixXf E(2, 3);
    Eigen::Vector3f E0 = view_pos[1] - view_pos[0];
    Eigen::Vector3f E1 = view_pos[2] - view_pos[0];
    float u0 = t.tex_coords[0].x(), u1 = t.tex_coords[1].x(), u2 = t.tex_coords[2].x();
    float v0 = t.tex_coords[0].y(), v1 = t.tex_coords[1].y(), v2 = t.tex_coords[2].y();
    float t1 = u1 - u0, b1 = v1 - v0, t2 = u2 - u0, b2 = v2 - v0;
    UV << b2, -b1, -t2, t1;
    E << E0.x(), E0.y(), E0.z(), E1.x(), E1.y(), E1.z();
    TB = 1.0 / (t1 * b2 - b1 * t2) * UV * E;
    Eigen::Vector3f tangent = TB.row(0);
    for (int w = minx; w < std::min((float)width, maxx + 1); w++)
    {
        for (int h = miny; h < std::min((float)height, maxy + 1); h++)
        {
            Eigen::Vector2f point = Eigen::Vector2f(w, h);
            if (insideTriangle(w*1.0, h*1.0, t.v))
            {
                float _alpha, _beta, _gamma;
                std::tie(_alpha, _beta, _gamma) = computeBarycentric2D(w, h, t.v);
                float z_interpolated = 1.0 / (_alpha / v[0].w() + _beta / v[1].w() + _gamma / v[2].w());
                int ind = get_index(w, h);
                if (depth_buf[ind] <= z_interpolated)
                {
                    depth_buf[ind] = z_interpolated;
                    float z0 = v[0].w(), z1 = v[1].w(), z2 = v[2].w();
                    Eigen::Vector3f zc = t.color[2] / z2, zn = t.normal[2] / z2, zp = view_pos[2] / z2;
                    Eigen::Vector2f zt = t.tex_coords[2] / z2;
                    auto interpolated_color = interpolate(_alpha,_beta,_gamma, t.color[0] / z0, t.color[1] / z1, zc ,z_interpolated);        // 颜色插值
                    auto interpolated_normal = interpolate(_alpha, _beta, _gamma, t.normal[0] /z0 , t.normal[1]/z1, zn, z_interpolated);      // 法向量插值
                    auto interpolated_texcoords = interpolate(_alpha, _beta, _gamma, t.tex_coords[0] / z0, t.tex_coords[1] / z1, zt, z_interpolated);   // uv坐标插值
                    auto interpolated_shadingcoords = interpolate(_alpha, _beta, _gamma, view_pos[0] / z0, view_pos[1] / z1, zp, z_interpolated);
                    fragment_shader_payload payload(interpolated_color, interpolated_normal, interpolated_texcoords, texture ? &*texture : nullptr);
                    payload.view_pos = interpolated_shadingcoords;
                    payload.tangent = tangent.normalized();
                    auto pixel_color = fragment_shader(payload);
                    set_pixel(Eigen::Vector2i(w,h), pixel_color);
                }
                //std::cout << z_interpolated << std::endl;
                //MSAA(z_interpolated, point, t);
            }
        }
    }
    // float Z = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    // float zp = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    // zp *= Z;

    // TODO: Interpolate the attributes:
    // auto interpolated_color
    // auto interpolated_normal
    // auto interpolated_texcoords
    // auto interpolated_shadingcoords

    // Use: fragment_shader_payload payload( interpolated_color, interpolated_normal.normalized(), interpolated_texcoords, texture ? &*texture : nullptr);
    // Use: payload.view_pos = interpolated_shadingcoords;
    // Use: Instead of passing the triangle's color directly to the frame buffer, pass the color to the shaders first to get the final color;
    // Use: auto pixel_color = fragment_shader(payload);

 
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
// 清空缓冲区
void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), - 1 * std::numeric_limits<float>::infinity());
    }
}
// 构造函数
rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);

    texture = std::nullopt;
}
// (x，y)的索引
int rst::rasterizer::get_index(int x, int y)
{
    return (height - y - 1) * width + x;
}
// 设置像素点颜色
void rst::rasterizer::set_pixel(const Vector2i &point, const Eigen::Vector3f &color)
{
    //old index: auto ind = point.y() + point.x() * width;
    int ind = (height-point.y() - 1)*width + point.x();
    frame_buf[ind] = color;
    //std::cout << "point:\n" << point << "\n\ncolor:\n" << color <<"\n" << std::endl;
}
// 设置顶点着色方法
void rst::rasterizer::set_vertex_shader(std::function<Eigen::Vector3f(vertex_shader_payload)> vert_shader)
{
    vertex_shader = vert_shader;
}
// 设置片元着色方法
void rst::rasterizer::set_fragment_shader(std::function<Eigen::Vector3f(fragment_shader_payload)> frag_shader)
{
    fragment_shader = frag_shader;
}

void rst::rasterizer::drawNormal(std::vector<Eigen::Vector3f> NewNormals, std::vector<Eigen::Vector3f> DrawViewPos)
{
    //std::cout << "DrawNormal" << DrawViewPos.size() << std::endl;
    for (int i = 0; i < DrawViewPos.size(); i++)
    {
        Eigen::Vector4f begin = to_vec4(DrawViewPos[i]);
        Eigen::Vector4f end = to_vec4(DrawViewPos[i] + 0.5 * NewNormals[i]);
        float f1 = (50 - 0.1) / 2.0;
        float f2 = (50 + 0.1) / 2.0;
        begin = projection * begin;
        end = projection * end;
        
        begin /= begin.w();
        end /= end.w();

        begin.x() = 0.5 * width * (begin.x() + 1.0);
        begin.y() = 0.5 * height * (begin.y() + 1.0);
        begin.z() = begin.z() * f1 + f2;

        end.x() = 0.5 * width * (end.x() + 1.0);
        end.y() = 0.5 * height * (end.y() + 1.0);
        end.z() = end.z() * f1 + f2;
        //std::cout << "begin:"<< begin << "\nend:" <<end<< std::endl;
        draw_line(begin.head<3>(), end.head<3>());
    }
}

