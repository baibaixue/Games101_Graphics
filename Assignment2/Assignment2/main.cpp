// clang-format off
#include <iostream>
#include <opencv2/opencv.hpp>
#include "rasterizer.h"
#include "global.h"
#include "Triangle.h"
// constexpr c++11中新增关键词，在编译时计算，避免运行时计算的消耗
constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0],
        0, 1, 0, -eye_pos[1],
        0, 0, 1, -eye_pos[2],
        0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle, Eigen::Vector3f& scale)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f Rotate;
    float angle = rotation_angle / 180.0 * MY_PI;
    Rotate << std::cos(angle), -std::sin(angle), 0, 0,
        std::sin(angle), std::cos(angle), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

    Eigen::Matrix4f Scale;
    Scale << scale.x(), 0, 0, 0,
        0, scale.y(), 0, 0,
        0, 0, scale.z(), 0,
        0, 0, 0, 1;

    model = Rotate * Scale * model;
    return model;
}
//绕任意过原点的轴旋转变换矩阵，axis:过原点的轴的方向， angle: 角度
Eigen::Matrix4f get_rotation(Eigen::Vector3f axis, float angle)
{
    Eigen::Matrix4f Model = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f Identity = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f dualMatrix;

    Eigen::Vector4f axis4f = Eigen::Vector4f(axis.x(), axis.y(), axis.z(), 1.f);

    dualMatrix << 0, -axis.x(), axis.y(), 0,
        axis.z(), 0, -axis.x(), 0,
        -axis.y(), axis.x(), 0, 0,
        0, 0, 0, 1;
    Eigen::Matrix4f Rotation = std::cos(angle / 180.0 * MY_PI) * Identity + (1 - std::cos(angle / 180.0 * MY_PI)) * axis4f * axis4f.transpose() + std::sin(angle / 180.0 * MY_PI) * dualMatrix;

    Model = Rotation * Model;

    return Model;
}
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    // TODO: Copy-paste your implementation from the previous assignment.
    Eigen::Matrix4f projection;
    //将透视t投影转变为正交投影
    float fov_angle = (eye_fov / 180.0 * MY_PI) / 2;
    float zn = -zNear;
    float zf = -zFar;
    Eigen::Matrix4f presp;
    presp << -1 / (aspect_ratio * std::tan(fov_angle)), 0, 0, 0,
        0, -1 / std::tan(fov_angle), 0, 0,
        0, 0, (zn + zf) / (zn - zf), -2 * zn * zf / (zn - zf),
        0, 0, 1, 0;
    projection = presp * projection;
    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc == 2)
    {
        command_line = true;
        filename = std::string(argv[1]);
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = { 0,0,5 };

    Eigen::Vector3f Scale = { 2.f,2.f,2.f };
    // 转动轴
    Eigen::Vector3f axis = { 0, 0 ,1 };
    std::vector<Eigen::Vector3f> pos
    {
        {2, 0, -2},
        {0, 2, -10},
        {-2, 0, -2},
        {3.5, -1, -2},
        {2.5, 1.5, -5},
        {-1, 0.5, -5}
    };
    // 图元顶点集合，这里表示两个图元（三角形），每个图元三个顶点
    std::vector<Eigen::Vector3i> ind
    {
            {0, 1, 2},
            {3, 4, 5}
    };

    std::vector<Eigen::Vector3f> cols
    {
            {217.0, 238.0, 185.0},
            {217.0, 238.0, 185.0},
            {217.0, 238.0, 185.0},
            {185.0, 217.0, 238.0},
            {185.0, 217.0, 238.0},
            {185.0, 217.0, 238.0}
    };

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);
    auto col_id = r.load_colors(cols);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle, Scale));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        //r.set_model(get_rotation(axis, angle));
        r.set_model(get_model_matrix(angle, Scale));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        //std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') angle += 10;
        else if (key == 'd') angle -= 10;

        //std::cout << angle << std::endl;
    }

    return 0;
}
// clang-format on