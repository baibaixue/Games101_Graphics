#pragma warning(disable:4996)
#include "Triangle.h"
#include "rasterizer.h"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;
// ��ͼ��������ֻ����ƽ�Ʊ任������ͼ������eye_pos��ƽ����ԭ��
Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    //4X4��λ����
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}
// ģ�;���rotation_angleΪģ����z����ת�ĽǶ�(�Ƕ���)
Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    Eigen::Matrix4f Rotate;
    float angle = rotation_angle / 180.0 * MY_PI;
    Rotate << std::cos(angle), -std::sin(angle), 0, 0,
        std::sin(angle), std::cos(angle), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

    model = Rotate * model;

    return model;
}
//�������ԭ�������ת�任����axis:��ԭ�����ķ��� angle: �Ƕ�
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

// ͶӰ����eye_fov: ��Ұ�� aspect_ratio �� ��߱� �� zNear : ��ƽ����� �� zFar �� Զƽ����� (���벻�����꣬�����������ΪZ�Ḻ���ᣬ������������Ӧ��Ϊ��)
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
    float zNear, float zFar)
{
    // Students will implement this function

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.
    //��͸��tͶӰת��Ϊ����ͶӰ
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
    float angle = -40;
    bool command_line = false;
    std::string filename = "output.png";
    // ����������3ʱ������ͼƬ������
    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default stof()�������ַ���ת��Ϊ������
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }
    // ���ڴ�С700*700
    rst::rasterizer r(700, 700);
    // ��������꣨0��0��5��
    Eigen::Vector3f eye_pos = { 0, 0, 5 };
    // �����ε���������
    std::vector<Eigen::Vector3f> pos{ {2.0,0.0,-3.0},{0.0,2.0,-1.0},{-2.0,0.0,-2.0} };
    // ��Ⱦ˳��
    std::vector<Eigen::Vector3i> ind{ {0, 1, 2} };

    // ת����
    Eigen::Vector3f axis = { 1, 1 ,1 };
    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        // ��ʼ����Ȼ�������֡������
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        // ��Ұ45�ȣ��ݺ��Ϊ1 ����ƽ�����Ϊ0.1 , Զƽ�����Ϊ50
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));
        //����������
        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        //opencv����ͼ��
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        //r.set_model(get_model_matrix(angle));
        r.set_model(get_rotation(axis, angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}

