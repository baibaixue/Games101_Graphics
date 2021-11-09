//
// Created by goksu on 4/6/19.
//
// rasterizer : ��դ��
#pragma once //����һ���ļ������ö��

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
    //���Ƶ�ͼԪ���ͣ��߻�������
    enum class Primitive
    {
        Line,
        Triangle
    };

    /*
     * For the curious : The draw function takes two buffer id's as its arguments.
     * These two structs make sure that if you mix up with their orders, the
     * compiler won't compile it. Aka : Type safety
     * ���ƺ���ʹ��������������id��Ϊ�������������ṹ�屣֤��������������ߵ�˳�򣬱������ᱨ��
     * ��Ϊ�����ͱ�����
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
        // λ�û�����
        pos_buf_id load_positions(const std::vector<Eigen::Vector3f>& positions);
        // ����������(��Ⱦ˳��)
        ind_buf_id load_indices(const std::vector<Eigen::Vector3i>& indices);
        // ���ڲ���ģ�;�����Ϊ�������ݸ���դ����
        void set_model(const Eigen::Matrix4f& m);
        // ����ͼ�任��������Ϊ�ڲ���ͼ����
        void set_view(const Eigen::Matrix4f& v);
        // ���ڲ���ͶӰ������Ϊ��������p,�����ݸ���դ����
        void set_projection(const Eigen::Matrix4f& p);

        //����Ļ���ص㣨x,y����Ϊ��r,g,b������ɫ�� ��д����Ӧ��֡������λ��
        void set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color);

        // ��ջ�����
        void clear(Buffers buff);
        //���ƻ���������
        void draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, Primitive type);

        std::vector<Eigen::Vector3f>& frame_buffer() { return frame_buf; }

    private:
        // ����
        void draw_line(Eigen::Vector3f begin, Eigen::Vector3f end);
        // ��դ��һ�������ε�դ�򣬼����λ���
        void rasterize_wireframe(const Triangle& t);

    private:
        Eigen::Matrix4f model;  // ģ�;���
        Eigen::Matrix4f view;   // ��ͼ����
        Eigen::Matrix4f projection; // ͶӰ����
        // λ�û�����
        std::map<int, std::vector<Eigen::Vector3f>> pos_buf;
        std::map<int, std::vector<Eigen::Vector3i>> ind_buf;

        // ���ƴ����е�ÿһ�����ض���Ҫ��ɫ���ݺ��������
        // ֡����������ڴ洢��Ҫ����Ļ�ϻ��Ƶ���ɫ����
        std::vector<Eigen::Vector3f> frame_buf ;
        // ��Ȼ�����
        std::vector<float> depth_buf ;
        int get_index(int x, int y);

        // ��Ҫ���ƵĴ��ڵĿ�͸�
        int width, height;

        int next_id = 0;
        int get_next_id() { return next_id++; }
    };
} // namespace rst
