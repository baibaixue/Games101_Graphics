#pragma once
#include "Scene.hpp"
// ���ߺ������ཻʱ���صĽṹ
struct hit_payload
{
    float tNear;    // ���ߵĵ����ཻ����Ҫ�����ʱ�䣨��ʾ������룩
    uint32_t index;// �ཻ������������
    Vector2f uv;    // �ཻ��uv
    Object* hit_obj;    // �ཻ������
};

class Renderer
{
public:
    void Render(const Scene& scene);

private:
};