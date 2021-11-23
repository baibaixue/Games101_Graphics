#pragma once
#include "Scene.hpp"
// 光线和物体相交时返回的结构
struct hit_payload
{
    float tNear;    // 光线的到达相交点需要的最短时间（表示最近距离）
    uint32_t index;// 相交的三角形索引
    Vector2f uv;    // 相交点uv
    Object* hit_obj;    // 相交的物体
};

class Renderer
{
public:
    void Render(const Scene& scene);

private:
};