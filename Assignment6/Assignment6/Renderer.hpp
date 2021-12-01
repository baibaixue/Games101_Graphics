//
// Created by goksu on 2/25/20.
//
#include "Scene.hpp"

#pragma once
// 碰撞点结构体
struct hit_payload
{
    float tNear;
    uint32_t index;
    Vector2f uv;
    Object* hit_obj;
};
// 渲染器
class Renderer
{
public:
    void Render(const Scene& scene);

private:
};
