//
// Created by Göksu Güvendiren on 2019-05-14.
//

#pragma once

#include "Vector.hpp"
//光线
class Light
{
public:
    Light(const Vector3f &p, const Vector3f &i) : position(p), intensity(i) {}
    virtual ~Light() = default;
    // 光源位置
    Vector3f position;
    // 光照强度
    Vector3f intensity;
};
