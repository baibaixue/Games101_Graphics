#pragma once

#include "Vector.hpp"
// 光源
class Light
{
public:
    Light(const Vector3f& p, const Vector3f& i)
        : position(p)
        , intensity(i)
    {}
    virtual ~Light() = default;
    Vector3f position;  //光源位置
    Vector3f intensity; // 光源强度
};
