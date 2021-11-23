#pragma once

#include "Vector.hpp"
// ��Դ
class Light
{
public:
    Light(const Vector3f& p, const Vector3f& i)
        : position(p)
        , intensity(i)
    {}
    virtual ~Light() = default;
    Vector3f position;  //��Դλ��
    Vector3f intensity; // ��Դǿ��
};
