#pragma once

#include <vector>
#include <memory>
#include "Vector.hpp"
#include "Object.hpp"
#include "Light.hpp"
// 场景
class Scene
{
public:
    // setting up options
    // 分辨率 1280 * 960
    int width = 1280;
    int height = 960;
    
    double fov = 90;
    // 背景颜色
    Vector3f backgroundColor = Vector3f(0.235294, 0.67451, 0.843137);   
    //最大深度，递归的最深次数
    int maxDepth = 5;
    // 一个极小值，ε
    float epsilon = 0.00001;

    Scene(int w, int h) : width(w), height(h)
    {}
    // 场景中添加物体和光源
    void Add(std::unique_ptr<Object> object) { objects.push_back(std::move(object)); }
    void Add(std::unique_ptr<Light> light) { lights.push_back(std::move(light)); }
    // nodiscard 标记： 保存返回值，
    // 取场景中的物体和光源
    [[nodiscard]] const std::vector<std::unique_ptr<Object> >& get_objects() const { return objects; }
    [[nodiscard]] const std::vector<std::unique_ptr<Light> >&  get_lights() const { return lights; }

private:
    // creating the scene (adding objects and lights)
    std::vector<std::unique_ptr<Object> > objects;// 场景中创建物体
    std::vector<std::unique_ptr<Light> > lights;    // 场景中创建光源
};