#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include "Vector.hpp"
#include "Renderer.hpp"
#include "Scene.hpp"
#include <optional>
// 角度制转弧度制
inline float deg2rad(const float &deg)
{ return deg * M_PI/180.0; }

// Compute reflection direction
// 反射光方向     
Vector3f reflect(const Vector3f &I, const Vector3f &N)
{
    return I - 2 * dotProduct(I, N) * N;
}

// [comment]
// （这里计算折射方向）
// Compute refraction direction using Snell's law使用斯奈尔定律计算折射方向
//需要处理两种情况
// 1. 光线在物体内部
// 2. 光线在物体外部
// We need to handle with care the two possible situations:
//
//    - When the ray is inside the object
//
//    - When the ray is outside.
//
// If the ray is outside, you need to make cosi positive cosi = -N.I
//
// If the ray is inside, you need to invert the refractive indices and negate the normal N
// [/comment]
Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior)
{
    // cosi为入射方向余弦值，cosi<0时，光线在物体外，cosi>0时，光线在物体内
    float cosi = clamp(-1, 1, dotProduct(I, N));
    //(光线在物体外时) 入射角正弦值：折射角正弦值=etat:etai
    float etai = 1, etat = ior;
    Vector3f n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);// 折射角余弦值的平方
    return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
}

// [comment]
// （这里计算反射后的能量）
// Compute Fresnel equation 计算菲涅尔方程
//
// \param I is the incident view direction I为入射方向
//
// \param N is the normal at the intersection point N为入射点的法线
//
// \param ior is the material refractive index ior为材质的折射系数
// [/comment]
float fresnel(const Vector3f &I, const Vector3f &N, const float &ior)
{
    float cosi = clamp(-1, 1, dotProduct(I, N));
    float etai = 1, etat = ior;
    if (cosi > 0) {  std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        return 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) / 2;
    }
    //由于能量守恒，折射能量占比+反射能量占比=1
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

// [comment]
// Returns true if the ray intersects an object, false otherwise. 如果光线和物体相交，返回true，否则返回false
//
// \param orig is the ray origin            （参数）orig为光线的起点
// \param dir is the ray direction          （参数）dir为光线的方向
// \param objects is the list of objects the scene contains （参数）objects为场景中包含的所有物体
// \param[out] tNear contains the distance to the cloesest intersected object.  （输出） tNear ：光线达到最近物体的距离
// \param[out] index stores the index of the intersect triangle if the interesected object is a mesh. （输出）index : 如果相交的物体是一个三角形网格时，index返回相交三角形的索引
// \param[out] uv stores the u and v barycentric coordinates of the intersected point   （输出） uv:uv储存相交点的重心坐标
// \param[out] *hitObject stores the pointer to the intersected object (used to retrieve material information, etc.) （输出）*hitObject : 用于储存相交的一系列物体，用于获取材质信息等。
// \param isShadowRay is it a shadow ray. We can return from the function sooner as soon as we have found a hit. （返回）isShadowRay :是否阴影
// [/comment]
std::optional<hit_payload> trace(
        const Vector3f &orig, const Vector3f &dir,
        const std::vector<std::unique_ptr<Object> > &objects)
{
    float tNear = kInfinity;
    std::optional<hit_payload> payload;
    for (const auto & object : objects)
    {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (object->intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear)
        {
            payload.emplace();
            payload->hit_obj = object.get();
            payload->tNear = tNearK;
            payload->index = indexK;
            payload->uv = uvK;
            tNear = tNearK;
        }
    }

    return payload;
}

// [comment]
// Implementation of the Whitted-style light transport algorithm (E [S*] (D|G) L)
// Wihtted风格光线追踪算法的实现
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
// 递归调用求光线与物体交点处的颜色
// If the material of the intersected object is either reflective or reflective and refractive,
// then we compute the reflection/refraction direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refraction depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is diffuse/glossy we use the Phong illumation model to compute the color
// at the intersection point.
// 如果物体的材质为反射或反射和折射，计算反射和折射的方向，并递归调用castRay方法，
// 如果物体表面材质透明时，我们通过菲涅尔方程的结果混合反射和折射的颜色
// 如果物体表面漫反射/有光泽， 我们运用phong模型计算交点处的颜色
// （参数）depth : 递归深度
// [/comment]
Vector3f castRay(
        const Vector3f &orig, const Vector3f &dir, const Scene& scene,
        int depth)
{
    if (depth > scene.maxDepth) {
        return Vector3f(0.0,0.0,0.0);
    }

    Vector3f hitColor = scene.backgroundColor;
    if (auto payload = trace(orig, dir, scene.get_objects()); payload)
    {
        Vector3f hitPoint = orig + dir * payload->tNear;// 相交处点的坐标
        Vector3f N; // normal
        Vector2f st; // st coordinates
        payload->hit_obj->getSurfaceProperties(hitPoint, dir, payload->index, payload->uv, N, st);
        switch (payload->hit_obj->materialType) {
            case REFLECTION_AND_REFRACTION:// 折射和反射的情况
            {
                Vector3f reflectionDirection = normalize(reflect(dir, N));// 反射光方向
                Vector3f refractionDirection = normalize(refract(dir, N, payload->hit_obj->ior));// 折射光方向
                Vector3f reflectionRayOrig = (dotProduct(reflectionDirection,   N) < 0) ?// 反射光起点
                                             hitPoint - N * scene.epsilon :
                                             hitPoint + N * scene.epsilon;
                Vector3f refractionRayOrig = (dotProduct(refractionDirection, N) < 0) ?// 折射光起点
                                             hitPoint - N * scene.epsilon :
                                             hitPoint + N * scene.epsilon;
                Vector3f reflectionColor = castRay(reflectionRayOrig, reflectionDirection, scene, depth + 1);// 反射光线递归
                Vector3f refractionColor = castRay(refractionRayOrig, refractionDirection, scene, depth + 1);// 折射光线递归
                float kr = fresnel(dir, N, payload->hit_obj->ior);// 根据菲涅尔定律计算反射光能量率，折射光能量率为（1-kr）
                hitColor = reflectionColor * kr + refractionColor * (1 - kr);
                break;
            }
            case REFLECTION:// 只有反射的情况
            {
                float kr = fresnel(dir, N, payload->hit_obj->ior);// 反射光能量
                Vector3f reflectionDirection = reflect(dir, N);// 反射光方向
                Vector3f reflectionRayOrig = (dotProduct(reflectionDirection, N) < 0) ?// 反射光起点
                                             hitPoint + N * scene.epsilon :
                                             hitPoint - N * scene.epsilon;
                hitColor = castRay(reflectionRayOrig, reflectionDirection, scene, depth + 1) * kr;
                break;
            }
            default:// 默认运用phong模型计算光照
            {
                // [comment]
                // We use the Phong illumation model int the default case. The phong model
                // is composed of a diffuse and a specular reflection component.
                // 运行漫反射系数和镜面反射系数计算phong模型
                // [/comment]
                Vector3f lightAmt = 0, specularColor = 0;
                Vector3f shadowPointOrig = (dotProduct(dir, N) < 0) ?
                                           hitPoint + N * scene.epsilon :
                                           hitPoint - N * scene.epsilon;
                // [comment]
                // Loop over all lights in the scene and sum their contribution up
                // We also apply the lambert cosine law
                // 遍历所有光源进行着色，利用了lambert余弦定理
                // [/comment]
                for (auto& light : scene.get_lights()) {
                    Vector3f lightDir = light->position - hitPoint;
                    // square of the distance between hitPoint and the light
                    float lightDistance2 = dotProduct(lightDir, lightDir);// 光源与相交点距离的平方
                    lightDir = normalize(lightDir);// 光照方向
                    float LdotN = std::max(0.f, dotProduct(lightDir, N));//光照方向点乘法线方向
                    // is the point in shadow, and is the nearest occluding object closer to the object than the light itself?
                    // 判断点是否在阴影中，逆着光线的方向寻找，如果发现光线有相交点并且相交点到光源的距离平方小于当前点到光源距离的平方，证明在阴影中，否则不在阴影中
                    auto shadow_res = trace(shadowPointOrig, lightDir, scene.get_objects());
                    bool inShadow = shadow_res && (shadow_res->tNear * shadow_res->tNear < lightDistance2);
                    // 如果在阴影中，不进行着色，否则进行漫反射的着色
                    lightAmt += inShadow ? 0 : light->intensity * LdotN;
                    Vector3f reflectionDirection = reflect(-lightDir, N);
                    // 这里计算镜面反射的情况
                    specularColor += powf(std::max(0.f, -dotProduct(reflectionDirection, dir)),
                        payload->hit_obj->specularExponent) * light->intensity;
                }

                hitColor = lightAmt * payload->hit_obj->evalDiffuseColor(st) * payload->hit_obj->Kd + specularColor * payload->hit_obj->Ks;
                break;
            }
        }
    }

    return hitColor;
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = std::tan(deg2rad(scene.fov * 0.5f));
    float imageAspectRatio = scene.width / (float)scene.height;

    // Use this variable as the eye position to start your rays.
    Vector3f eye_pos(0);
    int m = 0;
    for (float j = 0.f; j < scene.height; j = j + 1.f)
    {
        for (float i = 0.f; i < scene.width; i = i + 1.f)
        {
            // generate primary ray direction
            //float x = (float)(i + 0.5f) / scene.width * 2.f -1.f;
            //float y = (float)(scene.height - j - 0.5f) / scene.height * 2.f - 1.f;
            //x *= scale;
            //y *= scale;
            //x *= imageAspectRatio;
            //Vector3f dir = Vector3f(x, y, -1);

            // TODO: Find the x and y positions of the current pixel to get the direction
            // vector that passes through it.
            // Also, don't forget to multiply both of them with the variable *scale*, and
            // x (horizontal) variable with the *imageAspectRatio*            
            Vector3f r = (((i + 0.5f) / scene.width * 2.f - 1.f) * scale) * imageAspectRatio * Vector3f(1.f, 0.f, 0.f);
            Vector3f u = ((1.f - (j + 0.5f) / scene.height) * 2.f - 1.f) * scale * Vector3f(0.f, 1.f, 0.f);
            Vector3f dir = Vector3f(0, 0, -1) + r + u;

            dir = normalize(dir);
            framebuffer[m++] = 255 *  castRay(eye_pos, dir, scene, 0);
        }
        UpdateProgress(j / (float)scene.height);
    }

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (char)(clamp(0, 255, framebuffer[i].x));
        color[1] = (char)(clamp(0, 255, framebuffer[i].y));
        color[2] = (char)(clamp(0, 255, framebuffer[i].z));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);   
    int key = 0;

    while (key != 27)
    {
        std::string filename = "output.png";
        cv::Mat image(scene.height, scene.width, CV_32FC3, framebuffer.data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imshow("image", image);
        cv::imwrite(filename, image);
        key = cv::waitKey(10);
    }
}
