//
// Created by LEI XU on 4/28/19.
//
//
// This loader is created by Robert Smith.
// https://github.com/Bly7/OBJ-Loader
// Use the MIT license.


#ifndef RASTERIZER_OBJ_LOADER_H
#define RASTERIZER_OBJ_LOADER_H

#pragma once


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>

// Print progress to console while loading (large models)
// 加载时将进度打印到控制台（对于大型模型）
#define OBJL_CONSOLE_OUTPUT

// Namespace: OBJL
//
// Description: The namespace that holds eveyrthing that
//	is needed and used for the OBJ Model Loader
namespace objl
{
    // Structure: Vector2 2维向量结构
    //
    // Description: A 2D Vector that Holds Positional Data保存位置数据的二维向量
    struct Vector2
    {
        // Default Constructor默认构造函数
        Vector2()
        {
            X = 0.0f;
            Y = 0.0f;
        }
        // Variable Set Constructor复制构造函数
        Vector2(float X_, float Y_)
        {
            X = X_;
            Y = Y_;
        }
        // Bool Equals Operator Overload 运算符重载
        bool operator==(const Vector2& other) const
        {
            return (this->X == other.X && this->Y == other.Y);
        }
        // Bool Not Equals Operator Overload
        bool operator!=(const Vector2& other) const
        {
            return !(this->X == other.X && this->Y == other.Y);
        }
        // Addition Operator Overload
        Vector2 operator+(const Vector2& right) const
        {
            return Vector2(this->X + right.X, this->Y + right.Y);
        }
        // Subtraction Operator Overload
        Vector2 operator-(const Vector2& right) const
        {
            return Vector2(this->X - right.X, this->Y - right.Y);
        }
        // Float Multiplication Operator Overload
        Vector2 operator*(const float& other) const
        {
            return Vector2(this->X *other, this->Y * other);
        }

        // Positional Variables
        float X;
        float Y;
    };

    // Structure: Vector3
    //
    // Description: A 3D Vector that Holds Positional Data 保存位置信息的三维向量
    struct Vector3
    {
        // Default Constructor
        Vector3()
        {
            X = 0.0f;
            Y = 0.0f;
            Z = 0.0f;
        }
        // Variable Set Constructor
        Vector3(float X_, float Y_, float Z_)
        {
            X = X_;
            Y = Y_;
            Z = Z_;
        }
        // Bool Equals Operator Overload
        bool operator==(const Vector3& other) const
        {
            return (this->X == other.X && this->Y == other.Y && this->Z == other.Z);
        }
        // Bool Not Equals Operator Overload
        bool operator!=(const Vector3& other) const
        {
            return !(this->X == other.X && this->Y == other.Y && this->Z == other.Z);
        }
        // Addition Operator Overload
        Vector3 operator+(const Vector3& right) const
        {
            return Vector3(this->X + right.X, this->Y + right.Y, this->Z + right.Z);
        }
        // Subtraction Operator Overload
        Vector3 operator-(const Vector3& right) const
        {
            return Vector3(this->X - right.X, this->Y - right.Y, this->Z - right.Z);
        }
        // Float Multiplication Operator Overload
        Vector3 operator*(const float& other) const
        {
            return Vector3(this->X * other, this->Y * other, this->Z * other);
        }
        // Float Division Operator Overload
        Vector3 operator/(const float& other) const
        {
            return Vector3(this->X / other, this->Y / other, this->Z / other);
        }

        // Positional Variables
        float X;
        float Y;
        float Z;
    };

    // Structure: Vertex 顶点
    //
    // Description: Model Vertex object that holds
    //	a Position, Normal, and Texture Coordinate
    // 对保存顶点，法线和纹理坐标的顶点进行建模
    struct Vertex
    {
        // Position Vector 位置向量
        Vector3 Position;

        // Normal Vector 法向量
        Vector3 Normal;

        // Texture Coordinate Vector 纹理坐标向量
        Vector2 TextureCoordinate;
    };
    // 材质属性
    struct Material
    {
        Material()
        {
            name;
            Ns = 0.0f;
            Ni = 0.0f;
            d = 0.0f;
            illum = 0;
        }

        // Material Name 材质名称
        std::string name;
        // Ambient Color 环境光系数
        Vector3 Ka;
        // Diffuse Color 漫反射系数
        Vector3 Kd;
        // Specular Color 镜面反射系数
        Vector3 Ks;
        // Specular Exponent 镜面反射指数
        float Ns;
        // Optical Density 光密度
        float Ni;
        // Dissolve 溶解度？
        float d;
        // Illumination 光照强度
        int illum;
        // Ambient Texture Map 环境光纹理贴图
        std::string map_Ka;
        // Diffuse Texture Map 漫反射纹理贴图
        std::string map_Kd;
        // Specular Texture Map 镜面反射纹理贴图
        std::string map_Ks;
        // Specular Hightlight Map 镜面反射高光贴图
        std::string map_Ns;
        // Alpha Texture Map alpha纹理贴图
        std::string map_d;
        // Bump Map 凹凸贴图
        std::string map_bump;
    };

    // Structure: Mesh 网格
    //
    // Description: A Simple Mesh Object that holds
    //	a name, a vertex list, and an index list
    // 包含名称，顶点列表和索引列表的简单网格对象
    struct Mesh
    {
        // Default Constructor
        Mesh()
        {

        }
        // Variable Set Constructor
        Mesh(std::vector<Vertex>& _Vertices, std::vector<unsigned int>& _Indices)
        {
            Vertices = _Vertices;
            Indices = _Indices;
        }
        // Mesh Name 网格名称
        std::string MeshName;
        // Vertex List 顶点列表
        std::vector<Vertex> Vertices;
        // Index List 索引列表
        std::vector<unsigned int> Indices;

        // Material 网格材质
        Material MeshMaterial;
    };

    // Namespace: Math
    //
    // Description: The namespace that holds all of the math
    //	functions need for OBJL 
    // 包含object_loader(OBJL)所需的所有数学函数的命名空间
    namespace math
    {
        // Vector3 Cross Product 三维向量叉乘
        Vector3 CrossV3(const Vector3 a, const Vector3 b)
        {
            return Vector3(a.Y * b.Z - a.Z * b.Y,
                           a.Z * b.X - a.X * b.Z,
                           a.X * b.Y - a.Y * b.X);
        }

        // Vector3 Magnitude Calculation 三位向量求模运算
        float MagnitudeV3(const Vector3 in)
        {
            return (sqrtf(powf(in.X, 2) + powf(in.Y, 2) + powf(in.Z, 2)));
        }

        // Vector3 DotProduct 三维向量点乘运算
        float DotV3(const Vector3 a, const Vector3 b)
        {
            return (a.X * b.X) + (a.Y * b.Y) + (a.Z * b.Z);
        }

        // Angle between 2 Vector3 Objects 两个三维向量之间的夹角
        float AngleBetweenV3(const Vector3 a, const Vector3 b)
        {
            float angle = DotV3(a, b);
            angle /= (MagnitudeV3(a) * MagnitudeV3(b));
            return angle = acosf(angle);
        }

        // Projection Calculation of a onto b 向量a在向量b上的投影计算
        Vector3 ProjV3(const Vector3 a, const Vector3 b)
        {
            Vector3 bn = b / MagnitudeV3(b);
            return bn * DotV3(a, bn);
        }
    }

    // Namespace: Algorithm算法
    //
    // Description: The namespace that holds all of the
    // Algorithms needed for OBJL 包含所有OBJL所需的算法命名空间
    namespace algorithm
    {
        // Vector3 Multiplication Opertor Overload 三维向量数乘的运算符重载
        Vector3 operator*(const float& left, const Vector3& right)
        {
            return Vector3(right.X * left, right.Y * left, right.Z * left);
        }

        // A test to see if P1 is on the same side as P2 of a line segment ab
        // 测试P1，P2 两个点是否在线段ab的同一侧
        bool SameSide(Vector3 p1, Vector3 p2, Vector3 a, Vector3 b)
        {
            Vector3 cp1 = math::CrossV3(b - a, p1 - a);
            Vector3 cp2 = math::CrossV3(b - a, p2 - a);

            if (math::DotV3(cp1, cp2) >= 0)
                return true;
            else
                return false;
        }

        // Generate a cross produect normal for a triangle
        // 生成三角的面的法向量
        Vector3 GenTriNormal(Vector3 t1, Vector3 t2, Vector3 t3)
        {
            Vector3 u = t2 - t1;
            Vector3 v = t3 - t1;

            Vector3 normal = math::CrossV3(u,v);

            return normal;
        }

        // Check to see if a Vector3 Point is within a 3 Vector3 Triangle
        // 判断一个空间点是否在空间三角形的内部
        bool inTriangle(Vector3 point, Vector3 tri1, Vector3 tri2, Vector3 tri3)
        {
            // Test to see if it is within an infinite prism that the triangle outlines.
            // 检测点point是否位于空间三角形的棱柱内
            bool within_tri_prisim = SameSide(point, tri1, tri2, tri3) && SameSide(point, tri2, tri1, tri3)
                                     && SameSide(point, tri3, tri1, tri2);

            // If it isn't it will never be on the triangle
            // 不在三棱柱内，就不会位于三角形内
            if (!within_tri_prisim)
                return false;

            // Calulate Triangle's Normal
            // 计算三角形的法向量
            Vector3 n = GenTriNormal(tri1, tri2, tri3);

            // Project the point onto this normal
            // 计算点point在法向量n上的投影
            Vector3 proj = math::ProjV3(point, n);

            // If the distance from the triangle to the point is 0
            //	it lies on the triangle
            // 如果投影模长为0，点p位于三角形所在的平面上，空间点p在空间三角形的内部
            if (math::MagnitudeV3(proj) == 0)
                return true;
            else
                return false;
        }

        // Split a String into a string array at a given token
        // 在给定的标记token处，将字符串拆分为字符串数组（给定的字符串in，输出的字符串out）
        inline void split(const std::string &in,
                          std::vector<std::string> &out,
                          std::string token)
        {
            out.clear();

            std::string temp;

            for (int i = 0; i < int(in.size()); i++)
            {
                std::string test = in.substr(i, token.size());

                if (test == token)
                {
                    if (!temp.empty())
                    {
                        out.push_back(temp);
                        temp.clear();
                        i += (int)token.size() - 1;
                    }
                    else
                    {
                        out.push_back("");
                    }
                }
                else if (i + token.size() >= in.size())
                {
                    temp += in.substr(i, token.size());
                    out.push_back(temp);
                    break;
                }
                else
                {
                    temp += in[i];
                }
            }
        }

        // Get tail of string after first token and possibly following spaces
        // 获取第一个标记（token）后的字符串尾部，后可能跟空格
        inline std::string tail(const std::string &in)
        {
            size_t token_start = in.find_first_not_of(" \t");
            size_t space_start = in.find_first_of(" \t", token_start);
            size_t tail_start = in.find_first_not_of(" \t", space_start);
            size_t tail_end = in.find_last_not_of(" \t");
            if (tail_start != std::string::npos && tail_end != std::string::npos)
            {
                return in.substr(tail_start, tail_end - tail_start + 1);
            }
            else if (tail_start != std::string::npos)
            {
                return in.substr(tail_start);
            }
            return "";
        }

        // Get first token of string
        // 得到字符串的第一个标记（token）处
        inline std::string firstToken(const std::string &in)
        {
            if (!in.empty())
            {
                size_t token_start = in.find_first_not_of(" \t");
                size_t token_end = in.find_first_of(" \t", token_start);
                if (token_start != std::string::npos && token_end != std::string::npos)
                {
                    return in.substr(token_start, token_end - token_start);
                }
                else if (token_start != std::string::npos)
                {
                    return in.substr(token_start);
                }
            }
            return "";
        }

        // Get element at given index position
        // 在给定的索引处获取元素
        template <class T>
        inline const T & getElement(const std::vector<T> &elements, std::string &index)
        {
            int idx = std::stoi(index); // stoi:数字转字符串方法
            if (idx < 0)
                idx = int(elements.size()) + idx;
            else
                idx--;
            return elements[idx];
        }
    }

    // Class: Loader
    //
    // Description: The OBJ Model Loader 3D模型加载器
    class Loader
    {
    public:
        // Default Constructor
        Loader()
        {

        }
        ~Loader()
        {
            LoadedMeshes.clear();
        }

        // Load a file into the loader 加载文件到加载器（Loader）中
        //
        // If file is loaded return true 加载成功返回true
        //
        // If the file is unable to be found
        // or unable to be loaded return false 
        // 如果找不到文件或者加载失败返回false
        bool LoadFile(std::string Path)
        {
            // If the file is not an .obj file return false
            if (Path.substr(Path.size() - 4, 4) != ".obj")
                return false;


            std::ifstream file(Path);

            if (!file.is_open())
                return false;

            // 清空一下加载器
            LoadedMeshes.clear();
            LoadedVertices.clear();
            LoadedIndices.clear();

            std::vector<Vector3> Positions; // 位置
            std::vector<Vector2> TCoords;   // 纹理坐标
            std::vector<Vector3> Normals;   // 法向量

            std::vector<Vertex> Vertices;   // 顶点
            std::vector<unsigned int> Indices;  // 索引

            std::vector<std::string> MeshMatNames;  // 网格名称

            bool listening = false;
            std::string meshname;

            Mesh tempMesh;

#ifdef OBJL_CONSOLE_OUTPUT
            const unsigned int outputEveryNth = 1000; // 每读取1000行打印一次
            unsigned int outputIndicator = outputEveryNth;
#endif

            std::string curline;
            while (std::getline(file, curline))
            {
#ifdef OBJL_CONSOLE_OUTPUT
                if ((outputIndicator = ((outputIndicator + 1) % outputEveryNth)) == 1)
                {
                    if (!meshname.empty())
                    {
                        // 打印读取的定点数，纹理坐标数，法向量个数，三角形数，材质不为空时，打印最后一个材质的材质名
                        std::cout
                                << "\r- " << meshname
                                << "\t| vertices > " << Positions.size()
                                << "\t| texcoords > " << TCoords.size()
                                << "\t| normals > " << Normals.size()
                                << "\t| triangles > " << (Vertices.size() / 3)
                                << (!MeshMatNames.empty() ? "\t| material: " + MeshMatNames.back() : "");
                    }
                }
#endif

                // Generate a Mesh Object or Prepare for an object to be created
                // 生成网格对象或准备要创建的对象
                if (algorithm::firstToken(curline) == "o" || algorithm::firstToken(curline) == "g" || curline[0] == 'g')
                {
                    // 读网格对象的名称，以“o”或“g”开头 后 接网格对象的名称
                    if (!listening)
                    {
                        listening = true;

                        if (algorithm::firstToken(curline) == "o" || algorithm::firstToken(curline) == "g")
                        {
                            meshname = algorithm::tail(curline);
                        }
                        else
                        {
                            meshname = "unnamed";
                        }
                    }
                    else
                    {
                        // Generate the mesh to put into the array
                        // 生成要放入array的网格
                        if (!Indices.empty() && !Vertices.empty())
                        {
                            // Create Mesh
                            tempMesh = Mesh(Vertices, Indices);
                            tempMesh.MeshName = meshname;

                            // Insert Mesh
                            LoadedMeshes.push_back(tempMesh);

                            // Cleanup
                            Vertices.clear();
                            Indices.clear();
                            meshname.clear();

                            meshname = algorithm::tail(curline);
                        }
                        else
                        {
                            if (algorithm::firstToken(curline) == "o" || algorithm::firstToken(curline) == "g")
                            {
                                meshname = algorithm::tail(curline);
                            }
                            else
                            {
                                meshname = "unnamed";
                            }
                        }
                    }
#ifdef OBJL_CONSOLE_OUTPUT
                    std::cout << std::endl;
                    outputIndicator = 0;
#endif
                }
                // Generate a Vertex Position
                // v : 顶点位置(vector3) (vertex position)
                if (algorithm::firstToken(curline) == "v")
                {
                    std::vector<std::string> spos;
                    Vector3 vpos;
                    algorithm::split(algorithm::tail(curline), spos, " ");

                    vpos.X = std::stof(spos[0]);
                    vpos.Y = std::stof(spos[1]);
                    vpos.Z = std::stof(spos[2]);

                    Positions.push_back(vpos);
                }
                // Generate a Vertex Texture Coordinate
                // vt : 顶点的纹理坐标（Vector2）(vertex texture)
                if (algorithm::firstToken(curline) == "vt")
                {
                    std::vector<std::string> stex;
                    Vector2 vtex;
                    algorithm::split(algorithm::tail(curline), stex, " ");

                    vtex.X = std::stof(stex[0]);
                    vtex.Y = std::stof(stex[1]);

                    TCoords.push_back(vtex);
                }
                // Generate a Vertex Normal;
                // vn:顶点的法向量（snor:string normal vnor : vector normal）
                if (algorithm::firstToken(curline) == "vn")
                {
                    std::vector<std::string> snor;
                    Vector3 vnor;
                    algorithm::split(algorithm::tail(curline), snor, " ");

                    vnor.X = std::stof(snor[0]);
                    vnor.Y = std::stof(snor[1]);
                    vnor.Z = std::stof(snor[2]);

                    Normals.push_back(vnor);
                }
                // Generate a Face (vertices & indices)
                // f: 面（face） （顶点和索引）
                if (algorithm::firstToken(curline) == "f")
                {
                    // Generate the vertices
                    std::vector<Vertex> vVerts;
                    GenVerticesFromRawOBJ(vVerts, Positions, TCoords, Normals, curline);
                    //(vVert, 顶点位置，纹理坐标，法向量，当前处理行)

                    // Add Vertices
                    // 添加顶点
                    for (int i = 0; i < int(vVerts.size()); i++)
                    {
                        Vertices.push_back(vVerts[i]);

                        LoadedVertices.push_back(vVerts[i]);
                    }

                    std::vector<unsigned int> iIndices;

                    VertexTriangluation(iIndices, vVerts);

                    // Add Indices
                    // 添加索引
                    for (int i = 0; i < int(iIndices.size()); i++)
                    {
                        unsigned int indnum = (unsigned int)((Vertices.size()) - vVerts.size()) + iIndices[i];
                        Indices.push_back(indnum);

                        indnum = (unsigned int)((LoadedVertices.size()) - vVerts.size()) + iIndices[i];
                        LoadedIndices.push_back(indnum);

                    }
                }
                // Get Mesh Material Name
                // 获取网格使用的材质名称 usemtl
                if (algorithm::firstToken(curline) == "usemtl")
                {
                    MeshMatNames.push_back(algorithm::tail(curline));

                    // Create new Mesh, if Material changes within a group
                    // 如果一组材质发生变化，则创建新的网格
                    if (!Indices.empty() && !Vertices.empty())
                    {
                        // Create Mesh
                        tempMesh = Mesh(Vertices, Indices);
                        tempMesh.MeshName = meshname;
                        int i = 2;
                        while(1) {
                            tempMesh.MeshName = meshname + "_" + std::to_string(i);

                            for (auto &m : LoadedMeshes)
                                if (m.MeshName == tempMesh.MeshName)
                                    continue;
                            break;
                        }

                        // Insert Mesh
                        LoadedMeshes.push_back(tempMesh);

                        // Cleanup
                        Vertices.clear();
                        Indices.clear();
                    }

#ifdef OBJL_CONSOLE_OUTPUT
                    outputIndicator = 0;
#endif
                }
                // Load Materials
                // 加载材质
                if (algorithm::firstToken(curline) == "mtllib")
                {
                    // Generate LoadedMaterial

                    // Generate a path to the material file
                    std::vector<std::string> temp;
                    algorithm::split(Path, temp, "/");

                    std::string pathtomat = "";

                    if (temp.size() != 1)
                    {
                        for (int i = 0; i < temp.size() - 1; i++)
                        {
                            pathtomat += temp[i] + "/";
                        }
                    }


                    pathtomat += algorithm::tail(curline);

#ifdef OBJL_CONSOLE_OUTPUT
                    std::cout << std::endl << "- find materials in: " << pathtomat << std::endl;
#endif

                    // Load Materials
                    LoadMaterials(pathtomat);
                }
            }

#ifdef OBJL_CONSOLE_OUTPUT
            std::cout << std::endl;
#endif

            // Deal with last mesh
            // 处理最后一个网格

            if (!Indices.empty() && !Vertices.empty())
            {
                // Create Mesh
                tempMesh = Mesh(Vertices, Indices);
                tempMesh.MeshName = meshname;

                // Insert Mesh
                LoadedMeshes.push_back(tempMesh);
            }

            file.close();

            // Set Materials for each Mesh
            for (int i = 0; i < MeshMatNames.size(); i++)
            {
                std::string matname = MeshMatNames[i];

                // Find corresponding material name in loaded materials
                // when found copy material variables into mesh material
                for (int j = 0; j < LoadedMaterials.size(); j++)
                {
                    if (LoadedMaterials[j].name == matname)
                    {
                        LoadedMeshes[i].MeshMaterial = LoadedMaterials[j];
                        break;
                    }
                }
            }

            if (LoadedMeshes.empty() && LoadedVertices.empty() && LoadedIndices.empty())
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        // Loaded Mesh Objects 加载网格对象
        std::vector<Mesh> LoadedMeshes;
        // Loaded Vertex Objects 加载顶点对象
        std::vector<Vertex> LoadedVertices;
        // Loaded Index Positions 加载索引位置
        std::vector<unsigned int> LoadedIndices;
        // Loaded Material Objects 加载材质对象
        std::vector<Material> LoadedMaterials;

    private:
        // Generate vertices from a list of positions,
        //	tcoords, normals and a face line
        // 从位置列表，纹理坐标，法向量和f生成顶点（相当于将每个顶点和纹理坐标，法向量做对应）
        //(overts ：是输出（out）,iPointions , iTCoords, iNormals, icurline :是传入参数（in）)
        void GenVerticesFromRawOBJ(std::vector<Vertex>& oVerts,
                                   const std::vector<Vector3>& iPositions,
                                   const std::vector<Vector2>& iTCoords,
                                   const std::vector<Vector3>& iNormals,
                                   std::string icurline)
        {
            std::vector<std::string> sface, svert;
            Vertex vVert;
            algorithm::split(algorithm::tail(icurline), sface, " ");

            bool noNormal = false;

            // For every given vertex do this
            // 对每个给出的顶点进行遍历
            // f aa/bb/cc : f 顶点索引/uv点索引/法线索引
            for (int i = 0; i < int(sface.size()); i++)
            {
                // See What type the vertex is.
                int vtype;

                algorithm::split(sface[i], svert, "/");

                // Check for just position - v1 只是顶点索引
                if (svert.size() == 1)
                {
                    // Only position
                    vtype = 1;
                }

                // Check for position & texture - v1/vt1 顶点和uv点索引
                if (svert.size() == 2)
                {
                    // Position & Texture
                    vtype = 2;
                }

                // Check for Position, Texture and Normal - v1/vt1/vn1 顶点索引/uv点索引/法向量索引
                // or if Position and Normal - v1//vn1  顶点索引//法向量索引
                if (svert.size() == 3)
                {
                    if (svert[1] != "")
                    {
                        // Position, Texture, and Normal
                        vtype = 4;
                    }
                    else
                    {
                        // Position & Normal
                        vtype = 3;
                    }
                }

                // Calculate and store the vertex
                switch (vtype)
                {
                    case 1: // P
                    {
                        // 读顶点位置，uv坐标为（0,0）,无法向量
                        vVert.Position = algorithm::getElement(iPositions, svert[0]);
                        vVert.TextureCoordinate = Vector2(0, 0);
                        noNormal = true;
                        oVerts.push_back(vVert);
                        break;
                    }
                    case 2: // P/T
                    {
                        // 读顶点位置，读uv坐标位置，无法向量
                        vVert.Position = algorithm::getElement(iPositions, svert[0]);
                        vVert.TextureCoordinate = algorithm::getElement(iTCoords, svert[1]);
                        noNormal = true;
                        oVerts.push_back(vVert);
                        break;
                    }
                    case 3: // P//N
                    {
                        // 读顶点位置，uv坐标位置置为（0，0）,读法向量位置
                        vVert.Position = algorithm::getElement(iPositions, svert[0]);
                        vVert.TextureCoordinate = Vector2(0, 0);
                        vVert.Normal = algorithm::getElement(iNormals, svert[2]);
                        oVerts.push_back(vVert);
                        break;
                    }
                    case 4: // P/T/N
                    {
                        // 读顶点位置，读uv坐标位置，读法向量位置
                        vVert.Position = algorithm::getElement(iPositions, svert[0]);
                        vVert.TextureCoordinate = algorithm::getElement(iTCoords, svert[1]);
                        vVert.Normal = algorithm::getElement(iNormals, svert[2]);
                        oVerts.push_back(vVert);
                        break;
                    }
                    default:
                    {
                        break;
                    }
                }
            }

            // take care of missing normals
            // these may not be truly acurate but it is the
            // best they get for not compiling a mesh with normals
            // 如果文件中没有给出顶点的法向量索引，通过顶点坐标算出面的法向量
            if (noNormal)
            {
                Vector3 A = oVerts[0].Position - oVerts[1].Position;
                Vector3 B = oVerts[2].Position - oVerts[1].Position;

                Vector3 normal = math::CrossV3(A, B);

                for (int i = 0; i < int(oVerts.size()); i++)
                {
                    oVerts[i].Normal = normal;
                }
            }
        }

        // Triangulate a list of vertices into a face by printing
        //	inducies corresponding with triangles within it
        // 把顶点转化为三角面（oIndices: 索引（输出） iVerts: 顶点及其属性值(输入的参数)）
        void VertexTriangluation(std::vector<unsigned int>& oIndices,
                                 const std::vector<Vertex>& iVerts)
        {
            // If there are 2 or less verts,
            // no triangle can be created,
            // so exit
            // 如果顶点数小于2,顶点不构成面，结束
            if (iVerts.size() < 3)
            {
                return;
            }
            // If it is a triangle no need to calculate it
            // 如果只有一个三角形，则不需要去计算
            if (iVerts.size() == 3)
            {
                oIndices.push_back(0);
                oIndices.push_back(1);
                oIndices.push_back(2);
                return;
            }

            // Create a list of vertices
            std::vector<Vertex> tVerts = iVerts;

            while (true)
            {
                // For every vertex
                // 遍历每一个顶点
                for (int i = 0; i < int(tVerts.size()); i++)
                {
                    // pPrev = the previous vertex in the list
                    // pPrev是列表中的前一个顶点，i为0时就是最后一个（循环）
                    Vertex pPrev;
                    if (i == 0)
                    {
                        pPrev = tVerts[tVerts.size() - 1];
                    }
                    else
                    {
                        pPrev = tVerts[i - 1];
                    }

                    // pCur = the current vertex;
                    // pCur 时当前的顶点
                    Vertex pCur = tVerts[i];

                    // pNext = the next vertex in the list
                    // pNext是当前列表中的下一个顶点
                    Vertex pNext;
                    if (i == tVerts.size() - 1)
                    {
                        pNext = tVerts[0];
                    }
                    else
                    {
                        pNext = tVerts[i + 1];
                    }

                    // Check to see if there are only 3 verts left
                    // if so this is the last triangle
                    // 检查下列表中是否只剩下3个顶点，如果是，这三个顶点组成的三角形就是最后一个三角形
                    if (tVerts.size() == 3)
                    {
                        // Create a triangle from pCur, pPrev, pNext
                        // 将pCur，pPrev, pNext的顶点索引分别push在oIndices中
                        for (int j = 0; j < int(tVerts.size()); j++)
                        {
                            if (iVerts[j].Position == pCur.Position)
                                oIndices.push_back(j);
                            if (iVerts[j].Position == pPrev.Position)
                                oIndices.push_back(j);
                            if (iVerts[j].Position == pNext.Position)
                                oIndices.push_back(j);
                        }

                        tVerts.clear();
                        break;
                    }
                    if (tVerts.size() == 4) // 如果剩下4个顶点
                    {
                        // Create a triangle from pCur, pPrev, pNext
                        // 可以由pCur, pPrev， pNext 构成一个三角形面，将三角形面的三个顶点索引记录下来
                        for (int j = 0; j < int(iVerts.size()); j++)
                        {
                            if (iVerts[j].Position == pCur.Position)
                                oIndices.push_back(j);
                            if (iVerts[j].Position == pPrev.Position)
                                oIndices.push_back(j);
                            if (iVerts[j].Position == pNext.Position)
                                oIndices.push_back(j);
                        }

                        Vector3 tempVec;
                        // 找到pCur, pPrev, pNext之外的剩下一个点
                        for (int j = 0; j < int(tVerts.size()); j++)
                        {
                            if (tVerts[j].Position != pCur.Position
                                && tVerts[j].Position != pPrev.Position
                                && tVerts[j].Position != pNext.Position)
                            {
                                tempVec = tVerts[j].Position;
                                break;
                            }
                        }

                        // Create a triangle from tempVec, pPrev, pNext
                        // 将tempVec, pPrev， pNext构成一个三角形面并记录三角形面的索引
                        for (int j = 0; j < int(iVerts.size()); j++)
                        {
                            if (iVerts[j].Position == pPrev.Position)
                                oIndices.push_back(j);
                            if (iVerts[j].Position == pNext.Position)
                                oIndices.push_back(j);
                            if (iVerts[j].Position == tempVec)
                                oIndices.push_back(j);
                        }

                        tVerts.clear();
                        break;
                    }

                    // If Vertex is not an interior vertex
                    // 当顶点不是内顶点时
                    // 计算两个三维向量之间的夹角
                    float angle = math::AngleBetweenV3(pPrev.Position - pCur.Position, pNext.Position - pCur.Position) * (180 / 3.14159265359);
                    if (angle <= 0 && angle >= 180)
                        continue;

                    // If any vertices are within this triangle
                    // 如果找到其他顶点在构建的三角形内部，当前三个点不构成三角形
                    bool inTri = false;
                    for (int j = 0; j < int(iVerts.size()); j++)
                    {
                        if (algorithm::inTriangle(iVerts[j].Position, pPrev.Position, pCur.Position, pNext.Position)
                            && iVerts[j].Position != pPrev.Position
                            && iVerts[j].Position != pCur.Position
                            && iVerts[j].Position != pNext.Position)
                        {
                            inTri = true;
                            break;
                        }
                    }
                    if (inTri)
                        continue;

                    // Create a triangle from pCur, pPrev, pNext
                    // 根据pCur,pPrev,pNext构建三角形并记录顶点索引
                    for (int j = 0; j < int(iVerts.size()); j++)
                    {
                        if (iVerts[j].Position == pCur.Position)
                            oIndices.push_back(j);
                        if (iVerts[j].Position == pPrev.Position)
                            oIndices.push_back(j);
                        if (iVerts[j].Position == pNext.Position)
                            oIndices.push_back(j);
                    }

                    // Delete pCur from the list
                    // 把列表中的当前顶点删去
                    for (int j = 0; j < int(tVerts.size()); j++)
                    {
                        if (tVerts[j].Position == pCur.Position)
                        {
                            tVerts.erase(tVerts.begin() + j);
                            break;
                        }
                    }

                    // reset i to the start
                    // -1 since loop will add 1 to it
                    // 将i置为-1,保证每次构建一个三角形面之后，循环都是从i=0开始的
                    i = -1;
                }

                // if no triangles were created
                if (oIndices.size() == 0)
                    break;

                // if no more vertices
                if (tVerts.size() == 0)
                    break;
            }
        }

        // Load Materials from .mtl file
        // 从.mtl文件加载材质
        bool LoadMaterials(std::string path)
        {
            // If the file is not a material file return false
            // 如果不是材质文件返回false
            if (path.substr(path.size() - 4, path.size()) != ".mtl")
                return false;

            std::ifstream file(path);

            // If the file is not found return false
            // 找不到文件返回false
            if (!file.is_open())
                return false;

            Material tempMaterial;

            // 记录当前是否读取一个材质信息
            bool listening = false;

            // Go through each line looking for material variables
            // 按行读取
            std::string curline;
            while (std::getline(file, curline))
            {
                // new material and material name
                // 新的材质和材质名称
                if (algorithm::firstToken(curline) == "newmtl")
                {
                    if (!listening)
                    {
                        listening = true;

                        if (curline.size() > 7)
                        {
                            // "newmtl "字符串长度为7,之后的内容为材质名称，用tail方法处理字符串得到
                            tempMaterial.name = algorithm::tail(curline);
                        }
                        else
                        {
                            tempMaterial.name = "none";
                        }
                    }
                    else
                    {
                        // Generate the material

                        // Push Back loaded Material
                        // 将当前加载的材质放在材质加载器中
                        LoadedMaterials.push_back(tempMaterial);

                        // Clear Loaded Material
                        // 清空当前加载的材质（tempMaterial）并加载下一个新材质
                        tempMaterial = Material();

                        if (curline.size() > 7)
                        {
                            tempMaterial.name = algorithm::tail(curline);
                        }
                        else
                        {
                            tempMaterial.name = "none";
                        }
                    }
                }
                // Ambient Color    读材质的环境光系数信息
                if (algorithm::firstToken(curline) == "Ka")
                {
                    std::vector<std::string> temp;
                    algorithm::split(algorithm::tail(curline), temp, " ");

                    if (temp.size() != 3)
                        continue;

                    tempMaterial.Ka.X = std::stof(temp[0]);
                    tempMaterial.Ka.Y = std::stof(temp[1]);
                    tempMaterial.Ka.Z = std::stof(temp[2]);
                }
                // Diffuse Color 读材质漫反射系数信息
                if (algorithm::firstToken(curline) == "Kd")
                {
                    std::vector<std::string> temp;
                    algorithm::split(algorithm::tail(curline), temp, " ");

                    if (temp.size() != 3)
                        continue;

                    tempMaterial.Kd.X = std::stof(temp[0]);
                    tempMaterial.Kd.Y = std::stof(temp[1]);
                    tempMaterial.Kd.Z = std::stof(temp[2]);
                }
                // Specular Color 读材质镜面反射系数信息
                if (algorithm::firstToken(curline) == "Ks")
                {
                    std::vector<std::string> temp;
                    algorithm::split(algorithm::tail(curline), temp, " ");

                    if (temp.size() != 3)
                        continue;

                    tempMaterial.Ks.X = std::stof(temp[0]);
                    tempMaterial.Ks.Y = std::stof(temp[1]);
                    tempMaterial.Ks.Z = std::stof(temp[2]);
                }
                // Specular Exponent 读取材质的镜面反射指数
                if (algorithm::firstToken(curline) == "Ns")
                {
                    tempMaterial.Ns = std::stof(algorithm::tail(curline));
                }
                // Optical Density  光密度信息
                if (algorithm::firstToken(curline) == "Ni")
                {
                    tempMaterial.Ni = std::stof(algorithm::tail(curline));
                }
                // Dissolve 溶解度？信息
                if (algorithm::firstToken(curline) == "d")
                {
                    tempMaterial.d = std::stof(algorithm::tail(curline));
                }
                // Illumination 光照强度
                if (algorithm::firstToken(curline) == "illum")
                {
                    tempMaterial.illum = std::stoi(algorithm::tail(curline));
                }
                // Ambient Texture Map 环境光贴图
                if (algorithm::firstToken(curline) == "map_Ka")
                {
                    tempMaterial.map_Ka = algorithm::tail(curline);
                }
                // Diffuse Texture Map 漫反射贴图
                if (algorithm::firstToken(curline) == "map_Kd")
                {
                    tempMaterial.map_Kd = algorithm::tail(curline);
                }
                // Specular Texture Map 镜面反射贴图
                if (algorithm::firstToken(curline) == "map_Ks")
                {
                    tempMaterial.map_Ks = algorithm::tail(curline);
                }
                // Specular Hightlight Map 高光贴图
                if (algorithm::firstToken(curline) == "map_Ns")
                {
                    tempMaterial.map_Ns = algorithm::tail(curline);
                }
                // Alpha Texture Map alpha贴图
                if (algorithm::firstToken(curline) == "map_d")
                {
                    tempMaterial.map_d = algorithm::tail(curline);
                }
                // Bump Map 凹凸贴图
                if (algorithm::firstToken(curline) == "map_Bump" || algorithm::firstToken(curline) == "map_bump" || algorithm::firstToken(curline) == "bump")
                {
                    tempMaterial.map_bump = algorithm::tail(curline);
                }
            }

            // Deal with last material

            // Push Back loaded Material
            // 处理完最后一个材质信息，将材质信息加入到材质加载器中
            LoadedMaterials.push_back(tempMaterial);

            // Test to see if anything was loaded
            // If not return false
            // 检查是否有成功加载的材质信息，有则返回true，否则返回false
            if (LoadedMaterials.empty())
                return false;
                // If so return true
            else
                return true;
        }
    };
}
#endif //RASTERIZER_OBJ_LOADER_H
