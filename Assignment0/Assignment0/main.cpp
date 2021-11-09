#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>
#include<typeinfo>
using namespace std;
int main() {

    // Basic Example of cpp C++的基本例子
    cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    cout << a << endl;  //输出a
    cout << a / b << endl;    // 输出a/b
    cout << sqrt(b) << endl;   // 输出b的开二次方
    cout << acos(-1) << endl;  // -1 的反正弦值，结果为π
    cout << sin(30.0 / 180.0 * acos(-1)) << endl;   // 30度角的正弦值，转换为弧度制，其中acos(-1)为π， 角度值*π为弧度值

    // Example of vector 向量的例子
    cout << "Example of vector \n";
    // vector definition
    Eigen::Vector3f v(1.0f, 2.0f, 3.0f);  //定义一个三维向量
    Eigen::Vector3f w(1.0f, 0.0f, 0.0f);  //定义一个三维的单位向量
    // vector output
    cout << "Example of output \n";
    cout << v << endl;  //输出三维向量
    // vector add
    cout << "Example of add \n";
    cout << v + w << endl;  // 两向量相加
    // vector scalar multiply
    cout << "Example of scalar multiply \n";
    cout << v * 3.0f << endl;   // 向量的数乘
    cout << 2.0f * v << endl;
    //cout << "cross result:" << typeid(a) << endl;

    // Example of matrix
    cout << "Example of matrix \n";
    // matrix definition
    Eigen::Matrix3f i, j;    //定义3X3的矩阵
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    cout << "Example of output \n";
    cout << i << endl;  // 输出矩阵
    // matrix add i + j 矩阵加法
    cout << i + j << endl;
    // matrix scalar multiply i * 2.0   矩阵数乘
    cout << i * 2.0 << endl;
    // matrix multiply i * j    矩阵乘法
    cout << i * j << endl;
    // matrix multiply vector i * v 矩阵X向量
    cout << i * v << endl;


    // 给定一个点p=(2,1),将该点绕原点先逆时针旋转45度，在平移（1，2），计算出变换后的坐标，要求用齐次坐标进行计算

    Eigen::Vector3f point(2, 1, 1);   // 点
    Eigen::Matrix3f R, T, M;    //旋转和平移矩阵
    float pi = acos(-1);
    float angle = 45.0 / 180.0 * pi;
    R << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
    T << 1, 0, 1, 0, 1, 2, 0, 0, 1;
    cout << R << endl;
    cout << T << endl;
    M = T * R;  // M是两种变换综合起来的齐次方程矩阵
    cout << M << endl;
    cout << "Result of transform:" << endl;
    cout << T * R * point << endl;
    return 0;
}