#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>
#include<typeinfo>
using namespace std;
int main() {

    // Basic Example of cpp C++�Ļ�������
    cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    cout << a << endl;  //���a
    cout << a / b << endl;    // ���a/b
    cout << sqrt(b) << endl;   // ���b�Ŀ����η�
    cout << acos(-1) << endl;  // -1 �ķ�����ֵ�����Ϊ��
    cout << sin(30.0 / 180.0 * acos(-1)) << endl;   // 30�Ƚǵ�����ֵ��ת��Ϊ�����ƣ�����acos(-1)Ϊ�У� �Ƕ�ֵ*��Ϊ����ֵ

    // Example of vector ����������
    cout << "Example of vector \n";
    // vector definition
    Eigen::Vector3f v(1.0f, 2.0f, 3.0f);  //����һ����ά����
    Eigen::Vector3f w(1.0f, 0.0f, 0.0f);  //����һ����ά�ĵ�λ����
    // vector output
    cout << "Example of output \n";
    cout << v << endl;  //�����ά����
    // vector add
    cout << "Example of add \n";
    cout << v + w << endl;  // ���������
    // vector scalar multiply
    cout << "Example of scalar multiply \n";
    cout << v * 3.0f << endl;   // ����������
    cout << 2.0f * v << endl;
    //cout << "cross result:" << typeid(a) << endl;

    // Example of matrix
    cout << "Example of matrix \n";
    // matrix definition
    Eigen::Matrix3f i, j;    //����3X3�ľ���
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    cout << "Example of output \n";
    cout << i << endl;  // �������
    // matrix add i + j ����ӷ�
    cout << i + j << endl;
    // matrix scalar multiply i * 2.0   ��������
    cout << i * 2.0 << endl;
    // matrix multiply i * j    ����˷�
    cout << i * j << endl;
    // matrix multiply vector i * v ����X����
    cout << i * v << endl;


    // ����һ����p=(2,1),���õ���ԭ������ʱ����ת45�ȣ���ƽ�ƣ�1��2����������任������꣬Ҫ�������������м���

    Eigen::Vector3f point(2, 1, 1);   // ��
    Eigen::Matrix3f R, T, M;    //��ת��ƽ�ƾ���
    float pi = acos(-1);
    float angle = 45.0 / 180.0 * pi;
    R << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
    T << 1, 0, 1, 0, 1, 2, 0, 0, 1;
    cout << R << endl;
    cout << T << endl;
    M = T * R;  // M�����ֱ任�ۺ���������η��̾���
    cout << M << endl;
    cout << "Result of transform:" << endl;
    cout << T * R * point << endl;
    return 0;
}