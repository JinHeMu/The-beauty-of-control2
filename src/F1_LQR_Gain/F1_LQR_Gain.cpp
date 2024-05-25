//
// Created by jinhemu on 24-5-23.
//
#include <iostream>
#include <limits>
// Eigen 核心部分
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

MatrixXd F1_LQR_Gain(MatrixXd A, MatrixXd B, MatrixXd Q, MatrixXd R, MatrixXd S) {

    int k = 1;
    int max_iter = 200;
    int n = A.rows();    // 计算系统矩阵维度 n
    int p = B.cols();    // 计算输入矩阵维度 p
    double tol = 1e-3;
    double diff = std::numeric_limits<double>::infinity();

    MatrixXd F; // 输出反馈F矩阵
    MatrixXd P0 = S;    //  系统终值代价权重矩阵，定义为P0
    MatrixXd P = MatrixXd::Zero(n, n * max_iter);// 初始化矩阵P为0矩阵，后续用于存放计算得到的一系列矩阵P[k]
    MatrixXd P_k_min_1 = P0;
    MatrixXd F_N_min_k = Eigen::MatrixXd::Constant(P_k_min_1.rows(), P_k_min_1.cols(), std::numeric_limits<double>::infinity());
    MatrixXd F_N_min_k_pre;

    P.block(0, 0, n, n) = P0;

    while(diff > tol) {
        F_N_min_k_pre = F_N_min_k;
        F_N_min_k = (R + B.transpose() * P_k_min_1 * B).inverse() * B.transpose() * P_k_min_1 * A;
        MatrixXd P_k = Q + F_N_min_k.transpose() * R * F_N_min_k + (A - B * F_N_min_k).transpose() * P_k_min_1 * (A-B*F_N_min_k);
        P.middleCols(n * k - n, n) = P_k;
        P_k_min_1 = P_k;
        diff = abs((F_N_min_k - F_N_min_k_pre).maxCoeff());
        k++;
        cout <<  "F_N_min_k\t\t" << F_N_min_k << "\tF_N_min_k_pre\t\t" << F_N_min_k_pre
        << "\tdiff\t\t" << diff << endl;
        if (k > max_iter){
            cerr << "Maximum Number of Iterations Exceeded" << endl;
            break;
        }

    }
    cout << "No. of Interation is " << k << endl;

    F = F_N_min_k;

    return F;
}

int main(int argc, char **argv) {

    Matrix <double,1,1> A, B, Q, R, S;
    A << 1;B << 1;Q << 1;R << 1;S << 1;

    int n = A.rows();    // 计算系统矩阵维度 n
    int p = B.cols();    // 计算输入矩阵维度 p


    // 初始化系统状态
    VectorXd x0 = VectorXd::Constant(n, 1);
    VectorXd x = x0;
    // 初始化系统输入
    VectorXd u0 = VectorXd::Zero(n);
    VectorXd u = u0;

    int k_steps = 20;

    MatrixXd x_history = MatrixXd::Zero(n, k_steps);
    x_history.col(0) = x;

    MatrixXd u_history = MatrixXd::Zero(p, k_steps);
    u_history.col(0) = u;

    int N = k_steps;

    MatrixXd F = F1_LQR_Gain(A, B, Q, R, S);
    for (int i = 1; i < k_steps; ++i) {

        u = -F * x;
        x = A * x + B * u;
        x_history.col(i) = x;
        u_history.col(i) = u;
    }

    cout << "status：" << x_history << endl;
    cout << "intput" << u_history << endl;

    return 0;


}