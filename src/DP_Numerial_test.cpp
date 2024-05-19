#include <iostream>
// Eigen 核心部分
#include <eigen3/Eigen/Core>
#include <pangolin/pangolin.h>


int main(int argc, char const *argv[])
{
    //divide
    int N_h = 5;  //height divide
    int N_v = 3;  //speed divide

    //min and max
    int v_min= 0;   //speed min
    int v_max = 3;  //speed max

    float h_min = 0;  //height min
    float h_max = 10; //height max

    //final
    int h_final = 10;
    Eigen::VectorXd v_final(N_v + 1);


    int h_init = 0;
    int v_init = 0;

    float u_min = -0.26; 
    float u_max = 2;

    float h_step = static_cast<float>(h_max - h_min) / N_h;
    // 使用 Eigen::VectorXd 来创建等间距向量
    Eigen::VectorXd Hd(N_h + 1);
    for(int i = 0; i <= N_h; ++i) {
        Hd[i] = h_min + i * h_step;
    }

    float v_step = static_cast<float>(v_max - v_min) / N_v;
    Eigen::VectorXd Vd(N_v + 1);
    for(int i = 0; i <= N_v; ++i) {
        Vd[i] = v_min + i * v_step;
    }

    std::cout << "Hd = " << Hd << std::endl;
    std::cout << "Vd = " << Vd << std::endl;


    Eigen::MatrixXd J_costtogo(N_h + 1, N_v + 1);
    Eigen::MatrixXd Input_acc(N_h + 1, N_v + 1);

    std::cout << "J_costtogo = \n" << J_costtogo << std::endl;
    std::cout << "Input_acc = \n" << Input_acc << std::endl;

    Eigen::VectorXd v_avg = 0.5 * (v_final + Vd);

    std::cout << "v_avg = \n" << v_avg << std::endl;

    Eigen::VectorXd T_delta = ((h_max - h_min) / N_h) / v_avg.array();

    std::cout << "T_delta = \n" << T_delta << std::endl;

    Eigen::VectorXd  acc = (v_final - Vd).array() / T_delta.array();

    std::cout << "acc = \n" << acc << std::endl;
    
    Eigen::VectorXd J_temp = T_delta;


    // 筛选超限的加速度
    Eigen::VectorXd::Index num_over_limit = 0;
    Eigen::VectorXd over_limit_indices(acc.size());
    for (int i = 0; i < acc.size(); ++i) {
        if (acc[i] < u_min || acc[i] > u_max) {
            over_limit_indices[num_over_limit] = i;
            ++num_over_limit;
        }
    }
    std::cout << "over_limit_indices = \n" << over_limit_indices << std::endl;


    for (int i = 0; i < over_limit_indices.size(); ++i) {
        int temp = over_limit_indices[i];
        if (i==0 & over_limit_indices[i]==0) {
            J_temp(0) = std::numeric_limits<double>::infinity();
        }else if(over_limit_indices[i]!=0)
            J_temp(temp) = std::numeric_limits<double>::infinity();
        }
    

    // 更新代价和系统输入矩阵
    J_costtogo.row(1) = J_temp;
    Input_acc.row(1) = acc;

    // 打印结果
    std::cout << "J_costtogo:\n" << J_costtogo << std::endl;
    std::cout << "Input_acc:\n" << Input_acc << std::endl;


//*****************倒数第二级至第二级*************************
    Eigen::MatrixXd V_Marix = Vd.replicate(1, 4);
    
    std::cout << "V_Marix:\n" << V_Marix << std::endl;
    Eigen::MatrixXd V_avg_Marix = 0.5 * (V_Marix + V_Marix.transpose());
    std::cout << "V_avg_Marix:\n" << V_avg_Marix << std::endl;

    Eigen::MatrixXd T_delta_Matrix = ((h_max - h_min) / N_h) / V_avg_Marix.array();

    std::cout << "T_delta_Matrix:\n" << T_delta_Matrix << std::endl;


    Eigen::MatrixXd acc_Matrix = (V_Marix.transpose() - V_Marix).array() / T_delta_Matrix.array();

    std::cout << "acc_Matrix:\n" << acc_Matrix << std::endl;

    Eigen::MatrixXd J_temp_Martix = T_delta_Matrix;



    return 0;
}

