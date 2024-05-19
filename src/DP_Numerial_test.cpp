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
    Eigen::VectorXd v_init(N_v + 1);

    float u_min = -3; 
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

    Eigen::MatrixXd J_costtogo(N_h + 1, N_v + 1);

    Eigen::MatrixXd Input_acc(N_h + 1, N_v + 1);

    Eigen::VectorXd v_avg = 0.5 * (v_final + Vd);


    Eigen::VectorXd T_delta = ((h_max - h_min) / N_h) / v_avg.array();


    Eigen::VectorXd  acc = (v_final - Vd).array() / T_delta.array();

    
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


// 从 N_h 开始循环，每次递减 1，直到循环到 1
    for (int k = 2; k < N_h; ++k) {

    Eigen::MatrixXd V_Marix = Vd.replicate(1, 4);
    
    Eigen::MatrixXd V_avg_Marix = 0.5 * (V_Marix + V_Marix.transpose());

    Eigen::MatrixXd T_delta_Matrix = ((h_max - h_min) / N_h) / V_avg_Marix.array();

    Eigen::MatrixXd acc_Matrix = (V_Marix- V_Marix.transpose()).array() / T_delta_Matrix.array();

    Eigen::MatrixXd J_temp_Martix = T_delta_Matrix;

    for (int i = 0; i < acc_Matrix.rows(); ++i) {
        for (int j = 0; j < acc_Matrix.cols(); ++j) {
            if (acc_Matrix(i, j) > u_max) {
                J_temp_Martix(i, j) = std::numeric_limits<double>::infinity();
            }else if (acc_Matrix(i, j) < u_min)
            {
                J_temp_Martix(i, j) = std::numeric_limits<double>::infinity();
            }
            
        }
    }
    std::cout << "J_temp_Martix:\n" << J_temp_Martix << std::endl;

    Eigen::RowVectorXd row = J_costtogo.row(k-1);
    Eigen::MatrixXd J_last_Martix = row.transpose().replicate(1,4);

    std::cout << "J_last_Martix:\n" << J_last_Martix << std::endl;

    J_temp_Martix = J_temp_Martix + J_last_Martix;

    std::cout << "J_temp_Martix:\n" << J_temp_Martix << std::endl;

    Eigen::VectorXd min_row(J_temp_Martix.cols());

    Eigen::VectorXd acc_last_Matrix(J_temp_Martix.cols());

    Eigen::VectorXd colMin = J_temp_Martix.colwise().minCoeff();

    std::cout << "colMin:\n" << colMin << std::endl;
    
    for (int i = 0,count =0; i < J_temp_Martix.rows(); ++i) {
        for (int j = 0; j < J_temp_Martix.cols(); ++j) {
            if(J_temp_Martix(i,j) == colMin(count))
            {
                min_row(count) = i;
                count++;
            }
        }
    }
    
    for (int i = 0; i < acc_last_Matrix.rows(); ++i) {
        acc_last_Matrix(i) = acc_Matrix((int)min_row(i),i);
    }


    std::cout << "Index of minimum value: \n" << min_row << std::endl;

    std::cout << "acc_last_Matrix: \n" << acc_last_Matrix << std::endl;



    // 更新代价和系统输入矩阵
    J_costtogo.row(k) = colMin;
    Input_acc.row(k) = acc_last_Matrix;


    std::cout << "J_costtogo:\n" << J_costtogo << std::endl;
    std::cout << "Input_acc:\n" << Input_acc << std::endl;

    }
    
    v_avg = 0.5 * (v_init + Vd);

    T_delta = ((h_max - h_min) / N_h) / v_avg.array();

    acc = (Vd - v_init).array() / T_delta.array();

    J_temp = T_delta;




    // 筛选超限的加速度
    num_over_limit = 0;
    for (int i = 0; i < acc.size(); ++i) {
        if (acc[i] < u_min || acc[i] > u_max) {
            over_limit_indices[num_over_limit] = i;
            ++num_over_limit;
        }
    }

    for (int i = 0; i < over_limit_indices.size(); ++i) {
        int temp = over_limit_indices[i];
        if (i==0 & over_limit_indices[i]==0) {
            J_temp(0) = std::numeric_limits<double>::infinity();
        }else if(over_limit_indices[i]!=0)
            J_temp(temp) = std::numeric_limits<double>::infinity();
        }

    J_temp = J_temp + J_costtogo.row(N_h-1).transpose();

    std::cout << "J_temp:\n" << J_temp << std::endl;


    Eigen::Index minIndex;
    double minValue;

    minValue = J_temp.minCoeff(&minIndex);

    
    std::cout << "minValue:\n" << minValue << std::endl;
    std::cout << "minIndex:\n" << minIndex << std::endl;


    Eigen::VectorXd J_last(N_v+1);                                         
    Eigen::VectorXd acc_last(N_v+1); 

    J_last.setZero();
    acc_last.setZero();

    J_last(0) = minValue;
    acc_last(0) = acc(minIndex);

    std::cout << "J_last:\n" << J_last << std::endl;
    std::cout << "acc_last:\n" << acc_last << std::endl;

    J_costtogo.row(N_h) = J_last;
    Input_acc.row(N_h) = acc_last;

    // 打印结果
    std::cout << "J_costtogo:\n" << J_costtogo << std::endl;
    std::cout << "Input_acc:\n" << Input_acc << std::endl;

    

    return 0;
}

