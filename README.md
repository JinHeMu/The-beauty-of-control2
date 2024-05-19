# 控制之美2 c++代码

最近在学习控制之美2和视觉SLAM，理论确实非常优美，但是将其实现出来才是画龙点睛之笔。

环境使用SLAM里面的一些库，主要用Eigen库和Pangolin来试图实现其中控制方法，来提升自己的理论到实践的水平，掌握其中库的用法，仅作参考。

## Eigen库

- 数除以向量

  ```cpp
  Eigen::VectorXd T_delta = ((h_max - h_min) / N_h) / v_avg.array()；
  ```

- 向量除以向量

  ```cpp
  Eigen::MatrixXd acc_Matrix = (V_Marix.transpose() - V_Marix).array() / T_delta_Matrix.array();
  ```

- 向量复制

  ```cpp
  Eigen::MatrixXd V_Marix = Vd.replicate(1, 4);
  ```

- 查找最小值和索引

  ```cpp
  minValue = J_temp.minCoeff(&minIndex);
  maxValue = J_temp.maxCoeff(&maxIndex);
  ```

- 按列（行）操作

  ```
  Eigen::VectorXd colMin = J_temp_Martix.colwise().minCoeff();
  ```

  
