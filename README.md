## 基于EKF的四旋翼无人机姿态估计

> 说明：题为作者本科毕业设计的选题，使用EKF(Extended Kalman Filter, 扩展卡尔曼滤波)算法来对四旋翼无人机的姿态进行滤波和估计，姿态包括：俯仰角、滚转角、偏航角的角度值和角速度值。前提：角度值无法直接通过传感器直接测得，角速度值可以测得。

### 代码说明

- test1.m：一维线性卡尔曼滤波的例子
- jaccsd.m：用于求解EKF算法中的雅克比矩阵
- EKF.m：EKF算法仿真程序

### 仿真结果

- 说明：
- 1.仿真软件采用MATLAB2010b
- 2.控制量和姿态角速度值采用随机生成的数据(使用实际数据更好)
- 3.仿真过程偶尔会出现错误结果，原因是EKF计算过程中有几率出现奇异矩阵，导致算法无法进行下去
