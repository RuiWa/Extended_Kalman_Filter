%状态量：X = [x, vx, y, vy, z, vz]
%控制量：U = [u1, u2, u3, u4]

close all;
clear all;

%常系数
L= 0.3875;  %单位(m)
Ix = 0.05887;  %单位(kg·m^2)
Iy = 0.05887;
Iz = 0.13151;
g = 9.81; %单位(N/kg)

%动力学方程的常系数
a1 = -(Iy - Iz)/Ix;
a2 = -(Iz - Ix)/Iy;
a3 = -(Ix - Iy)/Iz;  
b1 = L/Ix;
b2 = L/Iy;
b3 = 1/Iz;

Ts = 0.1;                    %采样时间
t = 5;                       %仿真时间
len = fix(t/Ts);            %仿真步数
n = 6;                        %状态维度
w = 0.1;                     %过程标准差
v = 0.5;                      %测量标准差
Q = w^2*eye(n);        %过程方差
R = v^2;                    %测量值的方差

h=@(x)[x(2);x(4);x(6)];                  %测量方程
s=[1;2;3;3;2;1];                            %初始状态
x=s+w*randn(6,1);                      %初始化状态
P = eye(6);                                 %初始化协方差矩阵
xV = zeros(6,len);                       %EKF估计值
sV = zeros(6,len);                       %真实值
zV = zeros(3,len);                       %测量值

for k=1:len
  %随机赋值控制量
  u2 = 0.1*randn(1,1);
  u3 = 0.1*randn(1,1);
  u4 = 0.1*randn(1,1);
  
  z = h(s) + v*randn;                     
  sV(:,k)= s;                             %实际状态
  zV(:,k) = z;                           %状态测量值
  
  %状态方程
  f=@(x)[x(1)+Ts*x(2);
           (a1*x(4)*x(6) +b1*u2)*Ts+x(2);
           x(3)+Ts*x(4);
           (a2*x(2)*x(6) +b2*u3)*Ts+x(4);
           x(5)+Ts*x(6);
           (a3*x(2)*x(4) +b3*u4)*Ts+x(6);];  
  
  %一步预测，同时计算f的雅可比矩阵A
  [x1,A]=jaccsd(f,x); 
  
  %过程方差预测
  P=A*P*A'+Q;         
  
  %状态预测，同时计算h的雅可比矩阵H
  [z1,H]=jaccsd(h,x1); 
  
  %计算卡尔曼增益
  K=P*H'/(H*P*H'+R); 
  
  %状态EKF估计值
  x=x1+K*(z-z1);        
  
  %协方差更新
  P=P-K*H*P;          
  
  xV(:,k) = x;          
  
  %更新状态
  s = f(s) + w*randn(6,1);  
end

%俯仰角、滚转角、偏航角度值
for k=1:2:5
  figure(); hold on; 
  plot(sV(k,:),'-.'); %画出真实值
  plot(xV(k,:)) %画出最优估计值
  plot(abs(sV(k,:)-xV(k,:)), '--'); %画出误差值
  legend('真实状态', 'EKF最优估计估计值', '误差值');
end

%俯仰角速度、滚转角速度、偏航角速度度值
for k=2:2:6
  figure(); hold on; 
  plot(sV(k,:),'-.'); %画出真实值
  plot(xV(k,:)) %画出最优估计值
  plot(zV(k/2,:),':'); %画出状态测量值
  plot(abs(sV(k,:)-xV(k,:)), '--'); %画出误差值
  legend('真实状态', 'EKF最优估计估计值', '状态测量值', '误差值');
end














