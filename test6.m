% 
N =150;         %计算连续N个时刻 
n=3;            %状态维度
q=0.1;          %过程标准差
r=0.2;          %测量标准差
Q=q^2*eye(n);   %过程方差
R=r^2;          %测量值的方差 
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  %状态方程
h=@(x)[x(1);x(2);x(3)];                   %测量方程
s=[0;0;1];                                %初始状态
%初始化状态
x=s+q*randn(3,1);                         
P = eye(n);                               
xV = zeros(n,N);          
sV = zeros(n,N);         
zV = zeros(n,N);
for k=1:N
  z = h(s) + r*randn;                     
  sV(:,k)= s;                             %实际状态
  zV(:,k)  = z;                           %状态测量值
  [x1,A]=jaccsd(f,x); %计算f的雅可比矩阵，其中x1对应黄金公式line2
  P=A*P*A'+Q;         %过程方差预测，对应line3
  [z1,H]=jaccsd(h,x1); %计算h的雅可比矩阵
  K=P*H'*inv(H*P*H'+R); %卡尔曼增益，对应line4
  x=x1+K*(z-z1);        %状态EKF估计值，对应line5
  P=P-K*H*P;            %EKF方差，对应line6
  xV(:,k) = x;          %save
  s = f(s) + q*randn(3,1);  %update process 
end

for k=1:3
  FontSize=14;
  LineWidth=1;
  figure();
  plot(sV(k,:),'g-'); %画出真实值
  hold on;
  plot(xV(k,:),'b-','LineWidth',LineWidth) %画出最优估计值
  hold on;
  plot(zV(k,:),'k'); %画出状态测量值
  hold on;
  legend('真实状态', 'EKF最优估计估计值','状态测量值');
  xl=xlabel('时间(分钟)');
  t=['状态 ',num2str(k)] ;
  yl=ylabel(t);
  set(xl,'fontsize',FontSize);
  set(yl,'fontsize',FontSize);
  hold off;
  set(gca,'FontSize',FontSize);
end

figure();
plot(abs(sV(3,:)-xV(3,:)));