clear all;
N=200;
bsx(1)=1;
p(1)=10;
Z=randn(1,N)+25;
R = std(Z).^2;
w=randn(1,N);
Q = std(w).^2;
for t=2:N;
    x(t)=bsx(t-1);
    
    p1(t)=p(t-1)+Q;
    
    kg(t)=p1(t)/(p1(t)+R);
    
    bsx(t)=x(t)+kg(t)*(Z(t)-x(t));
    
    p(t)=(1-kg(t))*p1(t);
end
t=1:N;
plot(t,bsx,'r', t,Z,'g', t,x,'b');              % 红色线最优化估算结果滤波后的值，%绿色线观测值，蓝色线预测值
legend('Kalman滤波结果','观测值','预测值');
