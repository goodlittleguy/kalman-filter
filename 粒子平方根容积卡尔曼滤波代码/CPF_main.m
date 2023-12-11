%匀加速直线运动预测模型，观测极坐标模型在无迹卡尔曼滤波中的应用
clc;clear all;
close all;
%% 得到真实数据,第一行距离，第二行方位角，第三行航速，第四行航向
filename = 'D:\WPS Office\sData.xlsx';
[rawData,~,~] = xlsread(filename);
rawData = rawData';
n_std = 4; %共有四个要滤波出来的状态量
N = 100; %共有100组观测数据
std_real = zeros(n_std,N);
for i = 1:n_std
    std_real(i,:) = rawData(i,:);
end

%% 得到状态量X_reals = [x,y,Vx,Vy]
n_x = 4;
x_reals = zeros(n_x,N);
for i=1:N
    x_reals(1,i) = rawData(1,i)*cosd(rawData(2,i)); 
    x_reals(2,i) = rawData(1,i)*sind(rawData(2,i));
    x_reals(3,i) = rawData(3,i)*cosd(rawData(4,i));
    x_reals(4,i) = rawData(3,i)*sind(rawData(4,i));
end

%% 观测量 Z_k = [R_k,she_ta]+v 观测量为距观测点的距离与角度
n_z = 2;
%保存观测值
Z = zeros(n_z,N);
for i=n_std+1:n_std+n_z
    Z(i-n_std,:) = rawData(i,:);
end

z_xy = zeros(2,N);
for i=1:N
    z_xy(1,i) = Z(1,i)*cosd(Z(2,i));
    z_xy(2,i) = Z(1,i)*sind(Z(2,i));
end


%用于计算观测噪声大小
%{
m = zeros(2,N);
for i = 1:N
    m(1,i) = Z(1,i)-std_real(1,i);
    m(2,i) = Z(2,i)-std_real(2,i);
end
%}
%% 得到各个参数
%设每次观测的时间间隔
t = 12.6;
%预测模型求解

%预测噪声量W_k = [W_x,W_y] W_k~N(0,Q)
%预测模型为：X_k+1 = x_k*F+W_k

Q = [0.01 0 0 0;
    0 0.01 0 0;
    0 0 1e-6 0;
    0 0 0  1e-6;];
f = @(x)[x(1)+x(3)*t;x(2)+x(4)*t;x(3);x(4)];

%观测模型求解

%观测噪声 v = [V_r,V_sheta] ～ N(0,R)
R = diag([280^2,4^2]);
n_z = 2 ;%观测模型有两个维度
h = @(x)[sqrt(x(1)^2+x(2)^2);atan2d(x(2),x(1))]; %观测更新函数
%进行ukf变换
Xupf = zeros(n_x,N);
Xupf(:,1) = [z_xy(1,1),z_xy(2,1),0,0]; %0时刻真实状态;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P0 = diag([1000^2,120^2,3^2,1^2]);
n = 100; %一共有100个粒子进行迭代
P = zeros(n_x,n_x,n);
xold = zeros(n_x,n); %存放上一次粒子及更新后的粒子
%生成初始滤波值的粒子
%% 滤波函数
for i=1:n
    P(:,:,i) = P0;
    xold(:,i) = Xupf(:,1) +  sqrtm(P(:,:,i))*randn(n_x,1);
end
xm = mean(xold,2)*ones(1,n);
num = 0;
for i=2:N
    [Xupf(:,i),P,xold,xm] = UPF_fun(xold,xm,P,Z(:,i),R,Q,f,h,n,n_x,n_z);
end
%% 画图
T = 1:N;
figure;
subplot(1,3,1);
plot(x_reals(1,T),x_reals(2,T));
hold on; box on;
plot(Xupf(1,T),Xupf(2,T),'.r');

subplot(1,3,2);
plot(x_reals(1,T),x_reals(2,T));
hold on; box on;
plot(z_xy(1,T),z_xy(2,T),'.r');

subplot(1,3,3);

plot(Xupf(3,T),Xupf(4,T),'.r');

%% 
function [x,P,xold,xm] = UPF_fun(xold,xm,P,z,R,Q,f,h,n,n_x,n_z)

    %存放无迹滤波后的数据预测
    %xold_in = zeros(n_x,n);
    z_pre = zeros(n_z,n);
    %% 进行粒子预测
    for i = 1:n
        %xold_in(:,i) = f(xold(:,i));
        [P(:,:,i),xm(:,i)] = sckf(f,xm(:,i),z,h,P(:,:,i),Q,R);
        xold(:,i) = xm(:,i) + sqrtm(P(:,:,i))*randn(n_x,1);
    end
   
    %% 进行人工蜂群优化算法
    %xold = rbc(xold,z,h,n,2,6);
    
    %% 计算粒子权重
    w = zeros(1,n);
    for i=1:n
        z_pre(:,i) = h(xold(:,i)); %每个粒子的观测
        w1 = inv(sqrt(2*pi*det(R)))*exp(-0.5*(z-z_pre(:,i))'/R*(z-z_pre(:,i)));
        %w2 = exp(-0.5*(xold_in(:,i)-xold(:,i))'/Q*(xold_in(:,i)-xold(:,i)));
        w3 = inv(sqrt(2*pi*det(P(:,:,i))))*exp(-0.5*(xold(:,i)-xm(:,i))'/P(:,:,i)*(xold(:,i)-xm(:,i)));
        w(i) = w1/w3;
    end
    w = w./sum(w); %归一化权值
    outdex = rand_sampling(1:n,w);
    xm = xold(:,outdex);
    P = P(:,:,outdex);
    x = mean(xm,2);

end
