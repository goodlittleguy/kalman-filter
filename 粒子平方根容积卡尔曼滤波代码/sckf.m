function [P,x] = sckf(f,x_old,z,h,P,Q,R)
%% 预设参数
n_x = length(x_old);
n_z = length(z);
m = 2*n_x;
ki = sqrt(m/2)*[eye(n_x),-eye(n_x)];
w = 1/m;
Sq = chol(Q,"lower");
Sr = chol(R,"lower");
S = chol(P);

%% 预设各个变量的空间
sig_pre = zeros(n_x,2*n_x); %存放上一时刻滤波值的预测sigma点
z_sig = zeros(n_z,2*n_x); %存放观测预测sigma点
%% 进行状态量与对应协方差预测
sig = x_old + S*ki; %得到上一时刻滤波值的'sigma'点
for i =1:2*n_x
    sig_pre(:,i) = f(sig(:,i));
end
%sig_pre = sig_pre + sqrtm(Q)*randn(n_x,m);
x_pre = mean(sig_pre,2);  %得到预测值
%Pxx = w*(sig_pre*sig_pre')-x_pre*x_pre'+Q;
[~,Sxx] = qr([sqrt(w)*(sig_pre-x_pre) Sq]',0);
Sxx = Sxx';

%% 更新预测sigma点
sig_pre_new = x_pre + Sxx*ki;

%% 进行观测量预测
for i=1:m
    z_sig(:,i) = h(sig_pre_new(:,i));
end
z_pre = mean(z_sig,2); 
%% 
%{
%% 进行自适应更新
e = z-z_pre;  %新息
e_x = [e(1)*cosd(e(2)),e(1)*sind(e(2))]; 
nmd = [0.95,0.9]*sqrtm(R); %新息各个维度的阈值
c = [0.9,0.9]; %修正系数
token = ones(n_z);
for i=1:n_z
    if abs(e(i)) > nmd(i)
        token(i) = 0;
    end
end

%% 检查是否发生更新,若发生更新则更新对应的预测误差协方差的维数
tmp = 1; %检查是否发生更新
for i=1:n_x
    if token(i)==0
        Sxx(i,i) = Sxx(i,i) + sqrt(R(i,i)*c(i));
        tmp = 0;
    end
end
%如果发生更新也要更新预测观测点
if tmp == 0
    sig_pre_new = x_pre + Sxx*ki;
    for i=1:m
        z_sig(:,i) = h(sig_pre_new(:,i));
    end
    z_pre = mean(z_sig,2); 
end
%}
%% 进行观测量协方差的预测
%Pzz = w*((z_sig-z_pre)*(z_sig-z_pre)')+R;
[~,Szz] = qr([sqrt(w)*(z_sig-z_pre),Sr]',0);
Szz = Szz';
%% 进行状态与观测协方差预测
%Pxz = w*((sig_pre_new-x_pre)*(z_sig-z_pre)');
Pxz = w*((sig_pre_new-x_pre)*(z_sig-z_pre)');


%% 状态与协方差更新
k = (Pxz/Szz')/Szz;
x = x_pre + k*(z-z_pre);
%P = Pxx - k*Pzz*k';
[~,S] = qr([Sxx,k*Sr]',0);
%[~,S] = qr([sqrt(w)*(sig_pre_new-x_pre)-k*sqrt(w)*(z_sig-z_pre),k*Sr]',0);
S = S';
P = S*S';









