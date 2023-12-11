function [Pxx,flag,C0,token] = new_Pxx(Pxx,x_pre_sigma,z_pre_sigma,x_pre,z_pre,z,R,flag,C0,token,Q,w)
%%%%%%%%%%%%%%对Pxx进行收敛性分析 如果收敛下一步，若不收敛，更新Pxx
dec = z - z_pre_sigma;

P_dec = w*(z_pre_sigma-z_pre)*(z_pre_sigma-z_pre)';

S = 590;
thre = S*trace(P_dec);
dec_sum = dec'*dec;
if dec_sum > thre
    disp(6);
    token = 0;
    %标记要传下来
    if flag==1
        C0 = dec*dec';
        flag = 0;
    else
        rou = 0.9;
        C0 = (rou*C0+dec*dec')/(1+rou);
    end
    N = trace((C0-R)');
    M = trace(w*(z_pre_sigma-z_pre)*(z_pre_sigma-z_pre)');
    rou2 = trace(N/trace(M));
    if rou2<1
        rou2 = 1;
    end
    Pxx = rou2*w*(x_pre-x_pre_sigma)*(x_pre-x_pre_sigma)';
    Pxx = Pxx+Q;
end
end
