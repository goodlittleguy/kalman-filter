function outdex = rand_sampling(index,w)
%进行转盘子类型的重新采样
    %首先生成盘子"区间"，如果生成的随机数在i到i+1中间，那么复制第i+1个粒子
    n = length(w);
    outdex = zeros(1,n);
    c = cumsum(w);
    c(n) =1; %一定要注意将盘子补完
    %开始转盘子
    for j=1:n
        a = unifrnd(0,1); %a为（0，1）均匀分布下的随机取值(转到某一个区间)
        %定位该区间
        for k=1:n
            if(a<c(k))
                outdex(j) = index(k);
                break;
            end
        end
    end