%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NPOP,time]=Mutationf(POP,bu,bd,pm,eta_m,ppm,it) %混合变异函数的定义
%--------------------------------------------------------------------------
tic;  %Mutationf函数运行开始的计时器启动
[N,C1]=size(POP); %将解集POP的行号和列号分别赋给N和C1
Fit=POP(:,C1); %将POP中的第C1列赋给Fit 感觉对后面没有起到什么作用
POP=POP(:,1:C1-1); %将POP中的第一列到第C1-1列赋给POP
%eta_m=20;
NPOP=POP; %同时将POP的数据保存在数组NPOP中
C=C1-1; %让C=C1-1
for i=1:N  % 进行循环从第一行到第N行的循环
    r2=rand(); %随机产生一个0到1的随机数
    ppm=(0.3-0.2*it/250)*(Fit(i)-min(Fit))/(max(Fit)-min(Fit)); %相当于论文中的aspi 作为一个阀值来控制两个变异算子的交换
    if r2>ppm %当随机数r2大于这个阀值时 则进行执行多项式变异
    for j=1:C %从POP的第一列到倒数第二列进行循环
        r1=rand; %随机产生一个0到1之间的数送给r1
        if r1<=pm %如果r1小于代码开头给出的变异概率pm
            y=POP(i,j);%将POP（i，j）中的数值送给y
            yd=bd(j);yu=bu(j);%将上界bu和下界bd的第j列赋值给yu和yd
            if  y>yd %如果y大于上界 则进入下面的运算
                if (y-yd)<(yu-y) %若y减去下界yd的值小于上界减去y的值
                    delta=(y-yd)/(yu-yd);%选择小的值除以上界和下界的差
                else
                    delta=(yu-y)/(yu-yd);%选择小的值除以上界和下界的差
                end
                r2=rand; %在随机产生0到1之间的数并将其值赋给r2
                indi=1/(eta_m+1); %eta_m作为变异分布参数
                if r2<=0.5 
                    xy=1-delta;
                    val=2*r2+(1-2*r2)*(xy^(eta_m+1));
                    %val=2*r2;
                    deltaq=val^indi-1;
                else
                    xy=1-delta;
                    val=2*(1-r2)+2*(r2-0.5)*(xy^(eta_m+1));
                    %val=2*(1-r2);
                    deltaq=1-val^indi;
                end
                %针对随机数的值在界点0.5的范围段来分别s是否用的多项式变异公式或者高斯变异 控制变异特征 平均来说eta_m值越大时产生的变异尺度越小
                y=y+deltaq*(yu-yd);%抗体根据变异概率与上下界之差相乘再与原来的抗体相加得到的值就是变异后的抗体
                NPOP(i,j)=min(y,yu); %变异后的抗体与上界进行比较 将小的值取代原来的NPOP（i，j）
                y=NPOP(i,j);
                NPOP(i,j)=max(y,yd);%变异后的抗体与下界进行比较 将大的值取代原来的NPOP（i，j）
            else
                NPOP(i,j)=rand*(yu-yd)+yd;%针对上面if y>yd不成立的情况 要重新产生随机数
            end
        end
    end
 else
    for j=1:C% 针对if r2>ppm不成立时 要进行高斯变异
        r1=rand;%随机产生一个0到1之间的数并将数送给r1
        if r1<=pm %判断若r1小于或等于变异概率pm
            y=POP(i,j); %将POP（i，j）中的数值送给y
            yd=bd(j); yu=bu(j); %将上界bu和下界bd的第j列赋值给yu和yd
            y=y+0.1*normrnd(0,1)*(yu-yd);%normrnd(0,1)为一个平均值为0，标准差为1得瑟高斯随机数 
            NPOP(i,j)=min(y,yu);%变异后的抗体与上界进行比较 将小的值取代原来的NPOP（i，j）
            y=NPOP(i,j);
            NPOP(i,j)=max(y,yd);%变异后的抗体与下界进行比较 将大的值取代原来的NPOP（i，j）
            %NPOP(i,j)=rand*(yu-yd)+yd;
        end
    end
    end
end
time=toc; 