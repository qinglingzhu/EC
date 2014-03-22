%% 对重组函数的定义
function [NPOP,time]=Recombinationf(POP,bu,bd,CPOP,it)  
tic;  
eta_c=15;  %相对应于论文中的控制父代个体交叉程度的数量级
[N,C1]=size(POP);
Fit=POP(:,C1); %将解集数组POP的第C1列赋值给Fit中
POP=POP(:,1:C1-1); %将解集POP中第一列到第C1-1列中的数据赋给POP 即将原来的POP去除最后一列的情况
C=C1-1; %用C来表示C1-1
NPOP=POP; %将POP的值同时保存在数组NPOP
for i=1:N  %进行对POP的每行进行循环
    r1=rand; %随机产生一个0到1之间的数送给r1
    if r1<=1.1 % 这一步是一个肯定成立的过程 为什么要这样设计 有点不明白
        aa=randperm(size(CPOP,1));%randperm是matlab中自带的函数 randperm（x）表示随机生成从1到X的所有整数
        bb=aa(1);%将aa数组的第一个数赋给bb
        for j=1:C %对POP的列号开始进行循环
            par1=POP(i,j);%将第i行第j列的数据赋给par1
            par2=CPOP(bb,j);%将第bb行第j列的值赋给par2
            yd=bd(j);yu=bu(j);%将上界bu和下界bd的第j列赋值给yu和yd
            r2=rand; %随机产生0到1的随机数并将值赋给r2
            if r2<=0.5 %对r2进行判断 若r2小于0.5进行下面步骤
                if abs(par1-par2)>10^(-14)%若par1和par2的差的绝对值大于10的负14次方
                    y1=min(par1,par2);y2=max(par1,par2); %找到par1和par2中小的赋给y1 大的值赋给y2
                    if (y1-yd)>(yu-y2)%如果y1减去下界yd的值大于上界yu减去y2
                        beta=1+2*(yu-y2)/(y2-y1); %则分子取yu减去y2
                    else
                        beta=1+2*(y1-yd)/(y2-y1);%否则分子y1减去yd
                    end
                    expp=eta_c+1;beta=1/beta;alpha=2.0-beta^(expp); %与论文的公式相对应 由此得出beta的值小于1 则相对应的beta^(expp)小于1
                    % 而得到alpha的范围是在1到2之间  怎1/alpha的范围是0.5到1之间
                    r3=rand; %再次随机产生一个0到1的数值赋给r3
                    if r3<=1/alpha  %跟产生一个平均分布的随机数有什么不同呢?  这里alpha的值要受到不同的上下界的变化而变化 但其所有的范围都在0.5到1之间
                        alpha=alpha*r3;expp=1/(eta_c+1.0);
                        betaq=alpha^(expp);  %即论文中的控制交叉参数
                    else
                        alpha=1/(2.0-alpha*r3);expp=1/(eta_c+1);
                        betaq=alpha^(expp);%在产生的随机数大于1/alpha时 则betaq等于alpha^(expp)上面几步以及参数的设置 可以得到结论是当参数eta_c的值越大时
                        %重组生成的子代就以较高的概率更加接近父代 因为根据所给出的公式 我们可以判断出当eta_c越大是
                        %则相应的alpha就越大 因此其betaq就越大 但是betaq的值是小于1的 eta_c越大则betaq就更加向1接近
                    end  %父代抗体之间相近则交叉生成的子代越相似,
                    chld1=0.5*((y1+y2)-betaq*(y2-y1));%可以化简做chid1=0.5*（（1+betaq）*y1+（1-betaq）*y2）
                    chld2=0.5*((y1+y2)+betaq*(y2-y1)); %可以化简做chid1=0.5*（（1+betaq）*y2+（1-betaq）*y1）
                    %上面几步以及参数的设置 可以得到结论是当参数eta_c的值越大时重组生成的子代就以较高的概率更加接近父代 因为根据所给出的公式 我们可以判断出当eta_c越大是
                        %则相应的alpha就越大因此其betaq就越大当betaq越大时则根据上面的化简可以明确chil1的值趋向于y1 chil2的值趋向于y2
                    if rand<=0.5 %随机产生一个0到1范围的饿随机数 看其是否小于或等于0.5 
                        aa=max(chld1,yd);%将重组产生的子代chil1与下界进行比较 将大的值送给aa
                        NPOP(i,j)=min(aa,yu);%将aa与上界进行比较 将小的那个值赋给NPOP的第i行第j列
                    else
                        aa=max(chld2,yd);%将重组产生的子代chil2与下界进行比较 将大的值送给aa
                        NPOP(i,j)=min(aa,yu);%将aa与上界进行比较 将小的那个值赋给NPOP的第i行第j列
                    end
                    %这几步的步骤主要是重组产生的子代根据随机的概率 看将产生的两个子代中那个子代来取代原来解集中的NPOP（i，j）
                end  
            end
        end
    end
end
NPOP=[NPOP Fit];
time=toc; 