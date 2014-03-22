%%性能1：coverage of two sets
function [CT21,ave1,fancha1,CT12,ave2,fancha2]=CoverTSets(NF,Trials1,Trials2,NS1,NS2, pare1,pare2)
%输入：NF:目标函数个数
       %Trials1:算法1运行次数
       %Trials2:算法2运行次数
       %NS1：算法1各次执行得到的非支配解个数
       %NS2：算法2各次执行得到的非支配解个数
       %pare1：算法1所有运行次数得到的非支配解目标函数值
       %pare2：算法2所有运行得到的非支配解目标函数值
%输出：CT12=Ic(1,2); CT21=Ic(2,1)
CT12=[];CT21=[]; %初始化两个数组
for time =0:min(Trials1,Trials2)-1%为了能运行..我改了    %对此处的修改有些不是很理解  为什么要做修改
%for time =1:min(Trials1,Trials2)   % 这是原来的
   pa1=[];pa2=[];zz1=[];zz2=[];%初始化四个数组
   sx1=0;sx2=0;  %定义两个变量
   aa1=sum(NS1(1:((time-1)+1)))+1; %NS1(1:((time-1)+1))表示为从第一次到（time-1）+1次运行产生的非支配解的数目 sum则作为将这些都相加 这里要注意的是关于这个time
   %time是一个循环的变量 所以每次循环的aa1是不同的 
   bb1=aa1+NS1(time+1)-1;%当前循环多次所得到的非支配的个数与下一次非支配解的个数相加
   aa2=sum(NS2(1:((time-1)+1)))+1; %NS2(1:((time-1)+1))表示为从第一次到（time-1）+1次运行产生的非支配解的数目 sum则作为将这些都相加 这里要注意的是关于这个time
   %time是一个循环的变量 所以每次循环的aa2是不同的 
   bb2=aa2+NS2(time+1)-1;%当前循环多次所得到的非支配的个数与下一次非支配解的个数相加
   pa1=pare1(aa1:bb1,:);%将第aa1行到bb1行的非支配解针对算法1所得到响应的目标函数解赋给pa1
   pa2=pare2(aa2:bb2,:);%将第aa2行到bb2行的非支配解针对算法1所得到响应的目标函数解赋给pa2
   No1dy2=0;No2dby1=0;%定义两个变量
   for i=1:NS1(time+1)%针对算法一每次产生非支配解的数目进行循环
       j=1;
       aa=1:NS2(time+1);%这句是否为 表示aa的值从1到NS2按顺序每次进行不同的赋值
       AdbyB=1; %将AdbyB的值赋值为1
       while j<=NF % NF是作为目标函数的个数 j作为对列向的一种循环
           zz1=find(pa2( aa,j)<=pa1(i,j)); %从pa2中第aa行第j列查找是否存在小于或等于pa1中第i行第j列的数值 并将其找到的数所对应的行数赋给zz1
           %而aa是表示第time+1次运行所产非支解所有个数中的第一个开始查找
           if length(zz1)>0%若找到符合条件的情况下
               j=j+1; %对下一列的目标函数值进行查找
               aa=aa(zz1);  %我改了这里,请老师看一下我改得对不对  %这里的改变主要目的是什么 没有特别的理解 有待下一步的思考
               %aa=zz1  这是原来的  
           else
               j=NF+1; %作为退出循环的设定 主要是在没有满足length(zz1)>0的情况下就要退出while循环
               AdbyB=0;%对不存在pa2( aa,j)<=pa1(i,j)的一种标记
           end
       end
       if AdbyB==1%这是对存在pa2( aa,j)<=pa1(i,j)的一种设定
           No1dy2=No1dy2+1;%作为表示pa2中的目标函数接支配或等于pa1目标函数解的个数
       end
   end
   for i=1:NS2(time+1)%针对算法二运行每次运行产生非支配解的数目
       j=1;
       aa=1:NS1(time+1);%这句是否为 表示aa的值从1到NS1按顺序每次进行不同的赋值
       BdbyA=1;
       while j<=NF
           zz1=find(pa1(aa,j)<=pa2(i,j)); %从pa1中第aa行第j列查找是否存在小于或等于pa2中第ai行第j列的数值 并将其找到的数所对应的行数赋给zz1
           if length(zz1)>0 %若找到符合条件的情况下
               j=j+1; %对下一列的目标函数值进行查找
               aa=aa(zz1);%我改了这里,请老师看一下我改得对不对
               %aa=zz1  这是原来的
           else
               j=NF+1;  %作为退出循环的设定 主要是在没有满足length(zz1)>0的情况下就要退出while循环
               BdbyA=0; ;%对不存在pa1( aa,j)<=pa2(i,j)的一种标记
           end
       end
       if BdbyA==1%这是对存在pa1( aa,j)<=pa2(i,j)的一种设定
           No2dby1=No2dby1+1;%作为表示pa1中的目标函数接支配或等于pa2目标函数解的个数
       end
   end
   CT12(time+1)=No2dby1/NS2(time+1);%与论文中的Ic（A,B）=(a>=b)/|B|  条件是b属于B 存在A的情况 这里就相当于
   CT21(time+1)=No1dy2/NS1(time+1);
end
ave2=mean(CT12);
fancha2=var(CT12);
ave1=mean(CT21);  %mean表示的是求矩阵或者向量的平均值
fancha1=var(CT21);
CT12=CT12';
CT21=CT21';
eval(['save ','CoverTSetsResult',' CT21 ave1 fancha1 CT12 ave2 fancha2 ']) ;