%% 求距离的函数 与前面选择下一代种群比较类似
function [POP,pa,padis]=CDAf(tPOP,tpa) 
[Ns,C]=size(tpa); %将tpa的行号和列号分别赋给Ns和C 而tpa是作为目标函数值
padis=[];  %初始化一个padis的空数组
for i=1:C   %从列号进行循环
    [k,l]=sort(tpa(:,i)); %将第i列从小到大快速排序将排序后的顺序赋给k 并将其所对应原来的行号赋给l
    N=[];M=[]; %初始化两个空数组
    N=tpa(l,:);M=tPOP(l,:);%将tpa按照l的顺序的值存储在N中  tPOP按照l的顺序的值存储在M
    tpa=[];tPOP=[];%将tpa和tPOP赋值空
    tpa=N;tPOP=M;%再将N和M分别赋给tpa和tPOP  
    %以上几个步骤的作用就是为了将tpa按从小到大排序 并将排序后的顺序保存在原来的数组tpa和PtOP中
    tpa(1,C+1)=Inf;tpa(Ns,C+1)=Inf; %将第一行和第NS行 第C+1列的数值赋值为无穷大
    tpai1=[];tpad1=[];%初始化两个数组tpai1和tpad1
    tpai1=tpa(3:Ns,i);tpad1=tpa(1:(Ns-2),i);%将第3行到第NS行的第i列的所有数值赋给tpai1 将第1行到第NS-2行的第i列的所有数值赋给tpad1
    fimin=min(tpa(:,i));fimax=max(tpa(:,i));%找到tpa第i列中的最大和最小 分别赋给max和min
    tpa(2:(Ns-1),C+1)=tpa(2:(Ns-1),C+1)+(tpai1-tpad1)/(0.0001+fimax-fimin);%将tpai1减去tpad1的差并且除以(0.0001+max-min)所得的值
    %与第二行到第NS-1行的第C+1列的数相加 并将值赋给第二行到第NS-1行的第C+1列即 tpa(2:(Ns-1),C+1)
    %所有循环结束后则表示对于每个目标函数得出的每列值中的一个I I的距离为其后一个减去其前面一个
end
pa=tpa(:,1:C);%将tpa的第一列到第C列中所有的数据送给pa
POP=tPOP;%将tPOP中所有的值送给POP
padis=tpa(:,C+1);  %将tpa的第C+!行数据送给padis  因为所有的距离都保存在第C+1列中