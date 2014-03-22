%% UDPf函数的定义 旨在选择下一代的种群
function [EPOP,Epa,time]=UDPf(POP,pa,NM)
tic;
[Ns,C]=size(pa);
padis=[];
i=1; 
%去除目标函数值相同的个体
while i<Ns 
    deltf=pa-ones(Ns,1)*pa(i,:);
    deltf(i,:)=inf;  
    aa=find(sum(abs(deltf),2)==0);%aa存储与第i个个体相同的个体
    POP(aa,:)=[];pa(aa,:)=[];
    [Ns,C]=size(pa);%重新计算个体数
    i=i+1;
end
%若除去了相同解的种群数目仍然大于NM（非支配群体个数）则再次根据亲和度来选择下一代种群
if Ns>NM 
    for i=1:C 
        [k,l]=sort(pa(:,i)); 
        N=[];M=[];
        N=pa(l,:);M=POP(l,:);
        pa=[];POP=[];
        pa=N;POP=M;%按照第i个目标函数排序之后的结果
        pa(1,C+1)=Inf;  pa(Ns,C+1)=Inf;
        pai1=[];pad1=[];
        pai1=pa(3:Ns,i);pad1=pa(1:(Ns-2),i);
        fimin=min(pa(:,i));fimax=max(pa(:,i));
        pa(2:(Ns-1),C+1)=pa(2:(Ns-1),C+1)+(pai1-pad1)/(0.0001+fimax-fimin);
    end
    padis=pa(:,C+1); %存储拥挤距离
    pa=pa(:,1:C);
    POP=POP; %此句个人觉得没有很大的意义
    %padis=sum(pa,2);
    [aa,ss]=sort(-padis);%降序排序
    EPOP=POP(ss(1:NM),:);
    Epa=pa(ss(1:NM),:);
else
    EPOP=POP;Epa=pa;  %若对于得到的群体数目不大于NM时 则可以直接全部选中作为下一代
end
time=toc; 