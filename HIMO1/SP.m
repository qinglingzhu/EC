%% 求某个pareto端面的多样性测量：Spacing metric 
% *2014/1/14
% *qlzhu
% *输入：pareto端面,行是解的个数，列是目标函数值
% *输出：多样性表示值，越小表示解在目标空间的分布越平均
function SP=SP(paretof)
[N,C]=size(paretof);
d=zeros(1,N);
for i=1:N
    temp=paretof;
    temp(i,:)=[];
    Tempi=ones(N-1,1)*paretof(i,:);%-1是去掉自身
    d(i)=min(sum(abs(Tempi-temp),2));
end
dave=sum(d)/length(d);
%var求样本方差/N-1和std标准方差/N,下面的SP=var(d)^2
SP=sqrt(sum((dave-d).^2,2)/(N-1));
clc;
disp SP;