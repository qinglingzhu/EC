%% ��ĳ��pareto����Ķ����Բ�����Spacing metric 
% *2014/1/14
% *qlzhu
% *���룺pareto����,���ǽ�ĸ���������Ŀ�꺯��ֵ
% *����������Ա�ʾֵ��ԽС��ʾ����Ŀ��ռ�ķֲ�Խƽ��
function SP=SP(paretof)
[N,C]=size(paretof);
d=zeros(1,N);
for i=1:N
    temp=paretof;
    temp(i,:)=[];
    Tempi=ones(N-1,1)*paretof(i,:);%-1��ȥ������
    d(i)=min(sum(abs(Tempi-temp),2));
end
dave=sum(d)/length(d);
%var����������/N-1��std��׼����/N,�����SP=var(d)^2
SP=sqrt(sum((dave-d).^2,2)/(N-1));
clc;
disp SP;