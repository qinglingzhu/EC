%% UDPf�����Ķ��� ּ��ѡ����һ������Ⱥ
function [EPOP,Epa,time]=UDPf(POP,pa,NM)
tic;
[Ns,C]=size(pa);
padis=[];
i=1; 
%ȥ��Ŀ�꺯��ֵ��ͬ�ĸ���
while i<Ns 
    deltf=pa-ones(Ns,1)*pa(i,:);
    deltf(i,:)=inf;  
    aa=find(sum(abs(deltf),2)==0);%aa�洢���i��������ͬ�ĸ���
    POP(aa,:)=[];pa(aa,:)=[];
    [Ns,C]=size(pa);%���¼��������
    i=i+1;
end
%����ȥ����ͬ�����Ⱥ��Ŀ��Ȼ����NM����֧��Ⱥ����������ٴθ����׺Ͷ���ѡ����һ����Ⱥ
if Ns>NM 
    for i=1:C 
        [k,l]=sort(pa(:,i)); 
        N=[];M=[];
        N=pa(l,:);M=POP(l,:);
        pa=[];POP=[];
        pa=N;POP=M;%���յ�i��Ŀ�꺯������֮��Ľ��
        pa(1,C+1)=Inf;  pa(Ns,C+1)=Inf;
        pai1=[];pad1=[];
        pai1=pa(3:Ns,i);pad1=pa(1:(Ns-2),i);
        fimin=min(pa(:,i));fimax=max(pa(:,i));
        pa(2:(Ns-1),C+1)=pa(2:(Ns-1),C+1)+(pai1-pad1)/(0.0001+fimax-fimin);
    end
    padis=pa(:,C+1); %�洢ӵ������
    pa=pa(:,1:C);
    POP=POP; %�˾���˾���û�кܴ������
    %padis=sum(pa,2);
    [aa,ss]=sort(-padis);%��������
    EPOP=POP(ss(1:NM),:);
    Epa=pa(ss(1:NM),:);
else
    EPOP=POP;Epa=pa;  %�����ڵõ���Ⱥ����Ŀ������NMʱ �����ֱ��ȫ��ѡ����Ϊ��һ��
end
time=toc; 