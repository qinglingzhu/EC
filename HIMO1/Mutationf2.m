%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NPOP,time]=Mutationf2(POP,bu,bd,pm,eta_m)%��ʱΪʲô�ֶ���һ����ϱ���ĺ��� ���������ĺ���û�в���ppm��it �����ʾ�����������Ľ��
%--------------------------------------------------------------------------
tic;  %Mutationf2�������п�ʼ�ļ�ʱ������
[N,C]=size(POP);%���⼯POP���кź��кŷֱ𸳸�N��C
%eta_m=20;
NPOP=POP;%��POP�����ݱ���������NPOP��
for i=1:N %�ӵ�һ�е���N�н����ⲿ��ѭ��
    for j=1:C %�ӵ�һ�е���C�н�����ѭ��
        r1=rand;%�������һ��0��1֮�������� ����ֵ����r1
        if r1<=pm %�ж���r1С�ڻ���ڱ������pm
            y=POP(i,j); %��POP��i��j���е���ֵ�͸�y
            yd=bd(j);yu=bu(j); %���Ͻ�bu���½�bd�ĵ�j�и�ֵ��yu��yd
            NPOP(i,j)=rand*(yu-yd)+yd; %����������彫�丳�� NPOP(i,j)
        end
    end
end
time=toc; 
%--------------------------------------------------------------------------