%% ��ʾ��ά����ά��Ŀ�꺯������
function Frontshow(Epa)
%--------------------------------------------------------------------------
if size(Epa,2)==2 %���Epa������Ϊ2  ��Ŀ�꺯���ĸ���������
    %%plot pareto fronts in 2-D
    f1=Epa(:,1);%��һ��Ŀ�꺯����ֵ����f1
    f2=Epa(:,2);%�ڶ���Ŀ�꺯����ֵ����f2
    plot(f1,f2,'r*'); %�Ժ�ɫ����������ʾĿ�꺯��ֵ�ķֲ����
    grid on; %��ʾ����
    xlabel('Function 1');%x������ı��ΪFunction 1
    ylabel('Function 2');%y������ı��ΪFunction 2
    hold on   %��һ��ͼ��ͬʱ��ʾ��������
elseif size(Epa,2)>=3%���Epa����������3  ��Ŀ�꺯��ֵ�ĸ�������������
    %%plot pareto fronts in 3-D
    f1=Epa(:,1);f2=Epa(:,2);f3=Epa(:,3);
    plot3(f1,f2,f3,'kd'); %�Ժ�ɫ��������ʾĿ�꺯��ֵ�ķֲ����
    grid on; %��ʾ����
    xlabel('Function 1');%x������ı��ΪFunction 1
    ylabel('Function 2');%y������ı��ΪFunction 2
    zlabel('Function 3');%z������ı��ΪFunction 3
    hold on     %��һ��ͼ��ͬʱ��ʾ��������
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%