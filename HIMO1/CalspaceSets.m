%%����1��space metric
function [Diversitydeta1,Diversitydeta2]=CalspaceSets(NF,Trials1,Trials2,NS1,NS2, pare1,pare2)
%���룺NF:Ŀ�꺯������
       %Trials1:�㷨1���д���
       %Trials2:�㷨2���д���
       %NS1���㷨1����ִ�еõ��ķ�֧������
       %NS2���㷨2����ִ�еõ��ķ�֧������
       %pare1���㷨1�������д����õ��ķ�֧���Ŀ�꺯��ֵ
       %pare2���㷨2�������еõ��ķ�֧���Ŀ�꺯��ֵ
%�����CT12=Ic(1,2); CT21=Ic(2,1)
Diversitydeta1=[];
Diversitydeta2=[];
j=1;
for time =0:min(Trials1,Trials2)-1%Ϊ��������..�Ҹ���
%for time =1:min(Trials1,Trials2)    ����ԭ����
   pa1=[];pa2=[];
   aa1=sum(NS1(1:((time-1)+1)))+1;bb1=aa1+NS1(time+1)-1;
   aa2=sum(NS2(1:((time-1)+1)))+1;bb2=aa2+NS2(time+1)-1;
   pa1=pare1(aa1:bb1,:);
   pa2=pare2(aa2:bb2,:);
   row=bb1-aa1+1;
   row2=bb2-aa2+1;
   if row>1&&row2>1
   for i=1:row
      temp=pa1;
      temp(i,:)=[];
      Tempi=ones(row-1,1)*pa1(i,:);
      d(i)=min(sum(abs(Tempi-temp),2));
   end
   dave=sum(d)/length(d);
   Diversitydeta1(j)=sqrt(sum((dave-d).^2,2)/(row-1));

   for i=1:row2
      temp=pa2;
      temp(i,:)=[];
      Tempi=ones(row2-1,1)*pa2(i,:);
      d(i)=min(sum(abs(Tempi-temp),2));
   end
   dave=sum(d)/length(d);
   Diversitydeta2(j)=sqrt(sum((dave-d).^2,2)/(row2-1));
   j=j+1;
   end

end
Diversitydeta1=Diversitydeta1';
Diversitydeta2=Diversitydeta2';


% function space()
% Trials1=30;
% Trials2=30;
% NS1
eval(['save ','CalspaceSets',' Diversitydeta1 Diversitydeta2 ']) ;


