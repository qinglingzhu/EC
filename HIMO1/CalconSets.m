%% convergence metric
function [CT1,CT2]=CalconSets(NF,Trials1,Trials2,NS1,NS2, pare1,pare2,Epa)
%���룺NF:Ŀ�꺯������
       %Trials1:�㷨1���д���
       %Trials2:�㷨2���д���
       %NS1���㷨1����ִ�еõ��ķ�֧������
       %NS2���㷨2����ִ�еõ��ķ�֧������
       %pare1���㷨1�������д����õ��ķ�֧���Ŀ�꺯��ֵ
       %pare2���㷨2�������еõ��ķ�֧���Ŀ�꺯��ֵ
%�����CT12=Ic(1,2); CT21=Ic(2,1)
CT1=[];CT2=[];
[row2,col2]=size(Epa);
for time =0:min(Trials1,Trials2)-1%Ϊ��������..�Ҹ���
%for time =1:min(Trials1,Trials2)    ����ԭ����
   pa1=[];pa2=[];
   aa1=sum(NS1(1:((time-1)+1)))+1;bb1=aa1+NS1(time+1)-1;
   aa2=sum(NS2(1:((time-1)+1)))+1;bb2=aa2+NS2(time+1)-1;
   pa1=pare1(aa1:bb1,:);
   pa2=pare2(aa2:bb2,:);
   a=[];
   row=bb1-aa1+1;
       for i=1:row
           Temp=ones(row2,1)*pa1(i,:);
           for j=1:col2
               max1(j)=max(Epa(:,j));
               min1(j)=min(Epa(:,j));
          end
          maxmin=ones(row2,1)*(max1-min1+0.000001);
          % a(i)=min(sqrt(sum(((Temp-Epa)./maxmin).^2,2)));
             a(i)=min(sqrt(sum(((Temp-Epa)).^2,2)));
       end
   CT1(time+1)=mean(a);
   %CT1(time+1)=sqrt(sum(a.^2))/size(a,2);
   a=[];
       row=bb2-aa2+1;
       for i=1:row
           Temp=ones(row2,1)*pa2(i,:);
           for j=1:col2
               max1(j)=max(Epa(:,j));
               min1(j)=min(Epa(:,j));
           end
           maxmin=ones(row2,1)*(max1-min1+0.000001);
         % a(i)=min(sqrt(sum(((Temp-Epa)./maxmin).^2,2)));
             a(i)=min(sqrt(sum(((Temp-Epa)).^2,2)));
       end
   CT2(time+1)=mean(a);
  % CT2(time+1)=sqrt(sum(a.^2))/size(a,2);
end
CT2=CT2';
CT1=CT1';
