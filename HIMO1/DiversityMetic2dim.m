%% Diversity metric delta
function DiversityMetic2dim(Trials,NS,PFknow,PFtrue)
if nargin<2
    error('输入参数错误');
end
%在PFture中找极端点
[Nk,Ck]=size(PFknow);[Nt,Ct]=size(PFtrue);
if(Ct~=2)
    error('Diversity metric delta只能计算2个目标函数');
end
[c d]=sort(PFtrue(:,1));pt=PFtrue(d,:);clear c d;
for time=0:Trials-1
    %取出该次实验存储的paretof的范围aa->bb
    aa=sum(NS(1:((time-1)+1)))+1;bb=aa+NS(time+1)-1;
    %获得该次实验的paretof
    pa=PFknow(aa:bb,:);
    %将PFknow和PFtrue按某个目标函数排序
    [a b]=sort(pa(:,1));pk=pa(b,:);clear a b;
    %先算df
    df=[];dl=[];d=[];
    df=sqrt(sum(abs(pk(1,:)-pt(1,:)).^2,2));
    index=1;
    for i=2:NS(time+1)
        d(index)=sqrt(sum(abs(pk(i,:)-pk(i-1,:)).^2,2));index=index+1;
    end
    %最后选dl
    dl=sqrt(sum(abs(pk(end,:)-pt(end,:)).^2,2));
    dave=sum(d)/length(d);
    delta(time+1)=(df+dl+sum(abs(d-dave)))/(df+dl+(NS(time+1)-1)*dave);
end   
disp(mean(delta));

% x=1:1:1000;r=rand;
% disp(r);
% hold on;plot(1-r.^((1-x.*0.001).^2));
% index=1;
% for i=1:2:size(ZDT1,1)
%     PFtrue(index,:)=ZDT1(i,:);index=index+1;
% end