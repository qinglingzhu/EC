%% PFknow与PFtrue的收敛程度
function Convergence_metic(Trials,NS,PFknow,PFtrue)
[N1,C]=size(PFknow);[N2,C]=size(PFtrue);clear C;
CT=[];
for time=0:Trials-1
    pa1=[];pa2=[];
    %NS存储每次运行所得到的最优端面的个体数，则aa1->aa2是第time次运行结果区间
    aa1=sum(NS(1:((time-1)+1)))+1;bb1=aa1+NS(time+1)-1;
    %第time次运行得到的pareto端面
    pa1=PFknow(aa1:bb1,:);
    a=[];
    for i=1:NS(time+1)
        Temp=ones(N2,1)*PFknow(i,:);
        a(i)=min(sqrt(sum(((Temp-PFtrue)).^2,2)));
    end
    CT(time+1)=mean(a);
end
disp(mean(CT));
% index=1;
% for i=1:2:1000
%     PFtrue(index,:)=ZDT1(i,:);index=index+1;
% end