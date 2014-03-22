%ÐÔÄÜ  Diversity metric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Diversitydeta,meadD,fancha]=Diversity(pare1,NS1,Trials)
%--------------------------------------------------------------------------
for time =0:Trials-1
   pa1=[];
   aa1=sum(NS1(1:((time-1)+1)))+1;
   bb1=aa1+NS1(time+1)-1;
   pa1=pare1(aa1:bb1,:);
   Frontshow(pa1);
   [row,col]=size(pa1);
   if row>1
   for i=1:row
      temp=pa1;
      temp(i,:)=[];
      Tempi=ones(row-1,1)*pa1(i,:);
      d(i)=min(sum(abs(Tempi-temp),2));
   end
   dave=sum(d)/length(d);
   Diversitydeta(time+1)=sqrt(sum((dave-d).^2,2)/(row-1));
   else
       Diversitydeta(time+1)=0;
   end
end
meadD=mean(Diversitydeta);
fancha=var(Diversitydeta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%