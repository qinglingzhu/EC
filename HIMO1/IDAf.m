%% 函数IDAf判断非劣解的定义
function [DA,time]=IDAf(pa) 
tic; 
[N,C]=size(pa);
DA=ones(N,1); 
for i=1:N 
    temppa=pa; 
    temppa(i,:)=[]; 
    %从种群中筛选支配个体i的所有个体
    for j=1:C  
        LessEqual= temppa(:,j)<=pa(i,j);
        tepa=temppa(LessEqual,:);
        temppa=[]; 
        temppa=tepa;  
    end
    if size(temppa,1)~=0 
       k=1;
       while k<=C
            Lessthan=[];  
            Lessthan=find(temppa(:,k)<pa(i,k));  
            if size(Lessthan,1)~=0 
                DA(i)=0;k=C+1; %有解支配i,则i为支配解（劣解）
            else
                k=k+1;
            end
        end
    end    
end
time=toc;