%% 克隆函数的定义
function  [NPOP,time]=Clonef(POP,pa,CS)  
tic   
NC=[];
[N,C]=size(POP); 
[POP,pa,padis]=CDAf(POP,pa);%调用下面的CDAf函数算出种群之间距离,根据距离进行克隆
aa=find(padis==inf); %查找padis等于无穷大的情况 即padis的第一行和第Ns行
bb=find(padis~=inf); %查找padis不等于无穷大的情况 
if length(bb)>0  %若存在padis不等于无穷大的情况下
    padis(aa)=2*max(max(padis(bb))); %把边界点的抗体亲合度重新设定为除边界点外最大亲合度的两倍 即边界点的距离设置为除边界点外最大距离的两倍
    NC=ceil(CS*padis./sum(padis));%ceil函数是matlab自带的函数 ceil（X）则表示大于X的最小整数 CS为预定义克隆总数 根据亲合度（即拥挤距离）
    %来定义每个解的克隆个数 与论文中相对应 亲合度与克隆个数成正比 亲合度越大 即（距离越大）则克隆的个数越多 以此很好的保持种群的多样性
else
    NC=ceil(CS/length(aa))+zeros(1,N); %如果不存在bb长度大于零时 则克隆的个数为总的克隆数目除以aa的长度 即平均分配每个抗体的克隆数目
end
NPOP=[]; %将初始化定义一个新数组NPOP
FitPop=[POP padis]; %将解集和距离结合赋给数组FitPop
for i=1:N %将数组循环
    NiPOP=ones(NC(i),1)*FitPop(i,:); %根据NC得到每个抗体的的克隆数目
    NPOP=[NPOP;NiPOP]; %将之前的解和克隆后的解按行组合在一起从NPOP的最后一行开始添加NiPOP 
end
time=toc; 