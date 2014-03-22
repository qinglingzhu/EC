function HIMO()
EMOinstruction;
%--------------------------------------------------------------------------
TestNO=input('press the enter key after inputting the serial number of test problem:');
Trial=input('input the number of independent runs:');
p=0.2;
NA=20;                                                      % 活动群体的个数 
CS=100;                                                     % 克隆总数  
NM=100;                                                     % 支配群体的个数(把较优的解放在这个群体里面)
gmax=250;                                                   % 最大进化代数 
[bu,bd,testfunction]=getbud(TestNO);                        % bu 代表变量的上界;bd 代表变量的下界 对另外一个函数getbud.m的调用 TestNO是作为用户输入的选择值
c=size(bu,2);                                               % 变量的个数
pmy=1/c;                                                    % 设置变量的变异概率 
p=0.2;                                                      % 预定义参数 给定的变异相关的参数
if bu==bd 
    return;
end                                        % 当上界和下界相等时停止算法  
%--------------------------------------------------------------------------
Datime=date;        %获取年月日                                      
Datime(size(Datime,2)-4:size(Datime,2))=[];%测试时间 与上一行结合得到月与日
TestTime=clock;%测试时间
TestTime=[num2str(TestTime(4)),'-',num2str(TestTime(5))];%与上一行结合得到时间的显示 即时与分
Method='HIMO';
paretof=[]; %存储求得的非劣解 作为存储机制
runtime=[];%首先将运行的次数初始化为零
for trial=1:Trial %作为一个大循环 即要进行多少次的运行与最上面的相对应Trial=input('input the number of independent runs:');
                   %Trial是HIMO主函数运行时用户输入的值 即运行时time的计数 即 进行trial次250次运行
    timerbegin=clock; %开始计时
%--------------------------------------------------------------------------
POP=rand(NM,c).*(ones(NM,1)*bu-ones(NM,1)*bd)+ones(NM,1)*bd;  

load Kursawe.pf.txt
figure;Frontshow(Kursawe_pf,'k.');axis([-20 -4 -15 25]);%plot the PFture

ME=[];%Initialization对数组ME初始化

%--------------------------------------------------------------------------
pa=OVcom(POP,TestNO); %返回多目标函数值 对函数OVcom.m的调用 即将随机生成的种群进行代入到测试函数中所返回的值赋给数组pa
hold on;Frontshow(pa,'b*');%plot particles

[DON,DONt]=IDAf(pa);  %判断是否为非劣解,0是较差解,1是非劣解 对后面IDAf函数的调用即将上面所得到的目标函数值pa进行判断 
                      %判断的返回值赋给DON 0表示劣解  1表示非劣解 

nodom=find(DON==1);  %求得所有非劣解  find这个函数是表示求满足条件DON=1的情况下的所在的行号赋给nodom
MEpa=pa(nodom,:);MEPOP=POP(nodom,:);%获得所有非劣解的各个参数值和目标函数值 
                                    %根据上面得到的行号可以将之前的参数值与目标函数值进行劣解和非裂解的划分 以及劣值与非劣质的划分
Frontshow(MEpa,'r*');%plot nd-particles

numnod=size(MEPOP,1); %非劣解个数 size(MEPOP,1)函数是指返回MEPOP的行号
[ClonePOP,Clonepa,RNDCDt]=UDPf(MEPOP,MEpa,NA);%按距离拥挤程度排列,组成活动子群 对后面UDPf函数的调用 
Frontshow(Clonepa,'co');%plot nd-particles
legend('PFture','POP','非支配解','被选中克隆');
 %即在判断完非劣解的情况下所得到的解的个数可能仍然大于NM(支配群体的总个数) 这是后要进行选择下一代进行克隆的数目
it=0;Cloneti=0;DASti=0;PNmti=0;DONjudti=0;RNDCDti=0;AbAbfitcomti=0;
%it 进化代数 Cloneti:克隆所花费的时间 DASti:交叉所花费的时间 
%PNmti:变异所花费的时间 DONjudti:判断支配关系所花费的时间 RNDCDti:选择所花费的时间
%AbAbfitcomti:后面没再出现,无用数据
while it<gmax %it作为运行代数  gmax作为最大的运行代数 作为限制循环次数
%--------------------------------------------------------------------------
eta_m=20; %类似于论文中在GP-HM模块中所设置的变异分布指数n
pm=-1*p*pmy*(it/125)+(1+p)*pmy;   %相当于论文中的pm=（1+p）*pmy-2*p*pmy*（it/250）  p是刚开始给出的变量的变异概率
%但是此句所有的问题是论文中有一条限制条件it<gmax/2   设置变量的变异概率
if pm<pmy
    pm=pmy;  %论文中的pmy所表示的是预定义最小变异概率
end
ppm=0.2; %这个变量是否应该是pmy 否则上面的pm就无法运算  % 似乎在后面的混合变异的运算中参数中出现了ppm 否定当初的想法
cloneover=[]; %对克隆种群的数组进行初始化
[cloneover,Clonet]=Clonef(ClonePOP,Clonepa,CS);  %对活动子群进化克隆,返回克隆子群和操作时间
[cloneover,DASt]=Recombinationf(cloneover,bu,bd,ClonePOP,it);%对克隆子群进化模拟二进制交叉,返回交叉后克隆子群和操作时间
[cloneover,PNmt]=Mutationf(cloneover,bu,bd,pm,eta_m,ppm,it); %对交叉后子群进行多项式变异,返回克隆子群和操作时间
%--------------------------------------------------------------------------
clonepa=OVcom(cloneover,TestNO); %对函数OVcom的调用获取变异后克隆子群的目标函数值

clf;Frontshow(Kursawe_pf,'k.');axis([-20 -4 -15 25]);%plot the PFture
hold on;Frontshow(clonepa,'b*');%plot particles
NPOP=[MEPOP;cloneover]; %变异后克隆子群与原来非劣解组成下一代种群
Npa=[MEpa;clonepa];%变异后克隆子群的目标函数值与原来非劣解的目标函数值组后在一起送给Npa

[NDON,DONjudt]=IDAf(Npa);%对函数IDAf的调用将结合了的目标函数值Npa进行判断非劣解
Nnodom=find(NDON==1);%将上一步找到的非劣解进行查到非劣解的行号并返回赋给Nnodom
NEpa=Npa(Nnodom,:);NEPOP=NPOP(Nnodom,:);%复制非劣群体  将非劣解和非劣解所对应的目标函数值分别赋给两个数组NEPOP和NEpa

Nnumnod=size(NEPOP,1);%获得非劣群体的个数 size(NEPOP,1)函数matlab自带返回NEPOP的列号
[MEPOP MEpa,RNDCDt]=UDPf(NEPOP,NEpa,NM);%对函数UDPf的调用进行选择下一代种群 %即在非劣解大于NM（支配群体的总个数）时要进行的进一步的选择
%[MEPOP MEpa]=clustering(NEPOP,NEpa);    将选择的下一代种群的非劣解与所对应的目标函数值分别赋给MEPOP和MEpa
Frontshow(MEpa,'r*');%plot nd-particles
numnod=size(MEPOP,1);%选择的非劣群体的个数 size(MEPOP,1）函数matlab自带 返回MEPOP的列号
[ClonePOP,Clonepa,RNDCDt]=UDPf(MEPOP,MEpa,NA);%选择克隆子群 NA作为活动种群的数目即选择NA以内的的种群数目进行克隆
Frontshow(Clonepa,'co');%plot nd-particles
legend('PFture','POP','非支配解','被选中克隆');
%Update Dominant Population
%--------------------------------------------------------------------------
it=it+1; %与上面的while循环相对应 在it<gmax时循环一直进行下去
Cloneti=Cloneti+Clonet;%将每一代克隆所花费的时间相加
DASti=DASti+DASt; %将每一代交叉所花费的时间相加
PNmti=PNmti+PNmt;%将每一代变异所花费的时间相加
DONjudti=DONjudti+DONjudt;%将每一代进行判断非劣解（即判断支配关系）的时间相加 
RNDCDti=RNDCDti+RNDCDt;%把所有每一代选择所花费的时间累加起来
disp(sprintf('time: %d   generation: %d    number of nodominate:  %d',trial,it,numnod));%time即之前用户输入的trial 
%generation是作为运行的代数 本论文设置为250代 nodominate为非支配解的数目 本论文中设置为NM（100）
end  %the end of iterations 结束循环
%--------------------------------------------------------------------------
%Save the output solutions
[NS(trial),NF]=size(MEpa);%size(MEpa)为求得MEpa的行号与列号分别赋给NS（trial）和NF
Trials=trial; 
paretof=[paretof;MEpa]; %paretof作为存储非劣解
runtime(trial)=etime(clock,timerbegin); % etime(t1,t2)用来计算两个日期向量t1和t2之间的时间差 即在此为计算clock和timerbegin之间的时间差 
%此处的时间为整个程序的运行时间
Clonetime(trial)=Cloneti;%将所有代克隆花费总的时间赋给Clonetime(trial)
DAStime(trial)=DASti;%将所有代交叉所花费总的时间赋给DAStime(trial)
PNmtime(trial)=PNmti; %将所有代变异所花费总的时间赋给PNmtime(trial)=PNmti;
DONjudtime(trial)=DONjudti;%将所有代判断非劣解（即判断支配关系）所花费的总的时间赋给DONjudtime(trial)
RNDCDtime(trial)=RNDCDti;%将所有代选择所花费的总的时间赋给RNDCDtime(trial)
eval(['save ', testfunction Method Datime TestTime ,' Method gmax NA CS paretof runtime Trials NS NF TestTime Datime testfunction ']) ;
% Frontshow(MEpa);
% pause(1);
% hold off;
%eval()matlab中自带的库函数 函数的功能就是将括号内的字符串视为语句并运行 相当于在命令行输入命令
end  %the end of runs  
%--------------------------------------------------------------------------
Frontshow(MEpa);% plot the Pareto fronts solved by the last run 对后面输出函数的调用 输出目标函数值的图形分布

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DA,time]=IDAf(pa) %函数IDAf判断非劣解的定义
%--------------------------------------------------------------------------
tic; %程序开始的计时器启动
[N,C]=size(pa); %N:种群个体数量，C目标函数空间维数
DA=ones(N,1); %每个个体初始为非支配解
for i=1:N 
    temppa=pa; 
    temppa(i,:)=[]; %去除本身？？
    LEsign=ones(N-1,1);%这句话在后面没有用到 个人感觉没有什么作用
    for j=1:C  %取出第一层非支配解->temppa
        LessEqual=find(temppa(:,j)<=pa(i,j));
        tepa=[];%初始化数组tepa
        tepa=temppa(LessEqual,:);%将上面满足条件数据行号赋给数组tepa
        temppa=[];  %将temppa重新初始化
        temppa=tepa;  %将数组tepa值赋给temppa
    end
    if size(temppa,1)~=0 
        k=1;  %定义一个变量k
        while k<=C   %让k从第一列开始循环到第C列
            Lessthan=[];  %初始化数组lessthan
            Lessthan=find(temppa(:,k)<pa(i,k));  %在满足上面条件的基础上（即temppa除去第i行数据的第j列所有数据小于或等于pa中的第一行第一列的数)加以限制使只小于的情况
            if size(Lessthan,1)~=0 %若找到只小于的情况的值 则返回的下标值不等于1
                DA(i)=0;k=C+1;
            else
                k=k+1;%若没有找到小于的情况 只满足于等于的情况 则将列号加1 进行下一列的比较
            end
        end
    end    
end
time=toc;%与函数开始设置的tic相对应 tic用在程序的开始，作用是启动一个计时器，然后在程序尾部放一个toc，表示终止计时器，并返回tic启动以来的总时s间

%%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EPOP,Epa,time]=UDPf(POP,pa,NM) %UDPf函数的定义 旨在选择下一代的种群
%--------------------------------------------------------------------------
tic; %程序开始的计时器启动
[Ns,C]=size(pa); %目标函数值pa的行号与列号分别赋值给Ns和C
padis=[]; %初始化一个空数组 主要存储拥挤距离 即可以表示论文中的抗体的亲合度
i=1; 
while i<Ns   %删除相同的解 
    deltf=pa-ones(Ns,1)*pa(i,:);%表示在不同的pop 可能会得到相同的目标函数值pa 则这里所做的就是删除pa相同行保留一行即可
    deltf(i,:)=inf;  %将本身的第i行赋值为无穷大
    aa=find(sum(abs(deltf),2)==0); %sum(x,2)为matlab的自带库函数 表示将x每行的值相加 这里是指将deltf每行数的绝对值相加查找等于0的情况 
   % 并返回其行号   也就是找到了pa中存在相同的值的行号并且赋给aa
    POP(aa,:)=[];pa(aa,:)=[];  %将第aa行所对应的解POP和目标函数值赋值为空
    [Ns,C]=size(pa);%将取出了相同值的pa进行size求出行号的列号
    i=i+1; %从第一行到NS行进行循环
end
if Ns>NM %若除去了相同解的种群数目仍然大于NM（非支配群体个数）则再次根据亲和度来选择下一代种群
    for i=1:C  %从列号开始进行循环
        [k,l]=sort(pa(:,i)); %将第i列从小到大快速排序将排序后的顺序赋给k 并将其所对应原来的行号赋给l
        N=[];M=[];%初始化两个空数组
        N=pa(l,:);M=POP(l,:);%将pa按照l的顺序的值存储在N中 POP按照l的顺序的值存储在M
        pa=[];POP=[];%将pa和POP赋值空
        pa=N;POP=M;%再将N和M分别赋给pa和POP  
        %以上几个步骤的作用就是为了将pa按从小到大排序 并将排序后的顺序保存在原来的数组pa和POP中
        pa(1,C+1)=Inf;  pa(Ns,C+1)=Inf;%将第一行和第NS行 第C+1列的数值赋值为无穷大
        pai1=[];pad1=[];%初始化两个数组pai1和pad1
        pai1=pa(3:Ns,i);pad1=pa(1:(Ns-2),i);%将第3行到第NS行的第i列的所有数值赋给pai1 将第1行到第NS-2行的第i列的所有数值赋给pad1
        fimin=min(pa(:,i));fimax=max(pa(:,i));%找到pa第i列中的最大和最小 分别赋给fimax和fimin
        pa(2:(Ns-1),C+1)=pa(2:(Ns-1),C+1)+(pai1-pad1)/(0.0001+fimax-fimin);%将pai1减去pad1的差并且除以(0.0001+fimax-fimin)所得的值
        %与第二行到第NS-1行的第C+1列的数相加 并将值赋给第二行到第NS-1行的第C+1列即 pa(2:(Ns-1),C+1)
        %所有循环结束后则表示对于每个目标函数得出的每列值中的一个I I的距离为其后一个减去其前一个
    end
    padis=pa(:,C+1); %将第C+1列的数值赋给
    pa=pa(:,1:C); %将第一列到第C的数值赋给pa
    POP=POP;%此句个人觉得没有很大的意义
    %padis=sum(pa,2);
    [aa,ss]=sort(-padis);%-padis进行从小到大排序排序则亦是对padis进行从大到小的排序 并且将排序好赋给aa 并将其原来所对应的位置赋给ss
    %[aa,ss]=sort(padis);
    EPOP=POP(ss(1:NM),:);%按照上面得到的padis从大到小的顺序号ss将POP进行重新排列取到第一行到第NM行 并且保存在数组EPOP中 作为下一代
    Epa=pa(ss(1:NM),:);%按照上面得到的padis从大到小的顺序号ss 将pa进行重新排列取到第一行到第NM行 并且保存在数组Epa中 作为下一代
else
    EPOP=POP;Epa=pa;  %若对于得到的群体数目不大于NM时 则可以直接全部选中作为下一代
end
time=toc; %与函数开始设置的tic相对应 tic用在程序的开始，作用是启动一个计时器，然后在程序尾部放一个toc，表示终止计时器，并返回tic启动以来的总时s间
%%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [NPOP,time]=Clonef(POP,pa,CS)  %克隆函数的定义
%--------------------------------------------------------------------------
tic   %Clonef函数运行开始的计时器启动
NC=[]; %初始化一个数组NC
[N,C]=size(POP); %将解数组POP中的行号与列号分别赋给N和C
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
NPOP=[]; % 将初始化定义一个新数组NPOP
FitPop=[POP padis];%将解集和距离结合赋给数组FitPop
for i=1:N %将数组循环
    NiPOP=ones(NC(i),1)*FitPop(i,:); %根据NC得到每个抗体的的克隆数目
    NPOP=[NPOP;NiPOP]; %将之前的解和克隆后的解按行组合在一起从NPOP的最后一行开始添加NiPOP 
end
time=toc;  %与函数开始设置的tic相对应 tic用在程序的开始，作用是启动一个计时器，然后在程序尾部放一个toc，表示终止计时器，并返回tic启动以来的总时s间
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [POP,pa,padis]=CDAf(tPOP,tpa) %求距离的函数 与前面选择下一代种群比较类似
%--------------------------------------------------------------------------
[Ns,C]=size(tpa); %将tpa的行号和列号分别赋给Ns和C 而tpa是作为目标函数值
padis=[];  %初始化一个padis的空数组
for i=1:C   %从列号进行循环
    [k,l]=sort(tpa(:,i)); %将第i列从小到大快速排序将排序后的顺序赋给k 并将其所对应原来的行号赋给l
    N=[];M=[]; %初始化两个空数组
    N=tpa(l,:);M=tPOP(l,:);%将tpa按照l的顺序的值存储在N中  tPOP按照l的顺序的值存储在M
    tpa=[];tPOP=[];%将tpa和tPOP赋值空
    tpa=N;tPOP=M;%再将N和M分别赋给tpa和tPOP  
    %以上几个步骤的作用就是为了将tpa按从小到大排序 并将排序后的顺序保存在原来的数组tpa和PtOP中
    tpa(1,C+1)=Inf;tpa(Ns,C+1)=Inf; %将第一行和第NS行 第C+1列的数值赋值为无穷大
    tpai1=[];tpad1=[];%初始化两个数组tpai1和tpad1
    tpai1=tpa(3:Ns,i);tpad1=tpa(1:(Ns-2),i);%将第3行到第NS行的第i列的所有数值赋给tpai1 将第1行到第NS-2行的第i列的所有数值赋给tpad1
    fimin=min(tpa(:,i));fimax=max(tpa(:,i));%找到tpa第i列中的最大和最小 分别赋给max和min
    tpa(2:(Ns-1),C+1)=tpa(2:(Ns-1),C+1)+(tpai1-tpad1)/(0.0001+fimax-fimin);%将tpai1减去tpad1的差并且除以(0.0001+max-min)所得的值
     %与第二行到第NS-1行的第C+1列的数相加 并将值赋给第二行到第NS-1行的第C+1列即 tpa(2:(Ns-1),C+1)
      %所有循环结束后则表示对于每个目标函数得出的每列值中的一个I I的距离为其后一个减去其前面一个
end
pa=tpa(:,1:C);%将tpa的第一列到第C列中所有的数据送给pa
POP=tPOP;%将tPOP中所有的值送给POP
padis=tpa(:,C+1);  %将tpa的第C+!行数据送给padis  因为所有的距离都保存在第C+1列中
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NPOP,time]=Recombinationf(POP,bu,bd,CPOP,it)  %对重组函数的定义
%--------------------------------------------------------------------------
tic;  %Recombinationf函数运行开始的计时器启动
eta_c=15;  %相对应于论文中的控制父代个体交叉程度的数量级
[N,C1]=size(POP);%将解集POP的行号与列号分别赋值给N和C1
Fit=POP(:,C1); %将解集数组POP的第C1列赋值给Fit中
POP=POP(:,1:C1-1); %将解集POP中第一列到第C1-1列中的数据赋给POP 即将原来的POP去除最后一列的情况
C=C1-1; %用C来表示C1-1
NPOP=POP; %将POP的值同时保存在数组NPOP
for i=1:N  %进行对POP的每行进行循环
    r1=rand; %随机产生一个0到1之间的数送给r1
    if r1<=1.1 % 这一步是一个肯定成立的过程 为什么要这样设计 有点不明白
        aa=randperm(size(CPOP,1));%randperm是matlab中自带的函数 randperm（x）表示随机生成从1到X的所有整数
        bb=aa(1);%将aa数组的第一个数赋给bb
        for j=1:C %对POP的列号开始进行循环
            par1=POP(i,j);%将第i行第j列的数据赋给par1
            par2=CPOP(bb,j);%将第bb行第j列的值赋给par2
            yd=bd(j);yu=bu(j);%将上界bu和下界bd的第j列赋值给yu和yd
            r2=rand; %随机产生0到1的随机数并将值赋给r2
            if r2<=0.5 %对r2进行判断 若r2小于0.5进行下面步骤
                if abs(par1-par2)>10^(-14)%若par1和par2的差的绝对值大于10的负14次方
                    y1=min(par1,par2);y2=max(par1,par2); %找到par1和par2中小的赋给y1 大的值赋给y2
                    if (y1-yd)>(yu-y2)%如果y1减去下界yd的值大于上界yu减去y2
                        beta=1+2*(yu-y2)/(y2-y1); %则分子取yu减去y2
                    else
                        beta=1+2*(y1-yd)/(y2-y1);%否则分子y1减去yd
                    end
                    expp=eta_c+1;beta=1/beta;alpha=2.0-beta^(expp); %与论文的公式相对应 由此得出beta的值小于1 则相对应的beta^(expp)小于1
                    % 而得到alpha的范围是在1到2之间  怎1/alpha的范围是0.5到1之间
                    r3=rand; %再次随机产生一个0到1的数值赋给r3
                    if r3<=1/alpha  %跟产生一个平均分布的随机数有什么不同呢?  这里alpha的值要受到不同的上下界的变化而变化 但其所有的范围都在0.5到1之间
                        alpha=alpha*r3;expp=1/(eta_c+1.0);
                        betaq=alpha^(expp);  %即论文中的控制交叉参数
                    else
                        alpha=1/(2.0-alpha*r3);expp=1/(eta_c+1);
                        betaq=alpha^(expp);%在产生的随机数大于1/alpha时 则betaq等于alpha^(expp)上面几步以及参数的设置 可以得到结论是当参数eta_c的值越大时
                        %重组生成的子代就以较高的概率更加接近父代 因为根据所给出的公式 我们可以判断出当eta_c越大是
                        %则相应的alpha就越大 因此其betaq就越大 但是betaq的值是小于1的 eta_c越大则betaq就更加向1接近
                    end  %父代抗体之间相近则交叉生成的子代越相似,
                    chld1=0.5*((y1+y2)-betaq*(y2-y1));%可以化简做chid1=0.5*（（1+betaq）*y1+（1-betaq）*y2）
                    chld2=0.5*((y1+y2)+betaq*(y2-y1)); %可以化简做chid1=0.5*（（1+betaq）*y2+（1-betaq）*y1）
                    %上面几步以及参数的设置 可以得到结论是当参数eta_c的值越大时重组生成的子代就以较高的概率更加接近父代 因为根据所给出的公式 我们可以判断出当eta_c越大是
                        %则相应的alpha就越大因此其betaq就越大当betaq越大时则根据上面的化简可以明确chil1的值趋向于y1 chil2的值趋向于y2
                    if rand<=0.5 %随机产生一个0到1范围的饿随机数 看其是否小于或等于0.5 
                        aa=max(chld1,yd);%将重组产生的子代chil1与下界进行比较 将大的值送给aa
                        NPOP(i,j)=min(aa,yu);%将aa与上界进行比较 将小的那个值赋给NPOP的第i行第j列
                    else
                        aa=max(chld2,yd);%将重组产生的子代chil2与下界进行比较 将大的值送给aa
                        NPOP(i,j)=min(aa,yu);%将aa与上界进行比较 将小的那个值赋给NPOP的第i行第j列
                    end
                    %这几步的步骤主要是重组产生的子代根据随机的概率 看将产生的两个子代中那个子代来取代原来解集中的NPOP（i，j）
                end  
            end
        end
    end
end
NPOP=[NPOP Fit];
time=toc;  %与函数开始设置的tic相对应 tic用在程序的开始，作用是启动一个计时器，然后在程序尾部放一个toc，表示终止计时器，并返回tic启动以来的总时s间

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NPOP,time]=Mutationf(POP,bu,bd,pm,eta_m,ppm,it) %混合变异函数的定义
%--------------------------------------------------------------------------
tic;  %Mutationf函数运行开始的计时器启动
[N,C1]=size(POP); %将解集POP的行号和列号分别赋给N和C1
Fit=POP(:,C1); %将POP中的第C1列赋给Fit 感觉对后面没有起到什么作用
POP=POP(:,1:C1-1); %将POP中的第一列到第C1-1列赋给POP
%eta_m=20;
NPOP=POP; %同时将POP的数据保存在数组NPOP中
C=C1-1; %让C=C1-1
for i=1:N  % 进行循环从第一行到第N行的循环
    r2=rand(); %随机产生一个0到1的随机数
    ppm=(0.3-0.2*it/250)*(Fit(i)-min(Fit))/(max(Fit)-min(Fit)); %相当于论文中的aspi 作为一个阀值来控制两个变异算子的交换
    if r2>ppm %当随机数r2大于这个阀值时 则进行执行多项式变异
    for j=1:C %从POP的第一列到倒数第二列进行循环
        r1=rand; %随机产生一个0到1之间的数送给r1
        if r1<=pm %如果r1小于代码开头给出的变异概率pm
            y=POP(i,j);%将POP（i，j）中的数值送给y
            yd=bd(j);yu=bu(j);%将上界bu和下界bd的第j列赋值给yu和yd
            if  y>yd %如果y大于上界 则进入下面的运算
                if (y-yd)<(yu-y) %若y减去下界yd的值小于上界减去y的值
                    delta=(y-yd)/(yu-yd);%选择小的值除以上界和下界的差
                else
                    delta=(yu-y)/(yu-yd);%选择小的值除以上界和下界的差
                end
                r2=rand; %在随机产生0到1之间的数并将其值赋给r2
                indi=1/(eta_m+1); %eta_m作为变异分布参数
                if r2<=0.5 
                    xy=1-delta;
                    val=2*r2+(1-2*r2)*(xy^(eta_m+1));
                    %val=2*r2;
                    deltaq=val^indi-1;
                else
                    xy=1-delta;
                    val=2*(1-r2)+2*(r2-0.5)*(xy^(eta_m+1));
                    %val=2*(1-r2);
                    deltaq=1-val^indi;
                end
                %针对随机数的值在界点0.5的范围段来分别s是否用的多项式变异公式或者高斯变异 控制变异特征 平均来说eta_m值越大时产生的变异尺度越小
                y=y+deltaq*(yu-yd);%抗体根据变异概率与上下界之差相乘再与原来的抗体相加得到的值就是变异后的抗体
                NPOP(i,j)=min(y,yu); %变异后的抗体与上界进行比较 将小的值取代原来的NPOP（i，j）
                y=NPOP(i,j);
                NPOP(i,j)=max(y,yd);%变异后的抗体与下界进行比较 将大的值取代原来的NPOP（i，j）
            else
                NPOP(i,j)=rand*(yu-yd)+yd;%针对上面if y>yd不成立的情况 要重新产生随机数
            end
        end
    end
 else
    for j=1:C% 针对if r2>ppm不成立时 要进行高斯变异
        r1=rand;%随机产生一个0到1之间的数并将数送给r1
        if r1<=pm %判断若r1小于或等于变异概率pm
            y=POP(i,j); %将POP（i，j）中的数值送给y
            yd=bd(j); yu=bu(j); %将上界bu和下界bd的第j列赋值给yu和yd
            y=y+0.1*normrnd(0,1)*(yu-yd);%normrnd(0,1)为一个平均值为0，标准差为1得瑟高斯随机数 
            NPOP(i,j)=min(y,yu);%变异后的抗体与上界进行比较 将小的值取代原来的NPOP（i，j）
            y=NPOP(i,j);
            NPOP(i,j)=max(y,yd);%变异后的抗体与下界进行比较 将大的值取代原来的NPOP（i，j）
            %NPOP(i,j)=rand*(yu-yd)+yd;
        end
    end
    end
end
time=toc; 
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NPOP,time]=Mutationf2(POP,bu,bd,pm,eta_m)%此时为什么又定义一个混合变异的函数 相对于上面的函数没有参数ppm和it 这里表示的是随机变异的结果
%--------------------------------------------------------------------------
tic;  %Mutationf2函数运行开始的计时器启动
[N,C]=size(POP);%将解集POP的行号和列号分别赋给N和C
%eta_m=20;
NPOP=POP;%将POP的数据保存在数组NPOP中
for i=1:N %从第一行到第N行进行外部的循环
    for j=1:C %从第一列到第C列进行内循环
        r1=rand;%随机产生一个0到1之间的随机数 并将值赋给r1
        if r1<=pm %判断若r1小于或等于变异概率pm
            y=POP(i,j); %将POP（i，j）中的数值送给y
            yd=bd(j);yu=bu(j); %将上界bu和下界bd的第j列赋值给yu和yd
            NPOP(i,j)=rand*(yu-yd)+yd; %随机产生抗体将其赋给 NPOP(i,j)
        end
    end
end
time=toc; 
