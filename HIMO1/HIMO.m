function HIMO()
EMOinstruction;%% display the instruction for running the programming.
%--------------------------------------------------------------------------
TestNO=input('press the enter key after inputting the serial number of test problem:');
Trial=input('input the number of independent runs:');
%p=0.2;
NA=20;                                                      % 活动群体的个数 
CS=100;                                                     % 克隆总数  
NM=100;                                                     % 支配群体的个数(把较优的解放在这个群体里面)
gmax=250;                                                   % 最大进化代数 
[bu,bd,testfunction]=getbud(TestNO);                        % bu 代表变量的上界;bd 代表变量的下界 对另外一个函数getbud.m的调用 TestNO是作为用户输入的选择值
c=size(bu,2);                                               % 获取变量的个数 size（x，2）表示就bu的列数即变量的元数 是一元函数 还是二元函数 或者是三元函数等
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
POP=rand(NM,c).*(ones(NM,1)*bu-ones(NM,1)*bd)+ones(NM,1)*bd;%随机生成种群  这种生成种群的范围是用户对于getbud.m函数中的选择值
%ones(NM,1)是作为生成NM行 一列“1” rand(NM,c)函数是作为随机生成NM行C列0到1之间的数                                                           
ME=[];%Initialization对数组ME初始化
%--------------------------------------------------------------------------
pa=OVcom(POP,TestNO); %返回多目标函数值 对函数OVcom.m的调用 即将随机生成的种群进行代入到测试函数中所返回的值赋给数组pa
[DON,DONt]=IDAf(pa); %对于求最小值来说1:支配解，0表示非支配解
 %得到第一层端面
nodom=find(DON==1); 
MEpa=pa(nodom,:);MEPOP=POP(nodom,:);
numnod=size(MEPOP,1); %非劣解个数
[ClonePOP,Clonepa,RNDCDt]=UDPf(MEPOP,MEpa,NA);%按距离拥挤程度排列,组成活动子群 对后面UDPf函数的调用 
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
NPOP=[MEPOP;cloneover]; %变异后克隆子群与原来非劣解组成下一代种群
Npa=[MEpa;clonepa];%变异后克隆子群的目标函数值与原来非劣解的目标函数值组后在一起送给Npa
[NDON,DONjudt]=IDAf(Npa);%对函数IDAf的调用将结合了的目标函数值Npa进行判断非劣解
Nnodom=find(NDON==1);%将上一步找到的非劣解进行查到非劣解的行号并返回赋给Nnodom
NEpa=Npa(Nnodom,:);NEPOP=NPOP(Nnodom,:);%复制非劣群体  将非劣解和非劣解所对应的目标函数值分别赋给两个数组NEPOP和NEpa
Nnumnod=size(NEPOP,1);%获得非劣群体的个数 size(NEPOP,1)函数matlab自带返回NEPOP的列号
[MEPOP MEpa,RNDCDt]=UDPf(NEPOP,NEpa,NM);%对函数UDPf的调用进行选择下一代种群 %即在非劣解大于NM（支配群体的总个数）时要进行的进一步的选择
%[MEPOP MEpa]=clustering(NEPOP,NEpa);    将选择的下一代种群的非劣解与所对应的目标函数值分别赋给MEPOP和MEpa
numnod=size(MEPOP,1);%选择的非劣群体的个数 size(MEPOP,1）函数matlab自带 返回MEPOP的列号
[ClonePOP,Clonepa,RNDCDt]=UDPf(MEPOP,MEpa,NA);%选择克隆子群 NA作为活动种群的数目即选择NA以内的的种群数目进行克隆
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
end  %the end of runs  
%--------------------------------------------------------------------------
Frontshow(MEpa);% plot the Pareto fronts solved by the last run 对后面输出函数的调用 输出目标函数值的图形分布










