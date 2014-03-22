function HIMO()
EMOinstruction;%% display the instruction for running the programming.
%--------------------------------------------------------------------------
TestNO=input('press the enter key after inputting the serial number of test problem:');
Trial=input('input the number of independent runs:');
%p=0.2;
NA=20;                                                      % �Ⱥ��ĸ��� 
CS=100;                                                     % ��¡����  
NM=100;                                                     % ֧��Ⱥ��ĸ���(�ѽ��ŵĽ�������Ⱥ������)
gmax=250;                                                   % ���������� 
[bu,bd,testfunction]=getbud(TestNO);                        % bu ����������Ͻ�;bd ����������½� ������һ������getbud.m�ĵ��� TestNO����Ϊ�û������ѡ��ֵ
c=size(bu,2);                                               % ��ȡ�����ĸ��� size��x��2����ʾ��bu��������������Ԫ�� ��һԪ���� ���Ƕ�Ԫ���� ��������Ԫ������
pmy=1/c;                                                    % ���ñ����ı������ 
p=0.2;                                                      % Ԥ������� �����ı�����صĲ���
if bu==bd 
    return;
end                                        % ���Ͻ���½����ʱֹͣ�㷨  
%--------------------------------------------------------------------------
Datime=date;        %��ȡ������                                      
Datime(size(Datime,2)-4:size(Datime,2))=[];%����ʱ�� ����һ�н�ϵõ�������
TestTime=clock;%����ʱ��
TestTime=[num2str(TestTime(4)),'-',num2str(TestTime(5))];%����һ�н�ϵõ�ʱ�����ʾ ��ʱ���
Method='HIMO';
paretof=[]; %�洢��õķ��ӽ� ��Ϊ�洢����
runtime=[];%���Ƚ����еĴ�����ʼ��Ϊ��
for trial=1:Trial %��Ϊһ����ѭ�� ��Ҫ���ж��ٴε�����������������ӦTrial=input('input the number of independent runs:');
                   %Trial��HIMO����������ʱ�û������ֵ ������ʱtime�ļ��� �� ����trial��250������
    timerbegin=clock; %��ʼ��ʱ
%--------------------------------------------------------------------------
POP=rand(NM,c).*(ones(NM,1)*bu-ones(NM,1)*bd)+ones(NM,1)*bd;%���������Ⱥ  ����������Ⱥ�ķ�Χ���û�����getbud.m�����е�ѡ��ֵ
%ones(NM,1)����Ϊ����NM�� һ�С�1�� rand(NM,c)��������Ϊ�������NM��C��0��1֮�����                                                           
ME=[];%Initialization������ME��ʼ��
%--------------------------------------------------------------------------
pa=OVcom(POP,TestNO); %���ض�Ŀ�꺯��ֵ �Ժ���OVcom.m�ĵ��� ����������ɵ���Ⱥ���д��뵽���Ժ����������ص�ֵ��������pa
[DON,DONt]=IDAf(pa); %��������Сֵ��˵1:֧��⣬0��ʾ��֧���
 %�õ���һ�����
nodom=find(DON==1); 
MEpa=pa(nodom,:);MEPOP=POP(nodom,:);
numnod=size(MEPOP,1); %���ӽ����
[ClonePOP,Clonepa,RNDCDt]=UDPf(MEPOP,MEpa,NA);%������ӵ���̶�����,��ɻ��Ⱥ �Ժ���UDPf�����ĵ��� 
 %�����ж�����ӽ����������õ��Ľ�ĸ���������Ȼ����NM(֧��Ⱥ����ܸ���) ���Ǻ�Ҫ����ѡ����һ�����п�¡����Ŀ
it=0;Cloneti=0;DASti=0;PNmti=0;DONjudti=0;RNDCDti=0;AbAbfitcomti=0;
%it �������� Cloneti:��¡�����ѵ�ʱ�� DASti:���������ѵ�ʱ�� 
%PNmti:���������ѵ�ʱ�� DONjudti:�ж�֧���ϵ�����ѵ�ʱ�� RNDCDti:ѡ�������ѵ�ʱ��
%AbAbfitcomti:����û�ٳ���,��������
while it<gmax %it��Ϊ���д���  gmax��Ϊ�������д��� ��Ϊ����ѭ������
%--------------------------------------------------------------------------
eta_m=20; %��������������GP-HMģ���������õı���ֲ�ָ��n
pm=-1*p*pmy*(it/125)+(1+p)*pmy;   %�൱�������е�pm=��1+p��*pmy-2*p*pmy*��it/250��  p�Ǹտ�ʼ�����ı����ı������
%���Ǵ˾����е���������������һ����������it<gmax/2   ���ñ����ı������
if pm<pmy
    pm=pmy;  %�����е�pmy����ʾ����Ԥ������С�������
end
ppm=0.2; %��������Ƿ�Ӧ����pmy ���������pm���޷�����  % �ƺ��ں���Ļ�ϱ���������в����г�����ppm �񶨵������뷨
cloneover=[]; %�Կ�¡��Ⱥ��������г�ʼ��
[cloneover,Clonet]=Clonef(ClonePOP,Clonepa,CS);  %�Ի��Ⱥ������¡,���ؿ�¡��Ⱥ�Ͳ���ʱ��
[cloneover,DASt]=Recombinationf(cloneover,bu,bd,ClonePOP,it);%�Կ�¡��Ⱥ����ģ������ƽ���,���ؽ�����¡��Ⱥ�Ͳ���ʱ��
[cloneover,PNmt]=Mutationf(cloneover,bu,bd,pm,eta_m,ppm,it); %�Խ������Ⱥ���ж���ʽ����,���ؿ�¡��Ⱥ�Ͳ���ʱ��
%--------------------------------------------------------------------------
clonepa=OVcom(cloneover,TestNO); %�Ժ���OVcom�ĵ��û�ȡ������¡��Ⱥ��Ŀ�꺯��ֵ
NPOP=[MEPOP;cloneover]; %������¡��Ⱥ��ԭ�����ӽ������һ����Ⱥ
Npa=[MEpa;clonepa];%������¡��Ⱥ��Ŀ�꺯��ֵ��ԭ�����ӽ��Ŀ�꺯��ֵ�����һ���͸�Npa
[NDON,DONjudt]=IDAf(Npa);%�Ժ���IDAf�ĵ��ý�����˵�Ŀ�꺯��ֵNpa�����жϷ��ӽ�
Nnodom=find(NDON==1);%����һ���ҵ��ķ��ӽ���в鵽���ӽ���кŲ����ظ���Nnodom
NEpa=Npa(Nnodom,:);NEPOP=NPOP(Nnodom,:);%���Ʒ���Ⱥ��  �����ӽ�ͷ��ӽ�����Ӧ��Ŀ�꺯��ֵ�ֱ𸳸���������NEPOP��NEpa
Nnumnod=size(NEPOP,1);%��÷���Ⱥ��ĸ��� size(NEPOP,1)����matlab�Դ�����NEPOP���к�
[MEPOP MEpa,RNDCDt]=UDPf(NEPOP,NEpa,NM);%�Ժ���UDPf�ĵ��ý���ѡ����һ����Ⱥ %���ڷ��ӽ����NM��֧��Ⱥ����ܸ�����ʱҪ���еĽ�һ����ѡ��
%[MEPOP MEpa]=clustering(NEPOP,NEpa);    ��ѡ�����һ����Ⱥ�ķ��ӽ�������Ӧ��Ŀ�꺯��ֵ�ֱ𸳸�MEPOP��MEpa
numnod=size(MEPOP,1);%ѡ��ķ���Ⱥ��ĸ��� size(MEPOP,1������matlab�Դ� ����MEPOP���к�
[ClonePOP,Clonepa,RNDCDt]=UDPf(MEPOP,MEpa,NA);%ѡ���¡��Ⱥ NA��Ϊ���Ⱥ����Ŀ��ѡ��NA���ڵĵ���Ⱥ��Ŀ���п�¡
%Update Dominant Population
%--------------------------------------------------------------------------
it=it+1; %�������whileѭ�����Ӧ ��it<gmaxʱѭ��һֱ������ȥ
Cloneti=Cloneti+Clonet;%��ÿһ����¡�����ѵ�ʱ�����
DASti=DASti+DASt; %��ÿһ�����������ѵ�ʱ�����
PNmti=PNmti+PNmt;%��ÿһ�����������ѵ�ʱ�����
DONjudti=DONjudti+DONjudt;%��ÿһ�������жϷ��ӽ⣨���ж�֧���ϵ����ʱ����� 
RNDCDti=RNDCDti+RNDCDt;%������ÿһ��ѡ�������ѵ�ʱ���ۼ�����
disp(sprintf('time: %d   generation: %d    number of nodominate:  %d',trial,it,numnod));%time��֮ǰ�û������trial 
%generation����Ϊ���еĴ��� ����������Ϊ250�� nodominateΪ��֧������Ŀ ������������ΪNM��100��
end  %the end of iterations ����ѭ��
%--------------------------------------------------------------------------
%Save the output solutions
[NS(trial),NF]=size(MEpa);%size(MEpa)Ϊ���MEpa���к����кŷֱ𸳸�NS��trial����NF
Trials=trial; 
paretof=[paretof;MEpa]; %paretof��Ϊ�洢���ӽ�
runtime(trial)=etime(clock,timerbegin); % etime(t1,t2)��������������������t1��t2֮���ʱ��� ���ڴ�Ϊ����clock��timerbegin֮���ʱ��� 
%�˴���ʱ��Ϊ�������������ʱ��
Clonetime(trial)=Cloneti;%�����д���¡�����ܵ�ʱ�丳��Clonetime(trial)
DAStime(trial)=DASti;%�����д������������ܵ�ʱ�丳��DAStime(trial)
PNmtime(trial)=PNmti; %�����д������������ܵ�ʱ�丳��PNmtime(trial)=PNmti;
DONjudtime(trial)=DONjudti;%�����д��жϷ��ӽ⣨���ж�֧���ϵ�������ѵ��ܵ�ʱ�丳��DONjudtime(trial)
RNDCDtime(trial)=RNDCDti;%�����д�ѡ�������ѵ��ܵ�ʱ�丳��RNDCDtime(trial)
eval(['save ', testfunction Method Datime TestTime ,' Method gmax NA CS paretof runtime Trials NS NF TestTime Datime testfunction ']) ;
end  %the end of runs  
%--------------------------------------------------------------------------
Frontshow(MEpa);% plot the Pareto fronts solved by the last run �Ժ�����������ĵ��� ���Ŀ�꺯��ֵ��ͼ�ηֲ�










