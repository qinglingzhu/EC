%% ��¡�����Ķ���
function  [NPOP,time]=Clonef(POP,pa,CS)  
tic   
NC=[];
[N,C]=size(POP); 
[POP,pa,padis]=CDAf(POP,pa);%���������CDAf���������Ⱥ֮�����,���ݾ�����п�¡
aa=find(padis==inf); %����padis������������� ��padis�ĵ�һ�к͵�Ns��
bb=find(padis~=inf); %����padis��������������� 
if length(bb)>0  %������padis�����������������
    padis(aa)=2*max(max(padis(bb))); %�ѱ߽��Ŀ����׺϶������趨Ϊ���߽��������׺϶ȵ����� ���߽��ľ�������Ϊ���߽���������������
    NC=ceil(CS*padis./sum(padis));%ceil������matlab�Դ��ĺ��� ceil��X�����ʾ����X����С���� CSΪԤ�����¡���� �����׺϶ȣ���ӵ�����룩
    %������ÿ����Ŀ�¡���� �����������Ӧ �׺϶����¡���������� �׺϶�Խ�� ��������Խ�����¡�ĸ���Խ�� �Դ˺ܺõı�����Ⱥ�Ķ�����
else
    NC=ceil(CS/length(aa))+zeros(1,N); %���������bb���ȴ�����ʱ ���¡�ĸ���Ϊ�ܵĿ�¡��Ŀ����aa�ĳ��� ��ƽ������ÿ������Ŀ�¡��Ŀ
end
NPOP=[]; %����ʼ������һ��������NPOP
FitPop=[POP padis]; %���⼯�;����ϸ�������FitPop
for i=1:N %������ѭ��
    NiPOP=ones(NC(i),1)*FitPop(i,:); %����NC�õ�ÿ������ĵĿ�¡��Ŀ
    NPOP=[NPOP;NiPOP]; %��֮ǰ�Ľ�Ϳ�¡��Ľⰴ�������һ���NPOP�����һ�п�ʼ���NiPOP 
end
time=toc; 