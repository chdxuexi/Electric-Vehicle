function [record,fym]=PSO_Car_3_13()
clc;
clear;
close all;
%% ��ʼ����Ⱥ�Ĳ�������
tic;%��ʱ��
carNum = 1200;
d = carNum;                      % 100������100ά
pileNum = 222;                   % 222�����׮
N = 150;                         % ��ʼ��Ⱥ������ÿ�����Ӷ���һ�����ȷ�����
ger = 300;                       % ����������
[EV_arrive_time,EV_charging_time] = DataSource;   %���������������׮��ʱ������
Wmax=0.9;                        %����Ȩ�����ֵ
Wmin=0.4;                        %����Ȩ����Сֵ
xlimit = [1, pileNum];           % ����λ�ò�������(�������ʽ���Զ�ά)
xlimit = repmat(xlimit,d,1);
vlimit = [-10, 10];              % �����ٶ�����
vlimit = repmat(vlimit,d,1);
c_2 = 1.49;                      % ����ѧϰ����
c_3 = 1.49;                      % Ⱥ��ѧϰ����
for i = 1:d
    x = xlimit(i, 1) + (xlimit(i, 2) - xlimit(i, 1)) * rand(N, d);    %��ʼ��Ⱥ��λ��
end
x = round(x);                    % �Գ�ʼ��Ⱥ��λ�ý���ȡ������������
v = rand(N, d);                  % ��ʼ��Ⱥ���ٶ�
xm = x;                          % ÿ���������ʷ���λ��
ym = zeros(1, d);                % ��Ⱥ����ʷ���λ��
fxm = zeros(N, 1);               % ÿ���������ʷ�����Ӧ��
fym = inf;                       % ��Ⱥ��ʷ�����Ӧ��
%% ����Ⱥ����
iter = 1;
times = 1;
record = zeros(ger, 1);          
record_syd = 0;
while iter <= ger 
    c_1 = Wmax-(Wmax-Wmin)*(iter/ger)^2;
    syd = [];     %ÿһ�����ӵ���Ӧ��ֵ
    for i=1:N     %NΪ��Ⱥ����
        DX=x(i,:);
        T=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);    %���ú�����������Ӧ��ֵ
        syd=[syd,T];
    end
    fx = syd;                       % ��ǰ������N�����ӵ���Ӧ��ֵ   
    %��iter=1��ʱ��Ҫ�����fxm��xm��ʼ����ֱ�Ӹ�ֵ������
    if iter ==1
        for i = 1:N                 % N = 150 ��Ⱥ����
            if fxm(i) < fx(i)
                fxm(i) = fx(i);     % ���¸�����ʷ�����Ӧ��
                xm(i,:) = x(i,:);   % ���¸�����ʷ���λ��
            end
        end
        %��iter>1��ʱ������е�ǰfx��fxm��xm���бȽϽ����Ĳ���
    else
        for i = 1:N                % N = 150 ��Ⱥ����
            if fxm(i) > fx(i)
                fxm(i) = fx(i);    % ���¸�����ʷ�����Ӧ��
                xm(i,:) = x(i,:);  % ���¸�����ʷ���λ��
            end
        end
    end 
    if fym > min(fxm)              %fym�ĳ�ʼֵ=inf
        [fym, nmax] = min(fxm);    % ����Ⱥ����ʷ�����Ӧ��
        ym = xm(nmax, :);          % ����Ⱥ����ʷ���λ��
    end   
    v =  c_1*v  + c_2 * rand *(xm - x) + c_3 * rand *(repmat(ym, N, 1) - x);% �ٶȸ���   
    % �߽��ٶȴ���
    for i=1:d
        for j=1:N
            if  v(j,i)>vlimit(i,2)
                v(j,i)=vlimit(i,2);
            end
            if  v(j,i) < vlimit(i,1)
                v(j,i)=vlimit(i,1);
            end
        end
    end
    v = round(v);
    % λ�ø���
    x = x + v;
    % �߽�λ�ô���
    for i=1:d
        for j=1:N
            if  x(j,i)>xlimit(i,2)
                x(j,i)=xlimit(i,2);
            end
            if  x(j,i) < xlimit(i,1)
                x(j,i)=xlimit(i,1);
            end
        end
    end  
    record(iter) = fym;            %��Сֵ��¼
    iter = iter+1;
    times=times+1;
end
toc;

%% ������
% figure(1)
% plot(record,'Linewidth',1.5);
% xlabel('Number of iterations','Linewidth',1.5) ;
% ylabel('Fitness value','Linewidth',1.5) ;
% % title('��Ӧ�Ƚ�������');
% 
% disp(['��С��ʱ��',num2str(fym)]);
% disp(['�������Ƚ����',num2str(ym)]);
end
