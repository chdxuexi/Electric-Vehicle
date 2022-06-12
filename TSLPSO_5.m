%% ������TSLPSO_5  100����600��
function [record,fym]=TSLPSO_5()
clc;
clear;
close all;
%% ��ʼ����Ⱥ�Ĳ�������
carNum = 600;
d = carNum;                      % 600������600ά
pileNum = 222;                   % 222�����׮
N = 150;                         % ��ʼ��Ⱥ������ÿ�����Ӷ���һ�����ȷ�����
ger = 300;                       % ����������
Wmax=0.9;                        %����Ȩ�����ֵ
Wmin=0.4;                        %����Ȩ����Сֵ
c_2 = 1.49;                       % ����ѧϰ����
c_3 = 1.49;                       % Ⱥ��ѧϰ����
[EV_arrive_time,EV_charging_time] = DataSource; %���������������׮��ʱ������
% ��ʼ���ٶ�v x
xmin = 1;                  % ����λ�ò�������(�������ʽ���Զ�ά)
xmax = pileNum;
vmin = -10;                % �����ٶ�����
vmax = 10;
for i0=1:N
    for j0 = 1:d
        x(i0,j0) = xmin + (xmax - xmin) * rand; %��ʼ��Ⱥ��λ��
        v(i0,j0) = xmin + (xmax - xmin) * rand; %��ʼ��Ⱥ���ٶ�
    end
end
x = round(x);                    % �Գ�ʼ��Ⱥ��λ�ý���ȡ������������
syd = [];                        %ÿһ�����ӵ���Ӧ��ֵ
for i=1:N                        %NΪ��Ⱥ����
    DX=x(i,:);
    T=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);    %���ú�����������Ӧ��ֵ
    syd=[syd,T];
end
fx = syd;                        % ��ǰ������N�����ӵ���Ӧ��ֵ
xm = x;                          % ÿ���������ʷ���λ��
ym = zeros(1, d);                % ��Ⱥ����ʷ���λ��
fxm = zeros(N, 1);               % ÿ���������ʷ�����Ӧ��
fym = inf;                       % ��Ⱥ��ʷ�����Ӧ��
for i = 1:N                      % N = 150 ��Ⱥ����
    fxm(i) = fx(i);     % ���¸�����ʷ�����Ӧ��
    xm(i,:) = x(i,:);   % ���¸�����ʷ���λ��
end
if fym > min(fxm)             %fym�ĳ�ʼֵ=inf
    [fym, nmax] = min(fxm);   % ����Ⱥ����ʷ�����Ӧ��
    ym = xm(nmax, :);         % ����Ⱥ����ʷ���λ��
end

%% ����ѭ��
iter = 1;
xdl=xm;%�����ֲ�����ֵ
record = zeros(ger, 1); 
while iter <= ger % ger��������
    %����Ⱥ1
    c_1 = Wmax-(Wmax-Wmin)*(iter/ger)^2;
    for i = 1:N
    v(i,:) =  c_1*v(i,:)+c_2 * rand *(xdl(i,:) - x(i,:)) + c_3 * rand *(ym - x(i,:));     
    % �߽��ٶȴ���
        for j=1:d
                if  v(i,j)>vmax
                    v(i,j)=vmax;
                end
                if  v(i,j) < vmin
                    v(i,j)=vmin;
                end
        end
        v = round(v);
        % λ�ø���
        x(i,:) = x(i,:) + v(i,:);
        % �߽�λ�ô���
        for j=1:d
                if  x(i,j)>xmax
                    x(i,j)=xmax;
                end
                if  x(i,j) < xmin
                    x(i,j)=xmin;
                end
        end       
        %��Ӧ�ȼ���
        DX=x(i,:);
        fx(i)=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);%���ú�����������Ӧ��ֵ   
        if fx(i)<fxm(i)
            fxm(i) = fx(i);     % ���¸�����ʷ�����Ӧ��
            xm(i,:) = x(i,:);   % ���¸�����ʷ���λ��
            flag(i)=1;
        else
            flag(i)=0;
        end
        if flag(i)==1
            %�㷨2��ʼ�����ɷ���xdl��ֵ
            xdl(i,:)=xm(i,:);
            for j=1:d
                temp=xdl(i,:);
                if temp(j)==ym(j)
                    continue;
                end
                temp(j)=ym(j);
                temp_f=Car_Fitness(temp,N,EV_arrive_time,EV_charging_time);
                xdl_f=Car_Fitness(xdl(i,:),N,EV_arrive_time,EV_charging_time);
                if temp_f<xdl_f
                    xdl(i,j)=ym(j);
                end
            end%�㷨2����
        end
    end    
    %����ȫ�����Ž�
    if fym > min(fxm)           %fym�ĳ�ʼֵ=inf
        [fym, nmax] = min(fxm);   % ����Ⱥ����ʷ�����Ӧ��
        ym = xm(nmax, :);      % ����Ⱥ����ʷ���λ��
    end
    record(iter) = fym;%��Сֵ��¼
    iter = iter+1;
end

%% ������
figure(1)
plot(record);
xlabel('Number of iterations') ;
ylabel('Fitness value') ;
% title('��Ӧ�Ƚ�������');

disp(['��С��ʱ��',num2str(fym)]);
disp(['�������Ƚ����',num2str(ym)]);
end
