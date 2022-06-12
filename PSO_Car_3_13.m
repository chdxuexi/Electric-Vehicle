function [record,fym]=PSO_Car_3_13()
clc;
clear;
close all;
%% 初始化种群的参数设置
tic;%计时用
carNum = 1200;
d = carNum;                      % 100辆车，100维
pileNum = 222;                   % 222个充电桩
N = 150;                         % 初始种群个数（每个粒子都是一个调度方案）
ger = 300;                       % 最大迭代次数
[EV_arrive_time,EV_charging_time] = DataSource;   %传过来车辆、充电桩的时间数据
Wmax=0.9;                        %惯性权重最大值
Wmin=0.4;                        %惯性权重最小值
xlimit = [1, pileNum];           % 设置位置参数限制(矩阵的形式可以多维)
xlimit = repmat(xlimit,d,1);
vlimit = [-10, 10];              % 设置速度限制
vlimit = repmat(vlimit,d,1);
c_2 = 1.49;                      % 自我学习因子
c_3 = 1.49;                      % 群体学习因子
for i = 1:d
    x = xlimit(i, 1) + (xlimit(i, 2) - xlimit(i, 1)) * rand(N, d);    %初始种群的位置
end
x = round(x);                    % 对初始种群的位置进行取整，四舍五入
v = rand(N, d);                  % 初始种群的速度
xm = x;                          % 每个个体的历史最佳位置
ym = zeros(1, d);                % 种群的历史最佳位置
fxm = zeros(N, 1);               % 每个个体的历史最佳适应度
fym = inf;                       % 种群历史最佳适应度
%% 粒子群工作
iter = 1;
times = 1;
record = zeros(ger, 1);          
record_syd = 0;
while iter <= ger 
    c_1 = Wmax-(Wmax-Wmin)*(iter/ger)^2;
    syd = [];     %每一个粒子的适应度值
    for i=1:N     %N为种群个数
        DX=x(i,:);
        T=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);    %调用函数，计算适应度值
        syd=[syd,T];
    end
    fx = syd;                       % 当前代数，N个粒子的适应度值   
    %当iter=1的时候，要将这个fxm和xm初始化，直接赋值给它们
    if iter ==1
        for i = 1:N                 % N = 150 种群个数
            if fxm(i) < fx(i)
                fxm(i) = fx(i);     % 更新个体历史最佳适应度
                xm(i,:) = x(i,:);   % 更新个体历史最佳位置
            end
        end
        %当iter>1的时候，则进行当前fx和fxm、xm进行比较交换的操作
    else
        for i = 1:N                % N = 150 种群个数
            if fxm(i) > fx(i)
                fxm(i) = fx(i);    % 更新个体历史最佳适应度
                xm(i,:) = x(i,:);  % 更新个体历史最佳位置
            end
        end
    end 
    if fym > min(fxm)              %fym的初始值=inf
        [fym, nmax] = min(fxm);    % 更新群体历史最佳适应度
        ym = xm(nmax, :);          % 更新群体历史最佳位置
    end   
    v =  c_1*v  + c_2 * rand *(xm - x) + c_3 * rand *(repmat(ym, N, 1) - x);% 速度更新   
    % 边界速度处理
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
    % 位置更新
    x = x + v;
    % 边界位置处理
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
    record(iter) = fym;            %最小值记录
    iter = iter+1;
    times=times+1;
end
toc;

%% 结果输出
% figure(1)
% plot(record,'Linewidth',1.5);
% xlabel('Number of iterations','Linewidth',1.5) ;
% ylabel('Fitness value','Linewidth',1.5) ;
% % title('适应度进化曲线');
% 
% disp(['最小耗时：',num2str(fym)]);
% disp(['车辆调度结果：',num2str(ym)]);
end
