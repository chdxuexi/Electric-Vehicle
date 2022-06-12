%% 主函数TSLPSO_5  100车改600车
function [record,fym]=TSLPSO_5()
clc;
clear;
close all;
%% 初始化种群的参数设置
carNum = 600;
d = carNum;                      % 600辆车，600维
pileNum = 222;                   % 222个充电桩
N = 150;                         % 初始种群个数（每个粒子都是一个调度方案）
ger = 300;                       % 最大迭代次数
Wmax=0.9;                        %惯性权重最大值
Wmin=0.4;                        %惯性权重最小值
c_2 = 1.49;                       % 自我学习因子
c_3 = 1.49;                       % 群体学习因子
[EV_arrive_time,EV_charging_time] = DataSource; %传过来车辆、充电桩的时间数据
% 初始化速度v x
xmin = 1;                  % 设置位置参数限制(矩阵的形式可以多维)
xmax = pileNum;
vmin = -10;                % 设置速度限制
vmax = 10;
for i0=1:N
    for j0 = 1:d
        x(i0,j0) = xmin + (xmax - xmin) * rand; %初始种群的位置
        v(i0,j0) = xmin + (xmax - xmin) * rand; %初始种群的速度
    end
end
x = round(x);                    % 对初始种群的位置进行取整，四舍五入
syd = [];                        %每一个粒子的适应度值
for i=1:N                        %N为种群个数
    DX=x(i,:);
    T=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);    %调用函数，计算适应度值
    syd=[syd,T];
end
fx = syd;                        % 当前代数，N个粒子的适应度值
xm = x;                          % 每个个体的历史最佳位置
ym = zeros(1, d);                % 种群的历史最佳位置
fxm = zeros(N, 1);               % 每个个体的历史最佳适应度
fym = inf;                       % 种群历史最佳适应度
for i = 1:N                      % N = 150 种群个数
    fxm(i) = fx(i);     % 更新个体历史最佳适应度
    xm(i,:) = x(i,:);   % 更新个体历史最佳位置
end
if fym > min(fxm)             %fym的初始值=inf
    [fym, nmax] = min(fxm);   % 更新群体历史最佳适应度
    ym = xm(nmax, :);         % 更新群体历史最佳位置
end

%% 进入循环
iter = 1;
xdl=xm;%更换局部最优值
record = zeros(ger, 1); 
while iter <= ger % ger迭代次数
    %对种群1
    c_1 = Wmax-(Wmax-Wmin)*(iter/ger)^2;
    for i = 1:N
    v(i,:) =  c_1*v(i,:)+c_2 * rand *(xdl(i,:) - x(i,:)) + c_3 * rand *(ym - x(i,:));     
    % 边界速度处理
        for j=1:d
                if  v(i,j)>vmax
                    v(i,j)=vmax;
                end
                if  v(i,j) < vmin
                    v(i,j)=vmin;
                end
        end
        v = round(v);
        % 位置更新
        x(i,:) = x(i,:) + v(i,:);
        % 边界位置处理
        for j=1:d
                if  x(i,j)>xmax
                    x(i,j)=xmax;
                end
                if  x(i,j) < xmin
                    x(i,j)=xmin;
                end
        end       
        %适应度计算
        DX=x(i,:);
        fx(i)=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);%调用函数，计算适应度值   
        if fx(i)<fxm(i)
            fxm(i) = fx(i);     % 更新个体历史最佳适应度
            xm(i,:) = x(i,:);   % 更新个体历史最佳位置
            flag(i)=1;
        else
            flag(i)=0;
        end
        if flag(i)==1
            %算法2开始；怀疑返回xdl的值
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
            end%算法2结束
        end
    end    
    %更新全局最优解
    if fym > min(fxm)           %fym的初始值=inf
        [fym, nmax] = min(fxm);   % 更新群体历史最佳适应度
        ym = xm(nmax, :);      % 更新群体历史最佳位置
    end
    record(iter) = fym;%最小值记录
    iter = iter+1;
end

%% 结果输出
figure(1)
plot(record);
xlabel('Number of iterations') ;
ylabel('Fitness value') ;
% title('适应度进化曲线');

disp(['最小耗时：',num2str(fym)]);
disp(['车辆调度结果：',num2str(ym)]);
end
