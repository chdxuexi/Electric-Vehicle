function [record,fym]=PSO_Car_3_13()
clc;
clear;
close all;
%% Initialize the parameter settings of the population
tic;%计时用
carNum = 1200;
d = carNum;                      % 100 vehicles, 100 dimensions
pileNum = 222;                   % 222 charging piles
N = 150;                         % initial population
ger = 300;                       % iterations
[EV_arrive_time,EV_charging_time] = DataSource;   % Data of vehicles and charging piles
Wmax=0.9;                        
Wmin=0.4;                       
xlimit = [1, pileNum];           % Set positional parameter limits
xlimit = repmat(xlimit,d,1);
vlimit = [-10, 10];              % set velocity limit
vlimit = repmat(vlimit,d,1);
c_2 = 1.49;                      % cognitive coefficient
c_3 = 1.49;                      % social coefficient
for i = 1:d
    x = xlimit(i, 1) + (xlimit(i, 2) - xlimit(i, 1)) * rand(N, d);    % The location of the initial population
end
x = round(x);                    % round the position of the initial population
v = rand(N, d);                  % The velocity of the initial population
xm = x;                          % pbest
ym = zeros(1, d);                % gbest
fxm = zeros(N, 1);               % f(pbest)
fym = inf;                       % f(gbest)
%% pevcs
iter = 1;
times = 1;
record = zeros(ger, 1);          
record_syd = 0;
while iter <= ger 
    c_1 = Wmax-(Wmax-Wmin)*(iter/ger)^2;
    syd = [];     % The fitness value of each particle
    for i=1:N     
        DX=x(i,:);
        T=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);    % Call the function to calculate the fitness value
        syd=[syd,T];
    end
    fx = syd;                       % fx is the current iteration 
    % When iter=1, initialize fxm and xm
    if iter ==1
        for i = 1:N                
            if fxm(i) < fx(i)
                fxm(i) = fx(i);     % update f(pbest)
                xm(i,:) = x(i,:);   % update pbest
            end
        end
        % When iter>1, compare fx with fxm and xm
    else
        for i = 1:N                
            if fxm(i) > fx(i)
                fxm(i) = fx(i);    % update f(pbest)
                xm(i,:) = x(i,:);  % update pbest
            end
        end
    end 
    if fym > min(fxm)              % Initial fym = inf
        [fym, nmax] = min(fxm);    % update f(gbest)
        ym = xm(nmax, :);          % update gbest
    end   
    v =  c_1*v  + c_2 * rand *(xm - x) + c_3 * rand *(repmat(ym, N, 1) - x); % update velocity  
    % Boundary Velocity 
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
    % update location
    x = x + v;
    % Boundary location
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
    record(iter) = fym;            
    iter = iter+1;
    times=times+1;
end
toc;

%% conclusion
% figure(1)
% plot(record,'Linewidth',1.5);
% xlabel('Number of iterations','Linewidth',1.5) ;
% ylabel('Fitness value','Linewidth',1.5) ;
% % title('适应度进化曲线');
% 
% disp(['最小耗时：',num2str(fym)]);
% disp(['车辆调度结果：',num2str(ym)]);
end
