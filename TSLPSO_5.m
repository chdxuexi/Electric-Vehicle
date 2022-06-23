function [record,fym]=TSLPSO_5()
clc;
clear;
close all;
%% Initialize the parameter settings of the population
carNum = 600;
d = carNum;                      % 600 vehicles, 600 dimensions
pileNum = 222;                   % 222 charging piles
N = 150;                         % initial population
ger = 300;                       % iterations
Wmax=0.9;                        
Wmin=0.4;                        
c_2 = 1.49;                       % cognitive coefficient
c_3 = 1.49;                       % social coefficient
[EV_arrive_time,EV_charging_time] = DataSource; 
% initial velocity and location
xmin = 1;                  
xmax = pileNum;
vmin = -10;                
vmax = 10;
for i0=1:N
    for j0 = 1:d
        x(i0,j0) = xmin + (xmax - xmin) * rand; 
        v(i0,j0) = xmin + (xmax - xmin) * rand; 
    end
end
x = round(x);                    % round the position of the initial population
syd = [];                       
for i=1:N                        
    DX=x(i,:);
    T=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time);    
    syd=[syd,T];
end
fx = syd;                        % fitness value of N particles
xm = x;                          % pbest
ym = zeros(1, d);                % gbest
fxm = zeros(N, 1);               % f(pbest)
fym = inf;                       % f(gbest)
for i = 1:N                      
    fxm(i) = fx(i);     % update f(pbest)
    xm(i,:) = x(i,:);   % update pbest
end
if fym > min(fxm)            
    [fym, nmax] = min(fxm);   % update f(gbest)
    ym = xm(nmax, :);         % update gbest
end

%% enter the loop
iter = 1;
xdl=xm; 
record = zeros(ger, 1); 
while iter <= ger  
    c_1 = Wmax-(Wmax-Wmin)*(iter/ger)^2;
    for i = 1:N
    v(i,:) =  c_1*v(i,:)+c_2 * rand *(xdl(i,:) - x(i,:)) + c_3 * rand *(ym - x(i,:));     
    % Boundary Velocity Handling
        for j=1:d
                if  v(i,j)>vmax
                    v(i,j)=vmax;
                end
                if  v(i,j) < vmin
                    v(i,j)=vmin;
                end
        end
        v = round(v);
        % update location 
        x(i,:) = x(i,:) + v(i,:);
        % Boundary location handling
        for j=1:d
                if  x(i,j)>xmax
                    x(i,j)=xmax;
                end
                if  x(i,j) < xmin
                    x(i,j)=xmin;
                end
        end       
        %Calculate fitness
        DX=x(i,:);
        fx(i)=Car_Fitness(DX,N,EV_arrive_time,EV_charging_time); 
        if fx(i)<fxm(i)
            fxm(i) = fx(i);     
            xm(i,:) = x(i,:);   
            flag(i)=1;
        else
            flag(i)=0;
        end
        if flag(i)==1
            % learning particle
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
            end
        end
    end    
    % Update the global optimal solution
    if fym > min(fxm)           
        [fym, nmax] = min(fxm);   
        ym = xm(nmax, :);      
    end
    record(iter) = fym;  
    iter = iter+1;
end

%% conclusion
figure(1)
plot(record);
xlabel('Number of iterations') ;
ylabel('Fitness value') ;
% title('适应度进化曲线');

disp(['最小耗时：',num2str(fym)]);
disp(['车辆调度结果：',num2str(ym)]);
end
