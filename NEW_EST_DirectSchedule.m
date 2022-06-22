clc;
clear;
%% Sort by the first column, and select the vehicle with the smallest arrival time to schedule charging
new_lie = sortArrayEachRow(carNum,new_lie);                        
new_lie = sortrows(new_lie);                                        
priority = new_lie(:,pileNum+1);                                   
%% vehicle2station[] stores the relationship between the vehicle and the pile
vehicle2pile= zeros(1,pileNum);                                     
vehicle2pile = [vehicle2pile;new_hangCharing(carNum+1:carNum+3,1:pileNum)];
SumTimeOfPlies = zeros(1,pileNum);                                 
for i=1:carNum
    TheCar = priority(i);                                         
    tempSingleCar = StoreEachsingleCar(:,:,TheCar);                 
        %先判断在j桩是否有车排队
        tempSingleCar(4,1);%这是该车TheCar在这个桩
        if SumTimeOfPlies(1,tempSingleCar(4,1)) == 0
            %如果无，则安排车辆到该桩，并修改时间累加矩阵
            vehicle2pile(1,tempSingleCar(4,1)) = tempSingleCar(1,pileNum+1);%这里是车辆编号=tempSingleCar(2,1)
            SumTimeOfPlies(1,tempSingleCar(4,1)) = SumTimeOfPlies(1,tempSingleCar(4,1)) + tempSingleCar(1,1) + tempSingleCar(5,1);%这里是该车在该桩的（到+冲）时间
        else
            %如果有，则遍历到所有桩的代价，找到最小的一个安排车辆到该桩，并修改时间累加矩阵
            tempCost = zeros(1,pileNum);%获得2号车到所有充电桩的代价
            for k=1:pileNum
                if SumTimeOfPlies(1,k) ~= 0 %如果该桩有车在充电
                    tempCost(1,k) = SumTimeOfPlies(1,k);%更新在该桩充电的代价，就为原有的
                else%如果无车在该桩充电
                    tempOfLie = find(tempSingleCar(4,:) == k);
                    tempCost(1,k) = tempSingleCar(1,tempOfLie);%更新在该桩充电的代价，即为到达该桩的时间（不包括充电时间）
                end
            end 
            %从代价矩阵中找到代价最小的桩 ，将该车指定到这个桩   
            [C,index] = min(tempCost); %[C,I] = min(A)找到A中那些最小值C的索引位置，将他们放在向量index中返回。如果这里有多个相同最小值时，返回的将是第一个的索引。
            tempPile = index(1);%tempPile为该车确定在该桩充电，桩的编号
            if vehicle2pile(1,tempPile) == 0 
                %如果无
                vehicle2pile(1,tempPile) =  tempSingleCar(1,pileNum+1);%tempSingleCar(1,pileNum+1)为该车的编号；
            else %如果有，则将该车存储在该桩的最后一列
                [mm,nn] = size(vehicle2pile(:,tempPile));
                if vehicle2pile(mm,tempPile) == 0
                    vehicle2pile(mm,tempPile) =  tempSingleCar(1,pileNum+1);%tempSingleCar(1,pileNum+1)为该车的编号；
                else
                    vehicle2pile(mm+1,tempPile) =  tempSingleCar(1,pileNum+1);%tempSingleCar(1,pileNum+1)为该车的编号；
                end
            end
            SumTimeOfPlies(1,tempPile) = SumTimeOfPlies(1,tempPile) + tempSingleCar(1,tempPile) + tempSingleCar(5,tempPile);%更新代价矩阵为原有矩阵接收该车充电的桩+该车到达+该车充电
        end
end
%% 可视化每个充电桩有多少车在这里充电
for i=1:pileNum
    if vehicle2pile(1,i) == 0%如果这个数为零，说明该桩无车充电,所以跳过该桩
        continue;
    end
    vehicle2pile(1,i);
    carInitPosition(vehicle2pile(1,i));
    fprintf("第 %3d 辆在节点 %3d，到充电站 %3d 的 %3d 号桩充电\n",vehicle2pile(1,i),carInitPosition(vehicle2pile(1,i)),vehicle2pile(2,i),vehicle2pile(4,i));
end
%% 由vehicle2pile，将哪些车在哪些充电站充电的关系保存到vehicle2station
vehicle2station = zeros(1,carNum);
for i=1:pileNum
    fprintf('i=%3d',i)
    if vehicle2pile(1,i) == 0%如果这个数为零，说明该桩无车充电,所以跳过该桩
        continue;
    end
    vehicle2pile(1,i);
    vehicle2pile(2,i);
    vehicle2station(vehicle2pile(1,i)) = vehicle2pile(2,i);
    [am,an] = size(vehicle2pile(:,i));
    for kk=5:am
        if vehicle2pile(kk,i) == 0
            continue;
        else
            vehicle2station(vehicle2pile(kk,i)) = vehicle2pile(2,i);
        end
    end
end

%% 可视化每个充电站有多少车在这里充电
pileStationNum = [9 13 20 30 33 35 38 42 47 58 60 62 63 66 75 80 84];%充电桩的编号情况：
pileStationNumValue = [2 9 22 21 2 13 6 8 5 9 3 51 7 12 6 28 18];%充电桩的编号情况：
stationNum = 17;
for i=1:stationNum
    pile = pileStationNum(i);%充电站序号
    pileValue = pileStationNumValue(i);%该充电站充电桩数目
    sum1 = length(find(vehicle2station==pile));
    fprintf('\n 充电站 %3d 有%3d个桩；  %3d 辆车在这充电',pile,pileValue,sum1);
end
%% 计算所有车的总充电时间           
za = unique(vehicle2station);%za保存的是充电站的序号，且充电站有车辆进行充电
station_pile = [0,0,0,0,0,0,0,0,2,0,0,0,9,0,0,0,0,0,0,22,0,0,0,0,0,0,0,0,0,21,0,0,2,0,13,0,0,6,0,0,0,8,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,9,0,3,0,51,7,0,0,12,0,0,0,0,0,0,0,0,6,0,0,0,0,28,0,0,0,18,0,0,0,0,0,0];
            %充电站和充电桩的对应关系/上面为充电站：下面对应的充电桩的个数
T = 0; %总的消耗的所有时间 
EST_car_wait = zeros(1,carNum); 
EST_car_arrive = zeros(1,carNum); 
EST_car_charge = zeros(1,carNum); 
for i=1:length(za)
    zai = za(i);
    [zmmm,zmvv] = find(vehicle2station == zai);
    if length(zmvv) <= station_pile(zai)
        for j=1:length(zmvv) 
            T = T + EV_arrive_time(zmvv(j),zai) + EV_charging_time(zmvv(j),zai); 
            EST_car_arrive(1,zmvv(j)) = EV_arrive_time(zmvv(j),zai); 
            EST_car_charge(1,zmvv(j)) = EV_charging_time(zmvv(j),zai);
        end
    else
        %排队
        medium_car_arrive = [];
        medium_car_charging = [];
        for j=1:length(zmvv)
             fprintf("\n第 %3d 辆车在充电站 %3d 充电 该车辆在充电站编号为：%3d \n" ,zmvv(j),zai,carInitPosition(zmvv(j)))
            medium_car_arrive = [medium_car_arrive,EV_arrive_time(zmvv(j),zai)];
            medium_car_charging = [medium_car_charging,EV_charging_time(zmvv(j),zai)];
        end
        for j=1:length(zmvv)
            medium_car_arrive(2,j) = zmvv(j);
            EST_car_arrive(1,zmvv(j)) = medium_car_arrive(1,j);%保存该车的到达时间
            medium_car_charging(2,j) = zmvv(j);
            EST_car_charge(1,zmvv(j)) = medium_car_charging(1,j);%保存该车的充电时间
        end
        %到达时间进行排序
        for ii=1:length(medium_car_arrive)
            for jj=1:length(medium_car_arrive)-ii
                if(medium_car_arrive(jj)>medium_car_arrive(jj+1))
                    [medium_car_arrive(1,jj),medium_car_arrive(1,jj+1)] = swap(medium_car_arrive(1,jj),medium_car_arrive(1,jj+1));
                    [medium_car_arrive(2,jj),medium_car_arrive(2,jj+1)] = swap(medium_car_arrive(2,jj),medium_car_arrive(2,jj+1));
                    [medium_car_charging(1,jj),medium_car_charging(1,jj+1)] = swap(medium_car_charging(1,jj),medium_car_charging(1,jj+1)); 
                    [medium_car_charging(2,jj),medium_car_charging(2,jj+1)] = swap(medium_car_charging(2,jj),medium_car_charging(2,jj+1)); 
                end
            end
        end
        [T_temp,wait] = getChargeTime(medium_car_arrive,medium_car_charging,station_pile(zai));
        T = T +  T_temp;
        for k=1:length(zmvv)
            EST_car_wait(wait(2,k)) = wait(1,k);
        end
    end
end
fprintf('\nNEW_EST最早开始充电策略调度总的时间是: %d 小时\n\n',T)
fprintf('\n总的时间是: %d 小时\n\n',sum(EST_car_arrive+EST_car_charge+EST_car_wait));
