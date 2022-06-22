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
        %First determine whether there is a vehicle queue at the CPj 
        tempSingleCar(4,1);    
        if SumTimeOfPlies(1,tempSingleCar(4,1)) == 0
            %If none, arrange the vehicle to the pile
            vehicle2pile(1,tempSingleCar(4,1)) = tempSingleCar(1,pileNum+1);  %tempSingleCar(2,1) is the EVi’id
            SumTimeOfPlies(1,tempSingleCar(4,1)) = SumTimeOfPlies(1,tempSingleCar(4,1)) + tempSingleCar(1,1) + tempSingleCar(5,1); %charging time+arrive time
        else
            %If there is, traverse to the cost of all piles, find the shortest time, and arrange the vehicle to the pile
            tempCost = zeros(1,pileNum);  
            for k=1:pileNum
                if SumTimeOfPlies(1,k) ~= 0   %If the pile has a vehicle charging
                    tempCost(1,k) = SumTimeOfPlies(1,k);   
                else   %If no vehicle is charging at the pile
                    tempOfLie = find(tempSingleCar(4,:) == k);
                    tempCost(1,k) = tempSingleCar(1,tempOfLie); 
                end
            end 
            %Find the pile with the shortest time cost and assign the vehicle to this pile   
            [C,index] = min(tempCost);   
            tempPile = index(1);   %tempPile is the CPj’s id
            if vehicle2pile(1,tempPile) == 0 
                % if no
                vehicle2pile(1,tempPile) =  tempSingleCar(1,pileNum+1); %tempSingleCar(1,pileNum+1) is the EVi’s id；
            else % If there is, store the vehicle in the last column of the pile
                [mm,nn] = size(vehicle2pile(:,tempPile));
                if vehicle2pile(mm,tempPile) == 0
                    vehicle2pile(mm,tempPile) =  tempSingleCar(1,pileNum+1); 
                else
                    vehicle2pile(mm+1,tempPile) =  tempSingleCar(1,pileNum+1);
                end
            end
            SumTimeOfPlies(1,tempPile) = SumTimeOfPlies(1,tempPile) + tempSingleCar(1,tempPile) + tempSingleCar(5,tempPile);
        end
end
%% Visualize how many vehicles are charging here at each charging pile
for i=1:pileNum
    if vehicle2pile(1,i) == 0 
        continue;
    end
    vehicle2pile(1,i);
    carInitPosition(vehicle2pile(1,i));
    fprintf("第 %3d 辆在节点 %3d，到充电站 %3d 的 %3d 号桩充电\n",vehicle2pile(1,i),carInitPosition(vehicle2pile(1,i)),vehicle2pile(2,i),vehicle2pile(4,i));
end
%% vehicle2station[] stores the relationship between the vehicle and the pile
vehicle2station = zeros(1,carNum);
for i=1:pileNum
    fprintf('i=%3d',i)
    if vehicle2pile(1,i) == 0  
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

%% Visualize how many vehicles are charging here at each charging station
pileStationNum = [9 13 20 30 33 35 38 42 47 58 60 62 63 66 75 80 84]; %Number of charging pile nodes in the road network
pileStationNumValue = [2 9 22 21 2 13 6 8 5 9 3 51 7 12 6 28 18]; %the number of each charging pile
stationNum = 17;
for i=1:stationNum
    pile = pileStationNum(i); 
    pileValue = pileStationNumValue(i);
    sum1 = length(find(vehicle2station==pile));
    fprintf('\n 充电站 %3d 有%3d个桩；  %3d 辆车在这充电',pile,pileValue,sum1);
end
%% Calculate the total charging time for all cars           
za = unique(vehicle2station);  
station_pile = [0,0,0,0,0,0,0,0,2,0,0,0,9,0,0,0,0,0,0,22,0,0,0,0,0,0,0,0,0,21,0,0,2,0,13,0,0,6,0,0,0,8,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,9,0,3,0,51,7,0,0,12,0,0,0,0,0,0,0,0,6,0,0,0,0,28,0,0,0,18,0,0,0,0,0,0];
            
T = 0; %total charging time 
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
        %waiting
        medium_car_arrive = [];
        medium_car_charging = [];
        for j=1:length(zmvv)
             fprintf("\n第 %3d 辆车在充电站 %3d 充电 该车辆在充电站编号为：%3d \n" ,zmvv(j),zai,carInitPosition(zmvv(j)))
            medium_car_arrive = [medium_car_arrive,EV_arrive_time(zmvv(j),zai)];
            medium_car_charging = [medium_car_charging,EV_charging_time(zmvv(j),zai)];
        end
        for j=1:length(zmvv)
            medium_car_arrive(2,j) = zmvv(j);
            EST_car_arrive(1,zmvv(j)) = medium_car_arrive(1,j); %arrive time
            medium_car_charging(2,j) = zmvv(j);
            EST_car_charge(1,zmvv(j)) = medium_car_charging(1,j); %charging time
        end
        %sort the arrive time
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
