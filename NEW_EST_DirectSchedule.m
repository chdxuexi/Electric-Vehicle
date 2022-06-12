clc;
clear;
%% ����һ������ѡ�񵽴�ʱ����С�����������е��ȳ��
new_lie = sortArrayEachRow(carNum,new_lie);                         %����sortArrayEachRow�����������ȶ�ÿһ�зֱ�������򣬲���ÿһ�к����������Ӧ�ĳ������1��2,3....100��
new_lie = sortrows(new_lie);                                        %��new_lie���յ�һ�еĴ�С����ÿ�н�������
priority = new_lie(:,pileNum+1);                                    %�����͵õ������г�������˳���Ⱥ������
%% ����ÿ�������ĸ����׮���г�磬�����ϵ���浽 vehicle2station[]
vehicle2pile= zeros(1,pileNum);                                     %�����ǳ��׮��ţ�ֵ�ǳ�����ţ����·���㷨�����ﲻһ��
vehicle2pile = [vehicle2pile;new_hangCharing(carNum+1:carNum+3,1:pileNum)];%�����Ǹ�vehicle2pile�����м����˱��1,2,3��....pileNum
SumTimeOfPlies = zeros(1,pileNum);                                  %�������׮��ʱ���ۼӾ���
for i=1:carNum
    TheCar = priority(i);                                           %����ѭ�������ŵĳ������ = TheCar
    tempSingleCar = StoreEachsingleCar(:,:,TheCar);                 %����ѭ��,�ó���������Ϣȡ����tempSingleCar
        %���ж���j׮�Ƿ��г��Ŷ�
        tempSingleCar(4,1);%���Ǹó�TheCar�����׮
        if SumTimeOfPlies(1,tempSingleCar(4,1)) == 0
            %����ޣ����ų�������׮�����޸�ʱ���ۼӾ���
            vehicle2pile(1,tempSingleCar(4,1)) = tempSingleCar(1,pileNum+1);%�����ǳ������=tempSingleCar(2,1)
            SumTimeOfPlies(1,tempSingleCar(4,1)) = SumTimeOfPlies(1,tempSingleCar(4,1)) + tempSingleCar(1,1) + tempSingleCar(5,1);%�����Ǹó��ڸ�׮�ģ���+�壩ʱ��
        else
            %����У������������׮�Ĵ��ۣ��ҵ���С��һ�����ų�������׮�����޸�ʱ���ۼӾ���
            tempCost = zeros(1,pileNum);%���2�ų������г��׮�Ĵ���
            for k=1:pileNum
                if SumTimeOfPlies(1,k) ~= 0 %�����׮�г��ڳ��
                    tempCost(1,k) = SumTimeOfPlies(1,k);%�����ڸ�׮���Ĵ��ۣ���Ϊԭ�е�
                else%����޳��ڸ�׮���
                    tempOfLie = find(tempSingleCar(4,:) == k);
                    tempCost(1,k) = tempSingleCar(1,tempOfLie);%�����ڸ�׮���Ĵ��ۣ���Ϊ�����׮��ʱ�䣨���������ʱ�䣩
                end
            end 
            %�Ӵ��۾������ҵ�������С��׮ �����ó�ָ�������׮   
            [C,index] = min(tempCost); %[C,I] = min(A)�ҵ�A����Щ��СֵC������λ�ã������Ƿ�������index�з��ء���������ж����ͬ��Сֵʱ�����صĽ��ǵ�һ����������
            tempPile = index(1);%tempPileΪ�ó�ȷ���ڸ�׮��磬׮�ı��
            if vehicle2pile(1,tempPile) == 0 
                %�����
                vehicle2pile(1,tempPile) =  tempSingleCar(1,pileNum+1);%tempSingleCar(1,pileNum+1)Ϊ�ó��ı�ţ�
            else %����У��򽫸ó��洢�ڸ�׮�����һ��
                [mm,nn] = size(vehicle2pile(:,tempPile));
                if vehicle2pile(mm,tempPile) == 0
                    vehicle2pile(mm,tempPile) =  tempSingleCar(1,pileNum+1);%tempSingleCar(1,pileNum+1)Ϊ�ó��ı�ţ�
                else
                    vehicle2pile(mm+1,tempPile) =  tempSingleCar(1,pileNum+1);%tempSingleCar(1,pileNum+1)Ϊ�ó��ı�ţ�
                end
            end
            SumTimeOfPlies(1,tempPile) = SumTimeOfPlies(1,tempPile) + tempSingleCar(1,tempPile) + tempSingleCar(5,tempPile);%���´��۾���Ϊԭ�о�����ոó�����׮+�ó�����+�ó����
        end
end
%% ���ӻ�ÿ�����׮�ж��ٳ���������
for i=1:pileNum
    if vehicle2pile(1,i) == 0%��������Ϊ�㣬˵����׮�޳����,����������׮
        continue;
    end
    vehicle2pile(1,i);
    carInitPosition(vehicle2pile(1,i));
    fprintf("�� %3d ���ڽڵ� %3d�������վ %3d �� %3d ��׮���\n",vehicle2pile(1,i),carInitPosition(vehicle2pile(1,i)),vehicle2pile(2,i),vehicle2pile(4,i));
end
%% ��vehicle2pile������Щ������Щ���վ���Ĺ�ϵ���浽vehicle2station
vehicle2station = zeros(1,carNum);
for i=1:pileNum
    fprintf('i=%3d',i)
    if vehicle2pile(1,i) == 0%��������Ϊ�㣬˵����׮�޳����,����������׮
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

%% ���ӻ�ÿ�����վ�ж��ٳ���������
pileStationNum = [9 13 20 30 33 35 38 42 47 58 60 62 63 66 75 80 84];%���׮�ı�������
pileStationNumValue = [2 9 22 21 2 13 6 8 5 9 3 51 7 12 6 28 18];%���׮�ı�������
stationNum = 17;
for i=1:stationNum
    pile = pileStationNum(i);%���վ���
    pileValue = pileStationNumValue(i);%�ó��վ���׮��Ŀ
    sum1 = length(find(vehicle2station==pile));
    fprintf('\n ���վ %3d ��%3d��׮��  %3d ����������',pile,pileValue,sum1);
end
%% �������г����ܳ��ʱ��           
za = unique(vehicle2station);%za������ǳ��վ����ţ��ҳ��վ�г������г��
station_pile = [0,0,0,0,0,0,0,0,2,0,0,0,9,0,0,0,0,0,0,22,0,0,0,0,0,0,0,0,0,21,0,0,2,0,13,0,0,6,0,0,0,8,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,9,0,3,0,51,7,0,0,12,0,0,0,0,0,0,0,0,6,0,0,0,0,28,0,0,0,18,0,0,0,0,0,0];
            %���վ�ͳ��׮�Ķ�Ӧ��ϵ/����Ϊ���վ�������Ӧ�ĳ��׮�ĸ���
T = 0; %�ܵ����ĵ�����ʱ�� 
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
        %�Ŷ�
        medium_car_arrive = [];
        medium_car_charging = [];
        for j=1:length(zmvv)
             fprintf("\n�� %3d �����ڳ��վ %3d ��� �ó����ڳ��վ���Ϊ��%3d \n" ,zmvv(j),zai,carInitPosition(zmvv(j)))
            medium_car_arrive = [medium_car_arrive,EV_arrive_time(zmvv(j),zai)];
            medium_car_charging = [medium_car_charging,EV_charging_time(zmvv(j),zai)];
        end
        for j=1:length(zmvv)
            medium_car_arrive(2,j) = zmvv(j);
            EST_car_arrive(1,zmvv(j)) = medium_car_arrive(1,j);%����ó��ĵ���ʱ��
            medium_car_charging(2,j) = zmvv(j);
            EST_car_charge(1,zmvv(j)) = medium_car_charging(1,j);%����ó��ĳ��ʱ��
        end
        %����ʱ���������
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
fprintf('\nNEW_EST���翪ʼ�����Ե����ܵ�ʱ����: %d Сʱ\n\n',T)
fprintf('\n�ܵ�ʱ����: %d Сʱ\n\n',sum(EST_car_arrive+EST_car_charge+EST_car_wait));