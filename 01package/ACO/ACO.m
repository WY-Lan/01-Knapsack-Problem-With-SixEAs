function [his_best, best_pop] = ACO(NC_max,Ant_Quantity,lb, ub, dim, Values, Weights,maxw)

a = 1;
b = 0.1;
p = 0.8;
r = 0.8;
D = dim;

t = ones(1, Ant_Quantity);
dt = zeros(1, Ant_Quantity);

Ant_Position = zeros(D, Ant_Quantity);

%% 初始化
temp_Ant_Position = lb + rand(D, Ant_Quantity) * (ub - lb);
best_pop = 0;
his_best = zeros(NC_max,1);
%% 
for NC = 1:NC_max
    t = p * t + dt;
    Ant_Position = temp_Ant_Position;
    for i = 1:Ant_Quantity
        ETA_zero_num = 0;
        for ii = 1:Ant_Quantity
            Ant_Position(:,i) = constraint_population(Ant_Position(:,i),Values,Weights,maxw);
            Ant_Position(:,ii) = constraint_population(Ant_Position(:,ii),Values,Weights,maxw);
            temp = func(Ant_Position(:,i),Values) - func(Ant_Position(:,ii),Values);    %函数值差值赋给临时变量temp，避免重复计算
            if temp > 0
                Ant_ETA(i,ii) = temp;         %Ant_ETA是一个Ant_Quantity阶的方阵，其i行ii列表示从i到ii的启发量
            else
                Ant_ETA(i,ii) = 0;            %如果函数值不下降，启发量为0（会导致后面该方向概率为0）
                ETA_zero_num = ETA_zero_num + 1;    %统计启发量为0的个数    
            end
        end
        
        if ETA_zero_num == Ant_Quantity        %如果启发量全为0，则说明此蚂蚁是本轮蚂蚁中函数值最小的最优蚂蚁
            next = i;                         %它的下一个位置是他自身
        else
            %计算蚂蚁i向各点ii的移动概率
            sum_p = 0;
            for ii = 1:Ant_Quantity   %遍历每一只蚂蚁计算，从而获得移动概率分母上的求和
                sum_p = sum_p + t(ii)^a * Ant_ETA(i,ii)^b;
            end
            for ii = 1:Ant_Quantity   %遍历每一只蚂蚁计算，从而获得移动概率
                Ant_Possibility(i,ii) =  t(ii)^a * Ant_ETA(i,ii)^b / sum_p;    %Ant_Possibility是一个Ant_Quantity阶的方阵，其i行ii列表示从i到ii的概率
            end

            k = 0;
            for ii = 1:Ant_Quantity   %遍历每一只蚂蚁ii，找出蚂蚁i到ii概率不为0的，单独存储其编号和概率值
                if Ant_Possibility(i,ii) ~=0
                    k = k + 1;
                    K(k) = ii;                                   %存储概率不为零的蚂蚁编号
                    K_Possibility(k) = Ant_Possibility(i,ii);       %存储不为0的概率值
                end
            end
            K_Possibility = cumsum(K_Possibility);
            random_p = rand;
            next = 0;
            for ii = 1:k
                if random_p < K_Possibility(ii)
                    next = K(ii);
                end
                if next ~= 0
                    break
                end
            end
        end
        
        %蚂蚁移动，向next的某一邻域内移动
        Ant_Position(:,i) = constraint_population(Ant_Position(:,i),Values,Weights,maxw);
        Ant_Position(:,ii) = constraint_population(Ant_Position(:,ii),Values,Weights,maxw);
        dt(i) =  func(Ant_Position(:,i), Values) - func(Ant_Position(:,next), Values);        %留下本蚂蚁的信息量增量
        temp_Ant_Position(:,i) = Ant_Position(:,next) + (-1 + 2 * rand(D,1)) * r^NC;    %移动到next邻域中的某点
    end
    %计算移动之后所有蚂蚁中函数值最小的，以及最小函数值
    for i = 1:Ant_Quantity
        temp_Ant_Position(:,i) = constraint_population(temp_Ant_Position(:,i),Values,Weights,maxw);
        FUNC(i) = func(temp_Ant_Position(:,i), Values);
    end
    [FUNC_min(NC), FUNC_min_n] = min(FUNC);
    Position_min(:,NC) = temp_Ant_Position(:, FUNC_min_n);
    if NC == 1
        his_best(NC) = FUNC_min(NC);
        best_pop = temp_Ant_Position(:, FUNC_min_n);
    else
        if his_best(NC-1) > FUNC_min(NC)
            his_best(NC) = FUNC_min(NC);
            best_pop = temp_Ant_Position(:, FUNC_min_n);
        else
            his_best(NC) = his_best(NC-1);
        end
    end
end
end



function population = constraint_population(population, Values, Weights, maxw)
    for i = 1:size(population, 2)
        tmp_value = 0;
        tmp_weight = 0;
        for j = 1:size(population,1)
            if population(j, i) > 0.5
                tmp_value = tmp_value+Values(j);
                tmp_weight = tmp_weight+Weights(j);
            end
        end

        while tmp_weight > maxw
            indices_of_ones = find(population(:,i) > 0.5);
            random_index = indices_of_ones(randi(size(indices_of_ones,2)));
            population(random_index, i) = rand * 0.5;
            tmp_value = tmp_value - Values(random_index);
            tmp_weight = tmp_weight - Weights(random_index);
        end
    end

end

function fitness = func(X, Values)
% 求解最小值，故为负数
    fitness = 0;
    for i = 1:size(X,1)
        if X(i) > 0.5
            fitness = fitness - Values(i);
        end
    end
end