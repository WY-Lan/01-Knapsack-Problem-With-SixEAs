clear,clc;
fidin=fopen('..\Data\large_scale\knapPI_1_100_1000_1');    % 打开test2.txt文件             
tline = split(fgetl(fidin));
maxn = str2double(tline(1));
maxw = str2double(tline(2));
Values = zeros(maxn, 1);
Weights = zeros(maxn,1);
for i = 1:maxn
    tline = split(fgetl(fidin));
    Values(i) = str2double(tline(1));
    Weights(i) = str2double(tline(2));
end


%% 遗传算法，二进制编码，区分其它几种算法
N = 50;
Tmax = 1000;
crossRate = 0.8;
mutateRate = 0.05;
his_fitness = zeros(Tmax, 1);

%% 初始化个体
population = round(rand(N, maxn));
population = constraint_population(population,Values,Weights,maxw);
Fitness = zeros(N,1);
%% 迭代
for t = 1:Tmax
    population = constraint_population(population,Values,Weights,maxw);
    for i = 1:size(population, 1)
        Fitness(i) = cal_fitness(population(i,:),Values, Weights, maxw);
    end
    [current_fitness, current_pop] = max(Fitness);
    if t == 1
        his_fitness(t) = current_fitness;
        best_pop = population(current_pop,:);
    else
        if his_fitness(t-1) < current_fitness
            his_fitness(t) = current_fitness;
            best_pop = population(current_pop,:);
        else
            his_fitness(t) = his_fitness(t-1);
        end
    end

    % 选择 轮盘
    Prefix_Fit = cumsum(Fitness);
    Prefix_Fit = Prefix_Fit ./ sum(Fitness,1);
    for i = 1:size(population, 1)
        rd = rand;
        for j = 1:size(Prefix_Fit, 1)
            if rd < Prefix_Fit(j)
                population(i,:) = population(j,:);
                break;
            end
        end
    end

    % 交叉操作
    for i=1:2:size(population, 1)
        rd = rand;
        if rd < crossRate
            % 产生一个随机位置
            cut = floor( rand(1,1) * size(population,2)) + 1;
            for j = cut:size(population, 2)
                tmp = population(i,j);
                population(i,j) = population(i+1,j);
                population(i+1,j) = tmp;
            end
        end
    end

    % 变异操作
    for i=1:size(population, 1)
        rd = rand;
        if rd < mutateRate
            cut = floor( rand(1,1) * size(population,2)) + 1;
            if population(i, cut) == 1
                population(i, cut) = 0;
            else
                population(i, cut) = 1;
            end
        end
    end
end
plot(his_fitness);
best_weight = 0;
best_fitness = 0;
for i = 1:size(best_pop,2)
    if best_pop(i) == 1
        best_fitness = best_fitness + Values(i);
        best_weight = best_weight + Weights(i);
    end
end
best_fitness
best_weight


function population = constraint_population(population, Values, Weights, maxw)
    for i = 1:size(population,1)
        tmp_value = 0;
        tmp_weight = 0;
        for j = 1:size(population,2)
            if population(i,j) == 1
                tmp_value = tmp_value + Values(j);
                tmp_weight = tmp_weight + Weights(j);
            end
        end

        while tmp_weight > maxw
            indices_of_ones = find(population(i,:) == 1);
            random_index = indices_of_ones(randi(size(indices_of_ones,2)));
            population(i,random_index) = 0;
            tmp_value = tmp_value - Values(random_index);
            tmp_weight = tmp_weight - Weights(random_index);
        end
    end
end

function fitness = cal_fitness(X, Values, Weights, maxw)
% 寻找最大适应度
    current_weight = 0;
    fitness = 0;
    for i = 1:size(X,2)
        if X(i) == 1
            fitness = fitness + Values(i);
            current_weight = current_weight + Weights(i);
        end
    end
end