function [globalbest_x,globalbest_faval,his_best] = PSO(MaxNum,particlesize,lb,ub,dim,Values,Weights,maxw)
vmax = 0.5;
%% 种群初始化
X = lb+(ub-lb)*rand(particlesize,dim);     %粒子所在的位置
V = rand(particlesize,dim);         %粒子的飞翔速度

his_best = zeros(MaxNum,1);
fitness=zeros(particlesize,1);

X = constraint_population(X,Values,Weights,maxw);
for i=1:particlesize
    fitness(i)=fobj(X(i,:),Values);
end

personalbest_x=X;
personalbest_faval=fitness;
[globalbest_faval,i]=min(personalbest_faval);
globalbest_x=personalbest_x(i,:);
k=1;
while k<=MaxNum
    
    c1 = 0.5 + exp(-k/MaxNum)*k/(4*sqrt(2*pi)*MaxNum);
    c2 = 2.5 - exp(-k/MaxNum)*k/(8*sqrt(2*pi)*MaxNum);
    w = 0.9 - 0.005*k;
    for i=1:particlesize
        X(i,:) = constraint_population(X(i,:),Values,Weights,maxw);
        fitness(i)=fobj(X(i,:),Values);
        if fitness(i)<personalbest_faval(i) %判断当前位置是否是历史上最佳位置
            personalbest_faval(i)=fitness(i);
            personalbest_x(i,:)=X(i,:);
        end
    end
    [globalbest_faval,i]=min(personalbest_faval);
    globalbest_x=personalbest_x(i,:);
    his_best(k) = globalbest_faval;
    for i=1:particlesize %更新粒子群里每个个体的最新位置
        V(i,:)=w*V(i,:)+c1*rand*(personalbest_x(i,:)-X(i,:))...
            +c2*rand*(globalbest_x-X(i,:));
        for j=1:dim    %判断粒子的飞翔速度是否超过了最大飞翔速度
            if V(i,j)>vmax
                V(i,j)=vmax;
            elseif V(i,j)<-vmax
                V(i,j)=-vmax;
            end
        end
        X(i,:)=X(i,:)+V(i,:);
    end
    k=k+1;
end
end



function population = constraint_population(population, Values, Weights, maxw)
    for i = 1:size(population, 1)
        tmp_value = 0;
        tmp_weight = 0;
        for j = 1:size(population,2)
            if population(i, j) > 0.5
                tmp_value = tmp_value+Values(j);
                tmp_weight = tmp_weight+Weights(j);
            end
        end

        while tmp_weight > maxw
            indices_of_ones = find(population(i,:) > 0.5);
            random_index = indices_of_ones(randi(size(indices_of_ones,2)));
            population(i, random_index) = rand * 0.5;
            tmp_value = tmp_value - Values(random_index);
            tmp_weight = tmp_weight - Weights(random_index);
        end
    end

end

function fitness = fobj(X, Values)
% 求解最小值，故为负数
    fitness = 0;
    for i = 1:size(X,2)
        if X(i) > 0.5
            fitness = fitness - Values(i);
        end
    end
end

