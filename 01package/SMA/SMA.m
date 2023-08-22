function [Destination_fitness,bestPositions,Convergence_curve]=SMA(Max_iter,N,lb,ub,dim,Values,Weights,maxw)
    % initialize position
    bestPositions=zeros(1,dim);
    Destination_fitness=inf;%change this to -inf for maximization problems
    AllFitness = inf*ones(N,1);%record the fitness of all slime mold
    weight = ones(N,dim);%fitness weight of each slime mold

    %随机初始化种群，记住这里进行初始化的时候
    X=initialization(N,dim,ub,lb);

    Convergence_curve=zeros(1,Max_iter);
    it=1;  %Number of iterations
    lb=ones(1,dim).*lb; % lower boundary 
    ub=ones(1,dim).*ub; % upper boundary
    z=0.03; % parameter
    
    % Main loop
    while  it <= Max_iter
        
        %sort the fitness
        for i=1:N
            % Check if solutions go outside the search space and bring them back
            Flag4ub=X(i,:)>ub;
            Flag4lb=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;    %这个公式很常见
            X(i,:) = constraint_population(X(i,:), Values, Weights, maxw);
            AllFitness(i) = fobj(X(i,:),Values);
        end
        
        [SmellOrder,SmellIndex] = sort(AllFitness);  %Eq.(2.6)
        worstFitness = SmellOrder(N);
        bestFitness = SmellOrder(1);
    
        S=bestFitness-worstFitness+eps;  % plus eps to avoid denominator zero
        %eps  matlab中最小的正数
        %calculate the fitness weight of each slime mold
        for i=1:N
            for j=1:dim
                if i<=(N/2)  %Eq.(2.5)
                    weight(SmellIndex(i),j) = 1+rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
                else
                    weight(SmellIndex(i),j) = 1-rand()*log10((bestFitness-SmellOrder(i))/(S)+1);
                end
            end
        end
        
        %update the best fitness value and best position
        if bestFitness < Destination_fitness
            bestPositions=X(SmellIndex(1),:);
            Destination_fitness = bestFitness;
        end
        
        a = atanh(-(it/Max_iter)+1);   %Eq.(2.4)
        b = 1-it/Max_iter;     %这两个很类似MRFO和EO的参数，越到后面概率越大
        % Update the Position of search agents
        for i=1:N
            if rand<z     %Eq.(2.7)
                X(i,:) = (ub-lb)*rand+lb;   %变异，这里有问题，概率是死的，容易局部最优
            else
                p =tanh(abs(AllFitness(i)-Destination_fitness));  %Eq.(2.2)
                vb = -a + 2*a*rand(1,dim); %Eq.(2.3)
                vc = -b + 2*b*rand(1,dim);
                for j=1:dim
                    r = rand();
                    A = randi([1,N]);  % two positions randomly selected from population
                    B = randi([1,N]);
                    if r<p    %Eq.(2.1)
                        X(i,j) = bestPositions(j)+ vb(j)*(weight(i,j)*X(A,j)-X(B,j));
                    else
                        X(i,j) = vc(j)*X(i,j);
                    end
                end
            end
        end
        Convergence_curve(it)=Destination_fitness;
        it=it+1;
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




