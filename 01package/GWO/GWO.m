
function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,lb,ub,dim,Values,Weights,maxw)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
Positions = constraint_population(Positions, Values, Weights, maxw);

Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter
fitness=repmat(inf,size(Positions,1) ,1); 
% Main loop
while l<Max_iter
    Positions = constraint_population(Positions, Values, Weights, maxw);
    for i=1:size(Positions,1)  
       % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
        % Calculate objective function for each search agent
        fitness(i)=fobj(Positions(i,:),Values);
        
        % Update Alpha, Beta, and Delta
        if fitness(i)<Alpha_score 
            Alpha_score=fitness(i); % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness(i)>Alpha_score && fitness(i)<Beta_score 
            Beta_score=fitness(i); % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness(i)>Alpha_score && fitness(i)>Beta_score && fitness(i)<Delta_score 
            Delta_score=fitness(i); % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    new_Positions=zeros(1,dim);
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            new_Positions(j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
        new_fitness=fobj(new_Positions, Values);
        old_fitness=fobj(Positions(i,:), Values);
        if new_fitness<old_fitness
            Positions(i,:)=new_Positions;
        end
    end
    l=l+1; 
%     [~, SortOrder]=sort(fitness);
%     fitness=fitness(SortOrder);
%     PopPos_=Positions;
%     for i=1:size(Positions,1)
%         Positions(i,:)=PopPos_(SortOrder(i));
%     end

    Convergence_curve(l)=Alpha_score;
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
