function [Xprey, Food_Score,CNVG] = IHBA(tmax,N, lb, ub, dim, Values, Weights, maxw)
    beta= 6;     % the ability of HB to get the food  Eq.(2)
    w = 0.8;     %constant in Eq. (8)
    vec_flag=[1,-1];
    CNVG=zeros(tmax,1);
    %initialization
    X=initialization(N,dim,ub,lb);
    X = constraint_population(X, Values, Weights, maxw);
    %Evaluation
    Xnew=X;
    for i=1:N
        fitness(i) = fobj(X(i,:),Values);
    end

    [GYbest, gbest] = min(fitness);
    Xprey = X(gbest,:);

    for t = 1:tmax
        a=(cos(2*rand)+1)*(1-t/tmax);   %Eq. (5)
        alpha=a*(2*rand-1); 
        C=alpha;
        I=Intensity(N,Xprey,X); %intensity in Eq. (1)
        fAvg=mean(fitness);
        z=rand;
        for j=1:dim     %Eq. (4)
            if z<0.5
                z =2*z+rand/N;
                w1(j)=z;
            else
                z =2*(1-z)+rand/N; 
                w1(j)=z;
            end
        end
        %探索
        for i=1:N    %Eq. (7)
            F=vec_flag(floor(2*rand()+1));        
                di=((Xprey-X(i,:)));
                if rand<.5
                    Z = -a + 2*a*rand(1,dim);   %生成(连续)均匀分布的随机数   其中a是约束
                    H=Z.*X(i,:);   %定义是放大还是缩小
                    r3=rand;                
                    r4=rand;                
                    r5=rand;                
                    Xnew(i,:)=C.*H+(rand-a)*X(randi([1,N]),:)+F*beta*I(i)* Xprey+F*r3*alpha*(di)*abs(cos(2*pi*r4)*(1-cos(2*pi*r5)));
                else
                    r7=rand;
                    Xnew(i,:)=X(i,:)+C.*(C*(X(i,:)- Xnew(randi([1,N]),:))+rand*(X(i,:)-Xnew(randi([1,N]),:)))+F*r7*alpha*di;
                end 
                    
        end
        for i=1:N
            if fitness( i ) < fAvg
                Xnew1( i, : )=(randn+1).*X( i, : );
            else
                Xnew1( i, : )= repmat(lb,1,dim)  + (ub - lb)*w1;
                Xnew1( i, : )=(Xnew1( i, : )+X( i, : ))/2;
            end
            FU=Xnew1(i,:)>ub;FL=Xnew1(i,:)<lb;Xnew1(i,:)=(Xnew1(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;  
            
            Xnew1(i,:) = constraint_population(Xnew1(i,:), Values, Weights, maxw);
            tempFitness1 =fobj(Xnew1(i,:),Values); 
            FU=Xnew(i,:)>ub;FL=Xnew(i,:)<lb;Xnew(i,:)=(Xnew(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
            
            Xnew(i,:) = constraint_population(Xnew(i,:), Values, Weights, maxw);
            tempFitness = fobj(Xnew(i,:),Values); 
            if min([tempFitness,tempFitness1])<fitness(i)
                if tempFitness<fitness(i)
                    fitness(i)=tempFitness;
                    X(i,:)= Xnew(i,:);
                else tempFitness1 < fitness(i);            
                    fitness(i)=tempFitness1;
                    X(i,:)= Xnew1(i,:);
                end
            end
        end
        if min(fitness)<GYbest
            [GYbest,loc]=min(fitness);   
            GYbest=GYbest;
            Xprey=X(loc,:);
        end  
    
        %开发
        for i=1:N           %Eq. (8)
            if abs(a)>=w  
                delta= mean(Xnew,1);
                Xnew(i,:)=C*delta.*(X(i,:)-Xprey)+X(i,:);  %探索
            else           
               if rand>=0.5
                  h=randn(1,dim);%理论上最大值是正无穷，它产生的数均值是0，标准差是1，但实际上你不会看到很大的数产生（概率很小，几乎为0） 
               else
                  h=randn(1,1);%理论上最大值是正无穷，它产生的数均值是0，标准差是1，但实际上你不会看到很大的数产生（概率很小，几乎为0） 
                  % randn(m,n) 或 Y = randn([m n])     randn（random normal distribution）是一种产生标准正态分布的随机数或矩阵的函数
               end
               r1=rand; 
               Xnew(i,:)= Xprey-(Xprey*(2*r1-1)-X(i,:)*(2*r1-1)).*(h);  %2*r1-1=(-1,1)    Beta*h加强正态分布           
            end
          
        end
        for i=1:N
            FU=Xnew(i,:)>ub;FL=Xnew(i,:)<lb;Xnew(i,:)=(Xnew(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
            
            Xnew(i,:) = constraint_population(Xnew(i,:), Values, Weights, maxw);
            tempFitness = fobj(Xnew(i,:),Values); 
            if tempFitness<fitness(i)
                fitness(i)=tempFitness;
                X(i,:)= Xnew(i,:);
            end 
            if tempFitness<GYbest
                GYbest=tempFitness;
                Xprey=Xnew(i,:);
            end
        end
    
        %混沌    
        for i = 1 : N         
           if tempFitness<fitness(i)
                fitness(i)=tempFitness;
                X(i,:)= Xnew(i,:);
            end 
            if tempFitness<GYbest
                GYbest=tempFitness;
                Xprey=Xnew(i,:);
            end
        end
    
        [~, SortOrder]=sort(fitness);
        fitness=fitness(SortOrder);
        PopPos_=X;
        for i=1:N
            X(i,:)=PopPos_(SortOrder(i));
        end
    
    
    %     FU=X>ub;FL=X<lb;X=(X.*(~(FU+FL)))+ub.*FU+lb.*FL;
    %     [Ybest,index] = min(fitness);
        CNVG(t)=GYbest;
    %     if Ybest<GYbest
    %         GYbest=Ybest;
    %         Xprey = X(index,:);
    %     end    
    end
    Food_Score = GYbest;
end


function I=Intensity(N,Xprey,X)
    for i=1:N-1
        di(i) =( norm((X(i,:)-Xprey+eps))).^2;
        S(i)=( norm((X(i,:)-X(i+1,:)+eps))).^2;
    end
    di(N)=( norm((X(N,:)-Xprey+eps))).^2;
    S(N)=( norm((X(N,:)-X(1,:)+eps))).^2;
    for i=1:N
        r2=rand;
        I(i)=r2*S(i)/(4*pi*di(i));
    end
end

function [X]=initialization(N,dim,up,down)
    if size(up,2)==1
        X=rand(N,dim).*(up-down)+down;
    end
    if size(up,2)>1
        for i=1:dim
            high=up(i);low=down(i);
            X(:,i)=rand(N,1).*(high-low)+low;
        end
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
