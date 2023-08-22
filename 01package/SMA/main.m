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

%% 定义超参数
N=50;   % Number of particles
Max_iteration=1000; % Maximum number of iterations
lb = 0;
ub = 1;
dim = maxn;

[Xprey, Food_Score, CNVG]=SMA(Max_iteration,N,lb,ub,dim,Values,Weights,maxw);

Xprey

fitness = 0;
weight = 0;
for i = 1:size(Food_Score,2)
    if Food_Score(i) > 0.5
        fitness = fitness + Values(i);
        weight = weight + Weights(i);
    end
end
fitness
weight
plot(-1*CNVG)