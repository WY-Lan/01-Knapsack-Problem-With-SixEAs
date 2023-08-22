clc;
clear all;
close all;
%% 获取文件
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


MaxNum=1000;
dim=maxn;
particlesize=50;
lb = 0;
ub = 1;

[his_best, best_pop] = ACO(MaxNum,particlesize,lb,ub,dim,Values,Weights,maxw);

fitness = 0;
weight = 0;
for i = 1:size(best_pop,1)
    if best_pop(i) > 0.5
        fitness = fitness + Values(i);
        weight = weight + Weights(i);
    end
end
fitness
weight
plot(-1*his_best);



