clc;clear all;close all;
tic;                              %程序运行计时

%% 读取数据
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

%% 设定初始化的参数
MaxNum=1000;                    %粒子最大迭代次数
dim=maxn;                         %目标函数的自变量个数
particlesize=50;                    %粒子群规模
vmax=0.2;                        %粒子的最大飞翔速度
lb = 0;
ub = 1;

[globalbest_x,globalbest_faval,his_best] = PSO(MaxNum,particlesize,lb,ub,dim,Values,Weights,maxw);

plot(-1*his_best);
fitness = 0;
weight = 0;
for i = 1:size(globalbest_x,2)
    if globalbest_x(i) > 0.5
        fitness = fitness + Values(i);
        weight = weight + Weights(i);
    end
end
fitness
weight

toc;