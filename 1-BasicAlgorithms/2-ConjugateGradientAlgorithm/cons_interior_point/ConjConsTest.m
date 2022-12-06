clear
clc
close all

tic
x0 = [10; 10];
xLb = [-100; -100];
xUb = [100; 100];

% [xBest, yBest, info] = ConjGradeOptim(x0, 0, 0);
[xBest, yBest, info] = ConjGradeCons(x0, xLb, xUb, 0, 0);

toc