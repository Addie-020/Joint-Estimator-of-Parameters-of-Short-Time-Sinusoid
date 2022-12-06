clear
clc
close all

tic
x0 = [0.2; 0.5];
xLb = [-1; -1];
xUb = [1; 1];

% [xBest, yBest, info] = ConjGradeOptim(x0, 0, 0);
[xBest, yBest, info] = ConjGradeCons(x0, xLb, xUb, 0, 0);

toc