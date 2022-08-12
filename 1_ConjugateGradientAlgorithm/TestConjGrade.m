clear
clc
close all

tic
x0 = [0.2; 0.5; 0];

[xBest, yBest, info, dataLog] = ConjGradeOptim(x0, 0, 0);

toc