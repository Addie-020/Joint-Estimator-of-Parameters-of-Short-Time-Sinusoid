clear
clc
close all

tic
x0 = [0.2; 0.5];

[xBest, yBest, info] = ConjGradeOptim(x0, 0, 0);

toc