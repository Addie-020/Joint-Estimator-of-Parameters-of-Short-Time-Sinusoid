clear
clc
close all

tic
objFun = @TestFun;
x0 = [0; 0];
xLb = [-1; -1];
xUb = [1; 1];
rng default                 % For reproducibility
nvars = 2;                  % Number of variables
xFin = particleswarm(objFun, nvars, xLb, xUb);
% [xBest, fBset, info, dataLog] = ParticleSwarmOptim(objFun, x0, xLb, xUb);
toc