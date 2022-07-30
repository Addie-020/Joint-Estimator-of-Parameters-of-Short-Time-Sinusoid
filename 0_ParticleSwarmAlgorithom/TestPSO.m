clear
clc
close all

tic
objFun = @TestFun;
x0 = [0.2; 0.5];
xLb = [-2; -2];
xUb = [2; 2];
% rng default                 % For reproducibility
% nvars = 2;                  % Number of variables
% xFin = particleswarm(objFun, nvars, xLb, xUb);
[xBest, fBset, info, dataLog] = ParticleSwarmOptim(objFun, x0, xLb, xUb);
toc