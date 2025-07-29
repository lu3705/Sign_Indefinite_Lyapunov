clear
clc

debug = 1;

%% System

%System continu
Ac = [-0.05 1; -10 -0.5];
Bc = [0; 1];
wn = 4;

%% Discretize

%Frequence minimal
fn = wn/(2*pi);
Ts = 1/(2*fn);

[A,B] = Discretize(Ac,Bc,Ts,debug);

%% LMI

r = 0.4; %maximo 1
[feas,sol] = Synthesis_Global(A,B,r,debug, [1e-20, 1e6, 1e-20]);


%% Check if it's Schur Cohn

ItsSchurCohn(A, B, sol.K1,debug);

%% Plot traject

N=10;
x0 = [10;10];
PlotGlobal(Ac, Bc, Ts, sol, x0, N, [35,35]);
