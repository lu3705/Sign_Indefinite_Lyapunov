clear
clc

debug = 1;

%% System

Ra=1; La=0.5; J=0.01;
bm=0.1; Km=0.01; Kb=0.01;

Ac = [-bm/J  -Km/J; Kb/La -Ra/La];
Bc = [0;1/La];

wBW=11;

%% Discretize

%Frequence minimal
fBW = wBW/(2*pi);
Ts = 1/(2*fBW);  %Ts precisa ser menor do que isso


[A,B] = Discretize(Ac,Bc,Ts,debug);

%% LMI

r = 0.1; %maximo 1
[feas,sol] = YALMIP_Synthesis_Global(A,B,r,debug, [1e-4, 1e-4, 1e-4]);


%% Check if it's Schur Cohn

ItsSchurCohn(A, B, sol.K1,debug,r);

%% Plot traject

x0 = [30;-30];
PlotGlobal(Ac, Bc, Ts, sol, x0, [50,50]);
