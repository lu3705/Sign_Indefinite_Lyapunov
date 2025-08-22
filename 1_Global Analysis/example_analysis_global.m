clear
clc

debug = 1;

%% System

%motor 
Ra=1; La=0.5; J=0.01;
bm=0.1; Km=0.01; Kb=0.01;

Ac = [-bm/J  -Km/J; Kb/La -Ra/La];
Bc = [0;1/La];

wn=4;

K1= [-0.0009   -0.2294];
K2=-0.2246;

%% Discretize

%Frequence minimal
fn = wn/(2*pi);
Ts = 1/(2*fn);          %%max value of Ts

[A,B] = Discretize(Ac,Bc,Ts,debug);

%% Check if it is Schur Cohn

ItsSchurCohn(A, B, K1,debug);

%% LMI

[feas,sol] = YALMIP_Analysis_Global(A, B, K1, K2, debug);

%% Plot traject

x0 = [40;-40];
PlotGlobal(Ac, Bc, Ts, sol, x0, [50,50]);
