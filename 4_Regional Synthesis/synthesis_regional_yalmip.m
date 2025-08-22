clc
clear

debug=0;

% %Ex3 (Example 2 - [1])
Ac=[0 1 ; 1 0];
Bc=[0;-5];
wn=4;

% Ac=[0 1 ; 0 1];
% Bc=[-1;1];
% wn=8;

% Pendulum inverse
% g=9.81; l=0.5; m=0.15;
% mi=0.05; talb=5; dt=0.02;
% A=[1 dt 0; g*dt/l 1-(mi*dt/(m*l^2)) 0; 1 0 1];
% B=[0; talb*dt/(m*l^2);0];
% Ts=1;

%% Discretize

%Frequence minimal
fn=wn/(2*pi);
Ts=1/(2*fn);

[A,B]=Discretize(Ac,Bc,Ts,debug);

%% LMI

r=0.2; %max 1
[feas,sol]=YALMIP_Synthesis_Regional(A, B, r, debug, [1e-4, 1e-4, 5e-5, 5e-20]);

eig1=sol.eig1{1};
eig2=sol.eig2;
eig3=sol.eig3;
eig4=sol.eig4;

%% It is Schur Cohm?

ItsSchurCohn(A, B, sol.K1,debug,r);

%% Plot traject
x0=[-4; 0.8];
PlotRegional(A, B, Ts, sol, x0, r, [15,15,4],0);
