clc
clear

debug=1;

%Ex3 (Discretização do Exemplo 2 de [0]Sophie)
Ac=[0 1 ; 1 0];
Bc=[0;-5];
wn=4;

%B=nxmy
[n,m] = size(Bc);
Im = eye(m);
In = eye(n);

%% Discretize

%Frequence minimal
fn=wn/(2*pi);
Ts=1/(2*fn);

[A,B]=Discretize(Ac,Bc,Ts,debug);

%% LMI

r=0.2; %max 1
[feas,sol]=Synthesis_Regional(A,B,r, [1e-20, 1e-20, 1e-20, 1e-20],debug);

eig1=sol.eig1{1};
eig2=sol.eig2;
eig3=sol.eig3;
eig4=sol.eig4;

%% It is Schur Cohm?

ItsSchurCohn(A, B, sol.K1,debug,r);

%% Plot traject
x0=[-6; 9];
N=5;
PlotRegional(A, B, Ts, sol, x0, N, r, [25,25]);