clc
clear

debug=1;

%% 

%Ex3 (Discretização do Exemplo 2 de [0]Sophie)
Ac=[0 1 ; 1 0];
Bc=[0;-5];
wn=4;
K1=[0.4741    0.4073];
K2= -0.6060;

%B=nxm
[n,m] = size(Bc);
Im = eye(m);
In = eye(n);

%% Discretize

%Frequence minimal
fn=wn/(2*pi);
Ts=1/(2*fn);

[A,B]=Discretize(Ac,Bc,Ts,debug);

%% Check if is Schur Cohm

ItsSchurCohn(A, B, K1,debug);

%% LMI

tal=0.66;

if ~exist('tal')
    tal_min=0.01;
    tal_max=1;
    pas=0.01;
else
    tal_min=tal;
    tal_max=tal;
    pas=1;
end

for tal=tal_min:pas:tal_max
    pause(0.01)
    
    [feas,sol]=Analysis_Regional(A,B,K1,K2,tal, [1e-20, 1e-20, 1e-20],debug);
    
    %% Plot traject
    
    if round(feas,3) < 0
        x0=[-5; 5]; N=5;
        PlotRegional(A, B, Ts, sol, x0, N, tal, [35,35]);
    end

end