clc
clear

debug = 0;

%% System

%Example 2 de [1]
Ac = [0 1 ; 1 0];
Bc = [0;-5];
wn = 4;

r  = 0.95
K1 = [-0.0186   -1.3995];
K2 = 0.9963;

%% Discretize

%Frequence minimal
fn = wn/(2*pi);
Ts = 1/(2*fn);

[A,B]=Discretize(Ac,Bc,Ts,debug);

%% Check if it's Schur Cohm

ItsSchurCohn(A, B, K1,debug);

%% LMI

% tau = 0.12;

if ~exist('tau')
    tau_min = 0.01;
    tau_max = 0.5;
    pas = 0.01;
    debug_plot = 0;
else
    tau_min = tau;
    tau_max = tau;
    pas = 1;
    debug_plot = debug;
end

for tau=tau_min:pas:tau_max
    pause(0.01)
    
    [feas,sol] = YALMIP_Analysis_Regional(A,B,K1,K2,tau,debug, [1e-3, 1e-3,1e-20]);
    
    %% Plot traject
   
    if round(feas) == 0
        x0 = [-6; 10];
        PlotRegional(A, B, Ts, sol, x0, tau, [15,15],debug_plot);
    end

end
