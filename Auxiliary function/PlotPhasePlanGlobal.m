
clear
clc

debug =0;

%% System

%motor 
Ra=1; La=0.5; J=0.01;
bm=0.1; Km=0.01; Kb=0.01;

Ac = [-bm/J  -Km/J; Kb/La -Ra/La];
Bc = [0;1/La];

wn=4;

%% Discretize

%Frequence minimal
fn = wn/(2*pi);
Ts = 1/(2*fn);

[A,B] = Discretize(Ac,Bc,Ts,debug);

%% LMI

r = 0.4; %maximo 1
[feas,sol] = YALMIP_Synthesis_Global(A,B,r,debug, [1e-20, 1e-20, 1e-20]);


%% Check if it's Schur Cohn

ItsSchurCohn(A, B, sol.K1,debug);

%% Plot traject

x0 = [10;10];
range=[20,20];
% PlotGlobal(Ac, Bc, Ts, sol, x0, range);


%%
N=500;

range=[10,10];
range_x=range(1);
range_y=range(2);
    
K1=sol.K1;
K2=sol.K2;
P0=sol.P0;
P=sol.P;
H1=sol.H1;
H2=sol.H2;

[n,~] = size(B);

figure;
hold on;

[x1_grid, x2_grid] = meshgrid(linspace(-range_x, range_x, 24), ...
                              linspace(-range_y, range_y, 24));

for i = 1:size(x1_grid, 1)
    for j = 1:size(x1_grid, 2)
        x = [x1_grid(i,j); x2_grid(i,j); zeros(n-2,1)];

        [~, u_sat, ~] = Saturation(x, K1, K2);
        x_next = A*x + B*u_sat;

        v = x_next(1:2) - x(1:2);
        norm_v = norm(v) + 1e-9;  % evita divisão por zero
        v = v / norm_v;

        scale = 0.6;  % ajuste aqui conforme desejar
        dx = scale * v(1);
        dy = scale * v(2);

        h=quiver(x(1), x(2), dx, dy, 0, ...
            'Color', "#BB4CCF", ...
            'LineWidth', 1, ...
            'MaxHeadSize', 3, ...
            'HandleVisibility', 'off');
    end
end

x=[8;6];
traj=[];
V=0;
for k=1:N/Ts
    [u, u_sat, dz] = Saturation(x, K1, K2);

    Vant=V;
    V = [x; dz]'*P*[x; dz];
    dV=V-Vant;

    data = [x; V; dV; u; u_sat];
    traj = [traj data];

    x = A*x + B*u_sat;
end
        plot(traj(1,:) ,traj(2,:),'-.','LineWidth',1.3,'Color','b'); hold on
        plot(traj(1,1) ,traj(2,1),'xr','LineWidth',1.5,'MarkerSize',10,'Color','#001BBA'); hold on

        
x=[-8;-6];
traj=[];
V=0;
for k=1:N/Ts
    [u, u_sat, dz] = Saturation(x, K1, K2);

    Vant=V;
    V = [x; dz]'*P*[x; dz];
    dV=V-Vant;

    data = [x; V; dV; u; u_sat];
    traj = [traj data];

    x = A*x + B*u_sat;
end
        plot(traj(1,:) ,traj(2,:),'-.','LineWidth',1.6,'Color','#F5890B'); hold on
        plot(traj(1,1) ,traj(2,1),'xr','LineWidth',1.6,'MarkerSize',10,'Color','#DE7C0A'); hold on
        legend({'Trajectory 1','Start 1','Trajectory 2','Start 2'}, 'NumColumns', 2)



xlim([-range_x, range_x])
ylim([-range_y, range_y])
xlabel('$x_1$', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')