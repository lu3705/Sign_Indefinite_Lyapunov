clc
clear all
dt = 0.01;

A=[-0.05 1 ; -10 -0.5];
B=[0;1];
K1=[9.9 0.495];
K2=0;

% Ad=[0.7 0.1; 0.3 -0.4];
% Bd=[-0.5; 0];
% K1=[1 0.2];
% K2=0;
% 
% A=(Ad-eye(size(Ad)))/dt;
% B=Bd/dt;

%B=nxm
[n,m] = size(B);
Im = eye(m);
In = eye(n);
H2=zeros(m,m);
H1=zeros(m,n);

%%

setlmis([]) %Initialization

type=1;
%type = 1,P square and symmetrical
%     = 2,P rectangular [m,n]
%     = 3,P is of other type
%tipo,dim[l,c]
P0=lmivar(type,  [n,1]);          %quadrado e simétrico
P1=lmivar(type+1,[n,m]);          %retangular
P2=lmivar(type,  [m,1]);          %quadrado e simetrico 
T =lmivar(type,  [m,1]);          %quadrado e simétrico
T1=lmivar(type,  [m,1]);          %quadrado e simétrico
T2=lmivar(type,  [m,1]);          %quadrado e simétrico
T3=lmivar(type,  [m,1]);          %quadrado e simétrico

%Condição 1 (V>0)
% |P0   2[P1-(K1-H1)'T] | >0
% | 0   P2-2T(K2-Im-H2) |    
lmiterm([-1 1 1 P0],1, 1);  

lmiterm([-1 1 2 P1],1, 1);
lmiterm([-1 1 2 T],-1*(K1-H1)',1);

lmiterm([-1 2 2 P2],1,1);  
lmiterm([-1 2 2 T],-2,(K2-Im-H2)); 

%Condição 2 (V.<0)
% |P0*(A+BK1)  (A+B*K1)'*P1+P0*B*(K2-Im)+[K1(A+B*K1)]'*T3+(K1-H1)'*T1      P1+[K1(A+BK1)]'T2           | 
% | 0               T3*K1*B(K2-Im)+T1*(K2-Im-H2)+[B*(K2-Im)]'*P1       T3*(K2-Im)+P2+[K1*B(K2-Im)]'*T2 | <0
% | 0                             0                                           T2*(K2-Im)               |

ABK1=A+B*K1;
K2I=K2-Im;
lmiterm([2 1 1 P0],1, ABK1);  

lmiterm([2 1 2 P1],(ABK1)', 1);
lmiterm([2 1 2 P0],1,B*(K2I));
lmiterm([2 1 2 T3],(K1*(ABK1))',1);
lmiterm([2 1 2 T1],(K1-H1)',1);

lmiterm([2 2 2 P1],(B*(K2I))',1); 
lmiterm([2 2 2 T3],1,K1*B*(K2I));  
lmiterm([2 2 2 T1],1,K2I-H2); 

lmiterm([2 1 3 P1],1,1);  
lmiterm([2 1 3 T2],(K1*(ABK1))',1); 

lmiterm([2 2 3 T3],1,K2I);  
lmiterm([2 2 3 P2],1,1);
lmiterm([2 2 3 T2],(K1*B*(K2I))',1); 

lmiterm([2 3 3 T2],1,(K2I));  


%P0>0
lmiterm([-3 1 1 P0],1,1);

%T>0
lmiterm([-4 1 1 T],1,1);
lmiterm([-5 1 1 T1],1,1);
lmiterm([-6 1 1 T2],1,1);
lmiterm([-7 1 1 T3],1,1);

%declara sistema
LMIsys=getlmis;

options=[1e-5, 0, 0, 0, 0];
[feas,sol]=feasp(LMIsys, options, -0.5);                
%feas<0 (factivel)
%sol conjunto de soluções

%extrai valores da solução
P0opt=dec2mat(LMIsys,sol,P0);
P1opt=dec2mat(LMIsys,sol,P1);
P2opt=dec2mat(LMIsys,sol,P2);
Topt =dec2mat(LMIsys,sol,T );
T1opt=dec2mat(LMIsys,sol,T1);
T2opt=dec2mat(LMIsys,sol,T2);
T3opt=dec2mat(LMIsys,sol,T3);

P=[P0opt P1opt; P1opt' P2opt];

autovalores=eig(P);


%% Plot traject
x =[10;10];
N = 1000;
traj=[];
V_values=[];
dV_values=[];

for k=1:N
    traj=[traj x];

    u=K1*x;
    u_sat=min(1,max(-1,u));
    dz=u-u_sat;

    V=[x; dz]'*P*[x; dz]; 
    V_values=[V_values V];

    if k~=1
        dV=(V_values(k)-V_values(k-1))/dt;
        dV_values=[dV_values dV];
    end

    dx =A*x+B*u_sat;
    x=x+dt*dx;
end
figure
subplot 311
plot(traj(1,:), traj(2,:));
title("trajectory states");
subplot 312
plot(V_values);
title("trajectory V")
subplot 313
plot(dV_values,'.');
title("trajectory dV")

%%
[x1_grid,x2_grid]=meshgrid(linspace(-15,15,30), linspace(-35,30,30));

Vplot = zeros(size(x1_grid));
for i = 1:size(x1_grid, 1)
    for j = 1:size(x1_grid, 2)
        x = [x1_grid(i,j); x2_grid(i,j)];

        u = K1 * x;
        u_sat = min(1, max(-1, u));
        dz = u - u_sat;

        Vplot(i,j) = [x; dz]' * P * [x; dz];
    end
end

figure
plot3(traj(1,:), traj(2,:),V_values, '.r')
legend("trajectoire V","V")
hold on
mesh(x1_grid, x2_grid, Vplot)
shading interp; colorbar;

figure
plot(1:N,traj(1,:))
